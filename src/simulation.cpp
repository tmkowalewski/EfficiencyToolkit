#include "simulation.hpp"

#include <ROOT/TTreeProcessorMT.hxx>
#include <TBranch.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

constexpr int kAddbackThreshold = 150; // Energy threshold for add-back histograms (in keV)

////////////////////////////////////////////////////////////////////////////////
/// Return monotonic time in nanoseconds for progress estimates.
static long long GetSteadyTimeNs()
{ return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

////////////////////////////////////////////////////////////////////////////////
/// Format an estimated duration for the progress display.
static std::string FormatEstimatedTime(double seconds)
{
    if (!std::isfinite(seconds) || seconds < 0) return "--:--";
    long long s = static_cast<long long>(seconds + 0.5);
    char      buf[32];
    if (s >= 3600)
        std::snprintf(buf, sizeof buf, "%lld:%02lld:%02lld", s / 3600, (s / 60) % 60, s % 60);
    else
        std::snprintf(buf, sizeof buf, "%02lld:%02lld", s / 60, s % 60);
    return buf;
}

////////////////////////////////////////////////////////////////////////////////
/// Format one progress bar with a percentage and estimated remaining time.
static std::string FormatProgressBar(const std::string& label, Long64_t done, Long64_t total, long long startNs)
{
    const int barWidth = 40;
    double    frac     = total > 0 ? std::min(1.0, static_cast<double>(done) / total) : 0.0;
    int       filled   = static_cast<int>(frac * barWidth);
    double    elapsed  = (GetSteadyTimeNs() - startNs) * 1e-9;
    double    eta      = (done > 0 && startNs > 0) ? elapsed * (1.0 - frac) / frac : std::nan("");

    std::ostringstream os;
    os << std::left << std::setw(22) << label << " [" << std::string(filled, '=') << std::string(barWidth - filled, ' ') << "] " << std::right << std::fixed << std::setprecision(1) << std::setw(5) << frac * 100.0 << "%  ETA " << FormatEstimatedTime(eta);
    return os.str();
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw the current-file and total progress bars until sorting stops.
static void DisplaySortProgress(SortProgress_t& p)
{
    bool first = true;
    auto draw  = [&]() {
        std::string name;
        {
            std::lock_guard<std::mutex> lock(p.fNameMutex);
            name = p.fFileName + " (" + std::to_string(p.fFileIndex) + "/" + std::to_string(p.fFileCount) + ")";
        }
        if (!first) std::cout << "\033[2F";
        first = false;
        std::cout << "\033[2K" << FormatProgressBar(name, p.fFileDone, p.fFileTotal, p.fFileStartNs) << "\n"
                  << "\033[2K" << FormatProgressBar("Total", p.fAllDone, p.fAllTotal, p.fAllStartNs) << std::endl;
    };
    while (!p.fStop)
    {
        draw();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }
    draw();
}

////////////////////////////////////////////////////////////////////////////////
/// Sort one tree into detector histograms, updating optional progress counters.
void SortSimulationTree(const char* fileName, const char* treeName, const char* outFileName, SortProgress_t* progress)
{
    ROOT::EnableImplicitMT();
    TH1::AddDirectory(false); // histograms are owned/written explicitly; avoids thread-unsafe directory registration

    // Probe the tree for its branch layout
    std::unique_ptr<TFile> probeFile(TFile::Open(fileName, "READ"));
    TTree*                 tree = probeFile ? probeFile->Get<TTree>(treeName) : nullptr;
    if (!probeFile || probeFile->IsZombie() || !tree) { throw std::runtime_error(std::string("Could not read tree '") + treeName + "' from " + fileName); }
    const Long64_t numEntries = tree->GetEntries();

    std::map<std::string, std::vector<std::string>> detectorXtals; // Map detector name to list of crystal branches
    std::map<std::string, std::vector<std::string>> allDetectors;  // Map of all detectors (single and multi)
    std::map<std::string, TH1D*>                    xtalHists;
    std::map<std::string, TH1D*>                    addbackHists;
    std::vector<std::string>                        branchNames;

    // First pass: categorize branches and create histograms
    for (auto branch : *tree->GetListOfBranches())
    {
        std::string branchName = branch->GetName();

        // Skip event ID branch
        if (branchName == "evid") continue;

        // Determine if this is a multi-crystal detector and extract detector name
        std::string detectorName = branchName;

        std::smatch match;
        if (std::regex_search(branchName, match, std::regex("(.+)E\\d+$")))
        {
            // Multi-crystal detector: extract everything before E\d+
            detectorName = match[1].str();
            // Track which crystals belong to this detector
            detectorXtals[detectorName].push_back(branchName);
        }

        // Track all detectors
        allDetectors[detectorName].push_back(branchName);

        // Create crystal histogram for all branches (single and multi-crystal)
        std::string xtalHistName = branchName;
        TH1D*       xtalHist     = new TH1D(xtalHistName.c_str(), branchName.c_str(), 10000, 0, 10000); // Adjust range as needed
        xtalHists[branchName]    = xtalHist;

        branchNames.push_back(branchName);
    }

    // Create add-back histograms and sum histogram components for multi-crystal detectors
    for (auto& detPair : detectorXtals)
    {
        std::string               detectorName = detPair.first;
        std::vector<std::string>& xtals        = detPair.second;

        if (xtals.size() > 1) // Only if multiple crystals
        {
            // Create add-back histogram
            std::string addbackHistName = detectorName + "_addback";
            TH1D*       addbackHist     = new TH1D(addbackHistName.c_str(), ("Add-back: " + detectorName).c_str(), 10000, 0,
                                                   10000); // Adjust range as needed
            addbackHists[detectorName]  = addbackHist;
        }
    }
    probeFile->Close();

    // Second pass: fill histograms in parallel; each task fills private histograms and merges them under a lock
    if (progress)
    {
        progress->fFileDone    = 0;
        progress->fFileTotal   = numEntries;
        progress->fFileStartNs = GetSteadyTimeNs();
    }

    std::mutex mergeMutex;
    auto       makeHist = [](const std::string& name, const std::string& title) { return new TH1D(name.c_str(), title.c_str(), 10000, 0, 10000); };

    auto process = [&](TTreeReader& reader) {
        std::vector<std::unique_ptr<TTreeReaderValue<Double_t>>> values;
        std::map<std::string, size_t>                            index;
        std::vector<TH1D*>                                       localXtal;
        for (const auto& name : branchNames)
        {
            index[name] = values.size();
            values.emplace_back(new TTreeReaderValue<Double_t>(reader, name.c_str()));
            localXtal.push_back(makeHist(name, name));
        }
        std::map<std::string, TH1D*> localAddback;
        for (const auto& detPair : addbackHists) localAddback[detPair.first] = makeHist(detPair.second->GetName(), detPair.second->GetTitle());

        Long64_t pending = 0;
        auto     flush   = [&]() {
            if (!progress || pending == 0) return;
            progress->fFileDone += pending;
            progress->fAllDone += pending;
            pending = 0;
        };

        while (reader.Next())
        {
            for (size_t i = 0; i < values.size(); ++i) localXtal[i]->Fill(**values[i]);

            // Add-back histograms (event-by-event summation)
            for (auto& ab : localAddback)
            {
                Double_t addbackEnergy = 0;
                for (const auto& xtalName : detectorXtals[ab.first])
                {
                    Double_t v = **values[index[xtalName]];
                    if (std::isfinite(v) && v > kAddbackThreshold) addbackEnergy += v;
                }
                if (addbackEnergy > 0) ab.second->Fill(addbackEnergy); // Only fill if there's a valid energy
            }

            if (++pending >= 1000) flush();
        }
        flush();

        std::lock_guard<std::mutex> lock(mergeMutex);
        for (size_t i = 0; i < branchNames.size(); ++i) xtalHists[branchNames[i]]->Add(localXtal[i]);
        for (auto& ab : localAddback) addbackHists[ab.first]->Add(ab.second);
        for (auto h : localXtal) delete h;
        for (auto& ab : localAddback) delete ab.second;
    };

    ROOT::TTreeProcessorMT processor(fileName, treeName);
    processor.Process(process);

    if (progress) progress->fFileDone = numEntries;

    // Create an output ROOT file to store the results
    TFile* outFile = TFile::Open(outFileName, "RECREATE");
    if (!outFile || outFile->IsZombie()) throw std::runtime_error(std::string("Cannot create histogram file: ") + outFileName);

    // Create directories for each detector and write histograms
    for (auto& detPair : allDetectors)
    {
        std::string               detectorName = detPair.first;
        std::vector<std::string>& xtals        = detPair.second;

        // Create directory for this detector
        TDirectory* detDir = outFile->mkdir(detectorName.c_str());
        detDir->cd();

        // Write crystal histograms for this detector
        for (const auto& xtalName : xtals) { xtalHists[xtalName]->Write(); }

        // Write add-back histogram if this detector has multiple crystals
        if (detectorXtals.count(detectorName) && detectorXtals[detectorName].size() > 1)
        {
            addbackHists[detectorName]->Write();

            // Create and write sum histogram
            std::string sumHistName = detectorName + "_sum";
            TH1D*       sumHist     = new TH1D(sumHistName.c_str(), ("Sum: " + detectorName).c_str(), 10000, 0,
                                               10000); // Adjust range as needed

            // Add all crystal histograms together
            for (const auto& xtalName : xtals) { sumHist->Add(xtalHists[xtalName]); }

            sumHist->Write();
            delete sumHist;
        }
    }

    outFile->Write();
    if (outFile->TestBit(TFile::kWriteError)) throw std::runtime_error(std::string("Error writing histogram file: ") + outFileName);
    outFile->Close();
    if (outFile->TestBit(TFile::kWriteError)) throw std::runtime_error(std::string("Error closing histogram file: ") + outFileName);
    delete outFile;

    for (auto& h : xtalHists) delete h.second;
    for (auto& h : addbackHists) delete h.second;
}

////////////////////////////////////////////////////////////////////////////////
/// Sort all matching trees while displaying per-file and total progress.
void SortSimulationTrees(const char* directoryName, const char* outDirectoryName, const char* treeName)
{
    std::regex gammaTreePattern("^gamma\\.(\\d{5})\\.root$");

    std::filesystem::create_directories(outDirectoryName);

    std::vector<SortJob_t> jobs;

    for (const auto& entry : std::filesystem::directory_iterator(directoryName))
    {
        if (!entry.is_regular_file()) continue;

        std::string filename = entry.path().filename().string();
        std::smatch match;
        if (!std::regex_match(filename, match, gammaTreePattern)) continue;

        std::unique_ptr<TFile> file(TFile::Open(entry.path().c_str(), "READ"));
        if (!file || file->IsZombie()) { throw std::runtime_error("Failed to open file: " + entry.path().string()); }
        TTree* tree = file->Get<TTree>(treeName);
        if (!tree) { throw std::runtime_error(std::string("Tree '") + treeName + "' not found in " + entry.path().string()); }

        std::string outPath = (std::filesystem::path(outDirectoryName) / ("gamma." + match[1].str() + ".hists.root")).string();
        jobs.push_back({entry.path().string(), outPath, filename, tree->GetEntries()});
    }

    if (jobs.empty()) { throw std::runtime_error(std::string("No gamma.ENERGY.root files found in ") + directoryName); }
    std::sort(jobs.begin(), jobs.end(), [](const SortJob_t& a, const SortJob_t& b) { return a.fName < b.fName; });

    SortProgress_t progress;
    progress.fFileCount = static_cast<int>(jobs.size());
    for (const auto& job : jobs) progress.fAllTotal += job.fEntries;
    progress.fAllStartNs = GetSteadyTimeNs();

    std::thread barThread(DisplaySortProgress, std::ref(progress));

    try
    {
        for (size_t i = 0; i < jobs.size(); ++i)
        {
            {
                std::lock_guard<std::mutex> lock(progress.fNameMutex);
                progress.fFileName  = jobs[i].fName;
                progress.fFileIndex = static_cast<int>(i) + 1;
            }
            SortSimulationTree(jobs[i].fInputPath.c_str(), treeName, jobs[i].fOutputPath.c_str(), &progress);
        }
    } catch (...)
    {
        progress.fStop = true;
        barThread.join();
        throw;
    }

    progress.fStop = true;
    barThread.join();
}

// Statistical (sqrt(N)) and 5% systematic simulation uncertainty, added in quadrature
static const Double_t kSimulatedGammas = 1e7;

////////////////////////////////////////////////////////////////////////////////
/// Combine counting statistics and the simulation's 5% systematic error.
static Double_t CalculateSimulationCountError(Double_t counts)
{
    const Double_t systematic = 0.05 * counts;
    return std::sqrt(counts + systematic * systematic);
}

////////////////////////////////////////////////////////////////////////////////
/// Extract full-energy peaks and write simulation efficiency graphs and CSV data.
void CreateSimulationEfficiencyCurves(const char* directoryName, const char* outFileName)
{
    std::regex                                                        gammaFilePattern("gamma\\.([0-9.]+).hists.root$");
    std::map<std::string, std::vector<std::pair<Double_t, Double_t>>> histogramData; // hist_name -> (energy, bin_content)
    std::map<Double_t, std::string>                                   energyFileMap; // energy -> file path

    std::cout << "Scanning directory: " << directoryName << std::endl;

    // First pass: collect all gamma.ENERGY.root files
    for (const auto& entry : std::filesystem::directory_iterator(directoryName))
    {
        if (!entry.is_regular_file()) continue;

        std::string filename = entry.path().filename().string();
        std::smatch match;

        if (std::regex_search(filename, match, gammaFilePattern))
        {
            Double_t energy       = std::stod(match[1].str());
            energyFileMap[energy] = entry.path().string();
            std::cout << "Found: " << filename << " (Energy: " << energy << " keV)" << std::endl;
        }
    }

    if (energyFileMap.empty()) { throw std::runtime_error(std::string("No gamma.ENERGY.hists.root files found in ") + directoryName); }

    // Second pass: read each file and extract bin contents
    for (const auto& energyPair : energyFileMap)
    {
        Double_t    energy   = energyPair.first;
        std::string filePath = energyPair.second;

        TFile* file = TFile::Open(filePath.c_str(), "READ");
        if (!file || file->IsZombie())
        {
            delete file;
            throw std::runtime_error("Failed to open file: " + filePath);
        }

        // Recursively iterate through all histograms in the file
        std::function<void(TDirectory*)> processDirectory = [&](TDirectory* dir) {
            for (auto key : *dir->GetListOfKeys())
            {
                TObject* obj = dir->Get(key->GetName());

                if (obj->InheritsFrom(TDirectory::Class()))
                {
                    // Recursively process subdirectories
                    processDirectory((TDirectory*)obj);
                } else if (obj->InheritsFrom(TH1D::Class()))
                {
                    TH1D*       hist     = (TH1D*)obj;
                    std::string histName = hist->GetName();

                    // Sum the bin containing this energy and its two neighbors
                    Int_t    bin        = hist->FindBin(energy);
                    Double_t binContent = hist->Integral(std::max(bin - 1, 1), std::min(bin + 1, hist->GetNbinsX()));

                    histogramData[histName].push_back({energy, binContent});

                    std::cout << "  " << histName << " @ " << energy << " keV: " << binContent << std::endl;
                }
            }
        };

        processDirectory(file);
        file->Close();
        delete file;
    }

    // Third pass: create efficiency curves
    std::cout << "\nCreating efficiency curves..." << std::endl;

    auto outFile = TFile::Open(outFileName, "RECREATE"); // Create output ROOT file to store graphs
    if (!outFile || outFile->IsZombie()) throw std::runtime_error(std::string("Cannot create efficiency file: ") + outFileName);

    for (const auto& histPair : histogramData)
    {
        std::string histName   = histPair.first;
        const auto& dataPoints = histPair.second;

        if (dataPoints.empty()) continue;

        // Sort by energy
        std::vector<std::pair<Double_t, Double_t>> sortedData = dataPoints;
        std::sort(sortedData.begin(), sortedData.end());

        // Create graph
        TGraphErrors* graph = new TGraphErrors(sortedData.size());
        graph->SetTitle((histName + " Efficiency").c_str());
        graph->SetName((histName + "_eff").c_str());
        graph->GetXaxis()->SetTitle("Energy (keV)");
        graph->GetYaxis()->SetTitle("Efficiency (FED/Simulated)");

        for (size_t i = 0; i < sortedData.size(); ++i)
        {
            graph->SetPoint(i, sortedData[i].first, sortedData[i].second / kSimulatedGammas);
            graph->SetPointError(i, 0, CalculateSimulationCountError(sortedData[i].second) / kSimulatedGammas);
        }

        // Save TGraph to a ROOT file
        graph->Write();
        delete graph;
    }

    outFile->Write();
    if (outFile->TestBit(TFile::kWriteError)) throw std::runtime_error(std::string("Error writing efficiency file: ") + outFileName);
    outFile->Close();
    if (outFile->TestBit(TFile::kWriteError)) throw std::runtime_error(std::string("Error closing efficiency file: ") + outFileName);
    delete outFile;

    std::cout << "Efficiency curves created successfully!" << std::endl;

    // Write CSV output
    std::string csvFilename = std::string(outFileName);
    size_t      dotPos      = csvFilename.rfind('.');
    if (dotPos != std::string::npos)
        csvFilename = csvFilename.substr(0, dotPos) + ".csv";
    else
        csvFilename += ".csv";

    std::ofstream csvFile(csvFilename);
    if (!csvFile) throw std::runtime_error("Cannot create CSV: " + csvFilename);
    csvFile << "energy,detector,area_type,area,error,efficiency,efficiency_error\n";

    for (const auto& histPair : histogramData)
    {
        const std::string& histName = histPair.first;

        std::string detector, area_type;
        if (histName.size() > 8 && histName.substr(histName.size() - 8) == "_addback")
        {
            detector  = histName.substr(0, histName.size() - 8);
            area_type = "addback";
        } else if (histName.size() > 4 && histName.substr(histName.size() - 4) == "_sum")
        {
            detector  = histName.substr(0, histName.size() - 4);
            area_type = "sum";
        } else
        {
            detector  = histName;
            area_type = "single";
        }

        for (const auto& point : histPair.second) { csvFile << point.first << "," << detector << "," << area_type << "," << point.second << "," << CalculateSimulationCountError(point.second) << "," << point.second / kSimulatedGammas << "," << CalculateSimulationCountError(point.second) / kSimulatedGammas << "\n"; }
    }

    csvFile.close();
    if (!csvFile) throw std::runtime_error("Error writing CSV: " + csvFilename);
    std::cout << "CSV data written to: " << csvFilename << std::endl;
}