#include "efficiency_plots.hpp"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TStyle.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace EfficiencyPlots
{
////////////////////////////////////////////////////////////////////////////////
/// Strip leading and trailing whitespace.
std::string TrimWhitespace(const std::string& value)
{
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    return value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1);
}

////////////////////////////////////////////////////////////////////////////////
/// Split one CSV line, respecting quoted fields.
std::vector<std::string> ParseCSVRow(const std::string& line)
{
    std::vector<std::string> result;
    std::string              field;
    bool                     quoted = false;
    for (size_t i = 0; i < line.size(); ++i)
    {
        const char c = line[i];
        if (c == '"')
        {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"')
            {
                field += '"';
                ++i;
            } else
                quoted = !quoted;
        } else if (c == ',' && !quoted)
        {
            result.push_back(TrimWhitespace(field));
            field.clear();
        } else
            field += c;
    }
    if (quoted) throw std::runtime_error("Unterminated quoted CSV field");
    result.push_back(TrimWhitespace(field));
    return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Parse a finite numeric CSV value with a contextual error.
double ParseNumber(const std::string& text, const std::string& context)
{
    try
    {
        size_t       used  = 0;
        const double value = std::stod(text, &used);
        if (used != text.size() || !std::isfinite(value)) throw std::invalid_argument("number");
        return value;
    } catch (const std::exception&) { throw std::runtime_error("Invalid number '" + text + "' in " + context); }
}

////////////////////////////////////////////////////////////////////////////////
/// Index the CSV columns and reject duplicate names.
std::map<std::string, size_t> BuildColumnIndex(const std::vector<std::string>& fields)
{
    std::map<std::string, size_t> result;
    for (size_t i = 0; i < fields.size(); ++i)
        if (!result.emplace(fields[i], i).second) throw std::runtime_error("Duplicate CSV column: " + fields[i]);
    return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a required CSV column index.
size_t FindColumn(const std::map<std::string, size_t>& columns, const std::string& name)
{
    const auto it = columns.find(name);
    if (it == columns.end()) throw std::runtime_error("Missing CSV column: " + name);
    return it->second;
}

////////////////////////////////////////////////////////////////////////////////
/// Sort points by energy and reject duplicate energies.
void SortCurveByEnergy(Curve_t& curve)
{
    std::sort(curve.begin(), curve.end(), [](const Point_t& a, const Point_t& b) { return a.fEnergy < b.fEnergy; });
    for (size_t i = 1; i < curve.size(); ++i)
        if (curve[i - 1].fEnergy == curve[i].fEnergy) throw std::runtime_error("Duplicate energy in a detector/source curve: " + std::to_string(curve[i].fEnergy));
}

////////////////////////////////////////////////////////////////////////////////
/// Read detector counts and construct intensity-corrected source curves.
Experiment_t ReadExperimentalData(const std::string& path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open experimental CSV: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("Empty experimental CSV: " + path);
    const auto       fields       = ParseCSVRow(line);
    const auto       columns      = BuildColumnIndex(fields);
    const auto       energy       = FindColumn(columns, "Gamma Energy (keV)");
    const auto       energyUnc    = FindColumn(columns, "Energy Unc. (keV)");
    const auto       intensity    = FindColumn(columns, "Intensity (%)");
    const auto       intensityUnc = FindColumn(columns, "σ(I) (%)");
    const auto       source       = FindColumn(columns, "Source");
    const std::regex detectorPattern(R"(^([CB][0-9]+)(\+[CB][0-9]+)*$|^[A-Za-z]{1,2}$)");
    Experiment_t     result;
    for (size_t i = source + 1; i < fields.size(); ++i)
    {
        const auto& name = fields[i];
        if (std::regex_match(name, detectorPattern))
        {
            FindColumn(columns, name + " Unc.");
            result.fDetectors.push_back(name);
        }
    }
    if (result.fDetectors.empty()) throw std::runtime_error("No detector columns in " + path);

    size_t row = 1;
    while (std::getline(input, line))
    {
        ++row;
        if (TrimWhitespace(line).empty()) continue;
        const auto values  = ParseCSVRow(line);
        const auto context = path + ":" + std::to_string(row);
        if (values.size() != fields.size()) throw std::runtime_error("Wrong number of fields in " + context);
        const double e           = ParseNumber(values[energy], context);
        const double de          = ParseNumber(values[energyUnc], context);
        const double probability = ParseNumber(values[intensity], context) / 100.0;
        const double dp          = ParseNumber(values[intensityUnc], context) / 100.0;
        if (e <= 0 || de < 0 || probability <= 0 || probability > 1 || dp < 0 || values[source].empty()) throw std::runtime_error("Invalid energy, intensity, uncertainty, or source in " + context);
        for (const auto& detector : result.fDetectors)
        {
            const auto& count = values[FindColumn(columns, detector)];
            const auto& error = values[FindColumn(columns, detector + " Unc.")];
            if (count.empty() && error.empty()) continue;
            if (count.empty() || error.empty()) throw std::runtime_error("Incomplete detector measurement in " + context);
            const double n  = ParseNumber(count, context);
            const double dn = ParseNumber(error, context);
            if (n <= 0 || dn <= 0) throw std::runtime_error("Counts and uncertainties must be positive in " + context);
            result.fRaw[detector][values[source]].push_back({e, de, n, dn});
            const double y  = n / probability;
            const double dy = y * std::hypot(dn / n, dp / probability);
            result.fCorrected[detector][values[source]].push_back({e, de, y, dy});
        }
    }
    for (auto& detector : result.fRaw)
        for (auto& sourceCurve : detector.second) SortCurveByEnergy(sourceCurve.second);
    for (auto& detector : result.fCorrected)
        for (auto& sourceCurve : detector.second) SortCurveByEnergy(sourceCurve.second);
    return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Read simulated clover and CeBr efficiency curves.
Simulations_t ReadSimulationData(const std::string& path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open simulation CSV: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("Empty simulation CSV: " + path);
    const auto       fields      = ParseCSVRow(line);
    const auto       columns     = BuildColumnIndex(fields);
    const auto       energy      = FindColumn(columns, "energy");
    const auto       detector    = FindColumn(columns, "detector");
    const auto       type        = FindColumn(columns, "area_type");
    const auto       efficiency  = FindColumn(columns, "efficiency");
    const auto       uncertainty = FindColumn(columns, "efficiency_error");
    const std::regex clover(R"(^[CB][0-9]+$)");
    const std::regex cebr(R"(^[A-Za-z]{1,2}$)");
    Simulations_t    result;
    size_t           row = 1;
    while (std::getline(input, line))
    {
        ++row;
        if (TrimWhitespace(line).empty()) continue;
        const auto values  = ParseCSVRow(line);
        const auto context = path + ":" + std::to_string(row);
        if (values.size() != fields.size()) throw std::runtime_error("Wrong number of fields in " + context);
        const bool isClover = std::regex_match(values[detector], clover);
        const bool isCebr   = std::regex_match(values[detector], cebr);
        if ((!isClover || (values[type] != "sum" && values[type] != "addback")) && (!isCebr || values[type] != "single")) continue;
        const double e  = ParseNumber(values[energy], context);
        const double y  = ParseNumber(values[efficiency], context);
        const double dy = ParseNumber(values[uncertainty], context);
        if (e <= 0 || y < 0 || dy < 0) throw std::runtime_error("Invalid simulation value in " + context);
        result[values[detector]][values[type]].push_back({e, 0, y, dy});
    }
    if (result.empty()) throw std::runtime_error("No detector simulation curves in " + path);
    for (auto& detector : result)
        for (auto& type : detector.second) SortCurveByEnergy(type.second);
    return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Interpolate a curve and its uncertainty without extrapolating.
bool InterpolateCurve(const Curve_t& curve, double x, double& y, double& ey)
{
    if (curve.empty() || x < curve.front().fEnergy || x > curve.back().fEnergy) return false;
    const auto it = std::lower_bound(curve.begin(), curve.end(), x, [](const Point_t& p, double value) { return p.fEnergy < value; });
    if (it != curve.end() && it->fEnergy == x)
    {
        y  = it->fValue;
        ey = it->fValueError;
        return true;
    }
    if (it == curve.begin() || it == curve.end()) return false;
    const auto&  a = *(it - 1);
    const auto&  b = *it;
    const double t = (x - a.fEnergy) / (b.fEnergy - a.fEnergy);
    y              = (1 - t) * a.fValue + t * b.fValue;
    ey             = std::hypot((1 - t) * a.fValueError, t * b.fValueError);
    return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Split a combined detector name into its individual detectors.
std::vector<std::string> SplitDetectorComponents(const std::string& name)
{
    std::vector<std::string> result;
    std::istringstream       input(name);
    std::string              component;
    while (std::getline(input, component, '+')) result.push_back(component);
    return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Sum matching simulated detector curves and their independent errors.
std::map<std::string, Curve_t> CombineDetectorSimulations(const std::string& name, const Simulations_t& simulation)
{
    const auto                     parts = SplitDetectorComponents(name);
    std::map<std::string, Curve_t> result;
    const auto                     first = simulation.find(parts.front());
    if (first == simulation.end()) throw std::runtime_error("Missing simulation for " + parts.front());
    result = first->second;
    for (size_t i = 1; i < parts.size(); ++i)
    {
        const auto detector = simulation.find(parts[i]);
        if (detector == simulation.end()) throw std::runtime_error("Missing simulation for " + parts[i]);
        if (detector->second.size() != result.size()) throw std::runtime_error("Simulation types differ for " + name);
        for (auto& type : result)
        {
            const auto other = detector->second.find(type.first);
            if (other == detector->second.end() || other->second.size() != type.second.size()) throw std::runtime_error("Simulation curves do not match for " + name);
            for (size_t j = 0; j < type.second.size(); ++j)
            {
                auto&       p = type.second[j];
                const auto& q = other->second[j];
                if (p.fEnergy != q.fEnergy) throw std::runtime_error("Simulation energy grids differ for " + name);
                p.fValue += q.fValue;
                p.fValueError = std::hypot(p.fValueError, q.fValueError);
            }
        }
    }
    return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Align source curves using weighted interpolation over their common energies.
Alignment_t AlignSourceCurves(const Sources_t& sources)
{
    if (sources.empty()) throw std::runtime_error("No experimental sources to align");
    Alignment_t result;
    auto        reference = sources.begin();
    for (auto it = sources.begin(); it != sources.end(); ++it)
        if (it->second.size() > reference->second.size()) reference = it;
    result.fReference                 = reference->first;
    result.fFactors[reference->first] = 1;
    result.fCurves[reference->first]  = reference->second;
    for (const auto& source : sources)
    {
        if (source.first == reference->first) continue;
        double weightedLog = 0, weights = 0;
        size_t matches = 0;
        for (const auto& point : source.second)
        {
            double y = 0, ey = 0;
            if (!InterpolateCurve(reference->second, point.fEnergy, y, ey) || y <= 0) continue;
            const double variance = std::pow(point.fValueError / point.fValue, 2) + std::pow(ey / y, 2);
            if (variance <= 0 || !std::isfinite(variance)) continue;
            const double weight = 1 / variance;
            weightedLog += weight * std::log(y / point.fValue);
            weights += weight;
            ++matches;
        }
        if (matches == 0 || weights <= 0) throw std::runtime_error("No overlap for source " + source.first + " with reference " + reference->first);
        const double factor = std::exp(weightedLog / weights);
        if (!std::isfinite(factor) || factor <= 0) throw std::runtime_error("Invalid relative scale for " + source.first);
        result.fFactors[source.first] = factor;
        auto& curve                   = result.fCurves[source.first];
        for (const auto& p : source.second) curve.push_back({p.fEnergy, p.fEnergyError, p.fValue * factor, p.fValueError * factor});
    }
    return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Fit one simulation scale to the aligned experimental points.
double CalculateSimulationScale(const Sources_t& sources, const Curve_t& simulation)
{
    double numerator = 0, denominator = 0;
    for (const auto& source : sources)
        for (const auto& p : source.second)
        {
            double y = 0, ey = 0;
            if (!InterpolateCurve(simulation, p.fEnergy, y, ey) || p.fValueError <= 0) continue;
            const double weight = 1 / (p.fValueError * p.fValueError);
            numerator += weight * y * p.fValue;
            denominator += weight * y * y;
        }
    if (denominator <= 0 || numerator <= 0) throw std::runtime_error("Cannot scale simulation: no usable overlap");
    return numerator / denominator;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a styled ROOT graph with scaled values and uncertainties.
TGraphErrors* CreateGraph(const Curve_t& curve, const std::string& name, int color, int marker, double scale = 1)
{
    auto* result = new TGraphErrors(static_cast<int>(curve.size()));
    result->SetName(name.c_str());
    result->SetTitle(name.c_str());
    result->SetMarkerColor(color);
    result->SetLineColor(color);
    result->SetMarkerStyle(marker);
    result->SetMarkerSize(0.8);
    result->SetLineWidth(2);
    for (size_t i = 0; i < curve.size(); ++i)
    {
        const auto& p = curve[i];
        result->SetPoint(static_cast<int>(i), p.fEnergy, p.fValue * scale);
        result->SetPointError(static_cast<int>(i), p.fEnergyError, p.fValueError * scale);
    }
    return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Write source graphs and their canvas, optionally overlaying simulation curves.
void WriteSourcePlot(TDirectory* dir, const std::string& name, const std::string& title, const Sources_t& sources, const std::map<std::string, Curve_t>* simulation = nullptr, double scale = 1)
{
    dir->cd();
    TCanvas     canvas(name.c_str(), title.c_str(), 1000, 700);
    TMultiGraph multi;
    multi.SetTitle(title.c_str());
    TLegend   legend(0.7, 0.65, 0.9, 0.9);
    const int colors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2};
    int       index    = 0;
    for (const auto& source : sources)
    {
        auto* g = CreateGraph(source.second, name + "_" + source.first, colors[index % 6], 20 + index % 6);
        g->Write();
        multi.Add(g, "P");
        legend.AddEntry(g, source.first.c_str(), "pe");
        ++index;
    }
    if (simulation)
        for (const auto& type : *simulation)
        {
            const int color = type.first == "sum" ? kOrange + 7 : kBlack;
            auto*     g     = CreateGraph(type.second, name + "_simulation_" + type.first, color, 1, scale);
            g->SetLineStyle(type.first == "sum" ? 2 : 1);
            g->Write();
            multi.Add(g, "L");
            legend.AddEntry(g, ("simulation " + type.first).c_str(), "l");
        }
    multi.Draw("A");
    legend.Draw();
    canvas.Write();
}

////////////////////////////////////////////////////////////////////////////////
/// Write one detector's simulation, experimental graphs, and normalization.
void WriteDetectorPlots(TDirectory* parent, const std::string& name, const Simulations_t& simulation, const Experiment_t& experiment)
{
    auto* dir = parent->mkdir(name.c_str());
    if (!dir) throw std::runtime_error("Cannot create ROOT directory for " + name);
    const auto curves = CombineDetectorSimulations(name, simulation);
    if (curves.empty()) throw std::runtime_error("No simulation for " + name);
    dir->cd();
    TCanvas     canvas("simulation", (name + " simulation").c_str(), 1000, 700);
    TMultiGraph multi;
    multi.SetTitle((name + " simulated efficiency;Energy (keV);Efficiency (full-energy / simulated)").c_str());
    TLegend legend(0.7, 0.72, 0.9, 0.9);
    for (const auto& type : curves)
    {
        const int color = type.first == "sum" ? kOrange + 7 : kBlack;
        auto*     g     = CreateGraph(type.second, "simulation_" + type.first, color, 20);
        g->SetLineStyle(type.first == "sum" ? 2 : 1);
        g->Write();
        multi.Add(g, "LP");
        legend.AddEntry(g, type.first.c_str(), "lp");
    }
    multi.Draw("A");
    legend.Draw();
    canvas.Write();

    const auto raw = experiment.fRaw.find(name);
    if (raw == experiment.fRaw.end()) return;
    WriteSourcePlot(dir, "raw_counts", name + " raw counts;Energy (keV);Peak counts", raw->second);
    const auto corrected = experiment.fCorrected.find(name);
    if (corrected == experiment.fCorrected.end()) throw std::runtime_error("Missing corrected data for " + name);
    const auto aligned = AlignSourceCurves(corrected->second);
    const auto addback = curves.find("addback");
    const auto single  = curves.find("single");
    const auto fit     = addback != curves.end() ? addback : single;
    if (fit == curves.end()) throw std::runtime_error("No addback/single simulation to scale for " + name);
    const double scale = CalculateSimulationScale(aligned.fCurves, fit->second);
    WriteSourcePlot(dir, "relative", name + " relative intensity-corrected counts;Energy (keV);Relative counts / gamma emission probability", aligned.fCurves, &curves, scale);
    std::ostringstream metadata;
    metadata.precision(17);
    metadata << "reference=" << aligned.fReference << "; simulation_scale=" << scale << "; simulation_scale_curve=" << fit->first;
    for (const auto& factor : aligned.fFactors) metadata << "; source_scale_" << factor.first << "=" << factor.second;
    dir->cd();
    TNamed("normalization", metadata.str().c_str()).Write();
}
} // namespace EfficiencyPlots

////////////////////////////////////////////////////////////////////////////////
/// Validate the input curves and write all detector directories to a ROOT file.
void EfficiencyPlots::WriteEfficiencyPlots(const std::string& rawFile, const std::string& simulationFile, const std::string& outputFile)
{
    using namespace EfficiencyPlots;
    const auto experiment = ReadExperimentalData(rawFile);
    const auto simulation = ReadSimulationData(simulationFile);
    // Validate inputs and alignments before opening the output.
    for (const auto& detector : experiment.fDetectors)
    {
        const auto curves = CombineDetectorSimulations(detector, simulation);
        const auto source = experiment.fCorrected.find(detector);
        if (source != experiment.fCorrected.end())
        {
            const auto aligned = AlignSourceCurves(source->second);
            const auto addback = curves.find("addback");
            const auto single  = curves.find("single");
            const auto fit     = addback != curves.end() ? addback : single;
            if (fit == curves.end()) throw std::runtime_error("No addback/single simulation for " + detector);
            CalculateSimulationScale(aligned.fCurves, fit->second);
        }
    }
    std::unique_ptr<TFile> output(TFile::Open(outputFile.c_str(), "RECREATE"));
    if (!output || output->IsZombie()) throw std::runtime_error("Cannot create output: " + outputFile);
    const std::regex crossName(R"(^C[0-9]+$)");
    const std::regex backName(R"(^B[0-9]+$)");
    for (const auto& detector : simulation) { WriteDetectorPlots(output.get(), detector.first, simulation, experiment); }
    for (const auto& detector : experiment.fDetectors)
    {
        if (detector.find('+') == std::string::npos) continue;
        const auto parts   = SplitDetectorComponents(detector);
        const bool isCross = std::regex_match(parts.front(), crossName);
        for (const auto& part : parts)
            if (!std::regex_match(part, isCross ? crossName : backName)) throw std::runtime_error("Mixed or invalid detector group: " + detector);
        WriteDetectorPlots(output.get(), detector, simulation, experiment);
    }
    output->Write();
    if (output->TestBit(TFile::kWriteError)) throw std::runtime_error("Error writing output ROOT file");
    output->Close();
    if (output->TestBit(TFile::kWriteError)) throw std::runtime_error("Error closing output ROOT file");
}
