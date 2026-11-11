#include "simulation.hpp"

#include <TROOT.h>
#include <algorithm>
#include <boost/program_options.hpp>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <string>

//////////////////////////////////////S//////////////////////////////////////////
/// Parse options and run the requested simulation sorting and curve actions.
int main(int argc, char** argv)
{
    namespace po = boost::program_options;
    namespace fs = std::filesystem;
    try
    {
        std::string             trees, histograms, treeName, output;
        unsigned int            threads = 0;
        po::options_description options("Options");
        options.add_options()("help,h", "Show this help")("sort", "Sort trees into histogram files")("curves", "Create simulation efficiency ROOT and CSV files from histograms")("trees", po::value<std::string>(&trees)->default_value("trees"), "Input tree directory")("histograms", po::value<std::string>(&histograms)->default_value("hists"), "Histogram output/input directory")("tree-name", po::value<std::string>(&treeName)->default_value("edep"), "Tree name")(
            "output,o", po::value<std::string>(&output)->default_value("data/70Ge_cloverarray.eff.root"), "Efficiency ROOT output; CSV is written alongside it")("threads,j", po::value<unsigned int>(&threads)->default_value(0), "ROOT worker threads (0: ROOT default)");
        po::variables_map args;
        po::store(po::parse_command_line(argc, argv, options), args);
        po::notify(args);
        if (args.count("help"))
        {
            std::cout << "Usage: " << argv[0] << " [--sort] [--curves] [options]\n" << options << std::endl;
            return 0;
        }
        if (!args.count("sort") && !args.count("curves")) throw std::runtime_error("Specify --sort, --curves, or both; use --help for options");
        if (trees.empty() || histograms.empty() || treeName.empty() || output.empty()) throw std::runtime_error("Paths and tree name must not be empty");
        const fs::path rootOutput(output);
        if (args.count("curves") && rootOutput.extension() != ".root") throw std::runtime_error("Efficiency output must have a .root extension");
        auto csvOutput = rootOutput;
        csvOutput.replace_extension(".csv");

        bool overwrite = false;
        if (args.count("sort"))
        {
            if (!fs::is_directory(trees)) throw std::runtime_error("Tree directory does not exist: " + trees);
            const std::regex pattern(R"(^gamma\.(\d{5})\.root$)");
            bool             found = false;
            for (const auto& entry : fs::directory_iterator(trees))
            {
                std::smatch match;
                const auto  name = entry.path().filename().string();
                if (!entry.is_regular_file() || !std::regex_match(name, match, pattern)) continue;
                found             = true;
                const auto target = fs::path(histograms) / ("gamma." + match[1].str() + ".hists.root");
                if (fs::exists(target))
                {
                    if (fs::is_directory(target) || fs::equivalent(target, entry.path())) throw std::runtime_error("Invalid histogram output: " + target.string());
                    overwrite = true;
                }
            }
            if (!found) throw std::runtime_error("No gamma.ENERGY.root input trees in " + trees);
        }
        if (args.count("curves"))
        {
            if (!args.count("sort") && !fs::is_directory(histograms)) throw std::runtime_error("Histogram directory does not exist: " + histograms);
            if (rootOutput.filename().empty()) throw std::runtime_error("Efficiency output must name a file");
            if (std::regex_match(rootOutput.filename().string(), std::regex(R"(^gamma\.\d{5}\.hists\.root$)")) && fs::weakly_canonical(rootOutput.parent_path().empty() ? fs::path(".") : rootOutput.parent_path()) == fs::weakly_canonical(histograms)) throw std::runtime_error("Efficiency ROOT output must not overwrite a histogram input");
            if (fs::is_directory(histograms))
            {
                for (const auto& entry : fs::directory_iterator(histograms))
                    if (entry.is_regular_file() && ((fs::exists(rootOutput) && fs::equivalent(rootOutput, entry.path())) || (fs::exists(csvOutput) && fs::equivalent(csvOutput, entry.path())))) throw std::runtime_error("Efficiency output must not overwrite a histogram input");
            }
            if (fs::is_directory(trees))
            {
                for (const auto& entry : fs::directory_iterator(trees))
                    if (entry.is_regular_file() && ((fs::exists(rootOutput) && fs::equivalent(rootOutput, entry.path())) || (fs::exists(csvOutput) && fs::equivalent(csvOutput, entry.path())))) throw std::runtime_error("Efficiency output must not overwrite a tree input");
            }
            if (fs::exists(rootOutput) || fs::exists(csvOutput)) overwrite = true;
            if ((fs::exists(rootOutput) && fs::is_directory(rootOutput)) || (fs::exists(csvOutput) && fs::is_directory(csvOutput))) throw std::runtime_error("Efficiency output path is a directory");
        }
        if (overwrite)
        {
            std::cout << "One or more requested output files already exist. Overwrite? [y/N] " << std::flush;
            std::string answer;
            if (!std::getline(std::cin, answer)) throw std::runtime_error("No confirmation received; existing files were not changed");
            answer.erase(0, answer.find_first_not_of(" \t\r\n"));
            const auto last = answer.find_last_not_of(" \t\r\n");
            if (last != std::string::npos) answer.erase(last + 1);
            std::transform(answer.begin(), answer.end(), answer.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (answer != "y" && answer != "yes")
            {
                std::cout << "Existing files were not changed." << std::endl;
                return 1;
            }
        }
        gROOT->SetBatch(kTRUE);
        if (threads) ROOT::EnableImplicitMT(threads);
        if (args.count("sort")) SortSimulationTrees(trees.c_str(), histograms.c_str(), treeName.c_str());
        if (args.count("curves"))
        {
            if (rootOutput.has_parent_path()) fs::create_directories(rootOutput.parent_path());
            CreateSimulationEfficiencyCurves(histograms.c_str(), output.c_str());
        }
        return 0;
    } catch (const std::exception& error)
    {
        std::cerr << "Simulation efficiencies: " << error.what() << std::endl;
        return 1;
    }
}
