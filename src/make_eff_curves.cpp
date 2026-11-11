#include "efficiency_plots.hpp"

#include <TROOT.h>
#include <algorithm>
#include <boost/program_options.hpp>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

////////////////////////////////////////S////////////////////////////////////////
/// Parse options and safely write the experimental efficiency plots.
int main(int argc, char** argv)
{
    namespace po = boost::program_options;
    namespace fs = std::filesystem;
    try
    {
        std::string             raw, simulation, outputName, sourceInfo;
        EfficiencyPlots::RadwareFitOptions_t fitOptions;
        po::options_description options("Options");
        options.add_options()("help,h", "Show this help")("raw,r", po::value<std::string>(&raw)->default_value("data/156Eu-56Co-60Co_RawCounts.csv"), "Experimental raw-count CSV")("simulation,s", po::value<std::string>(&simulation)->default_value("data/70Ge_cloverarray.eff.csv"), "Simulated efficiency CSV")("source-info", po::value<std::string>(&sourceInfo), "Source exposure CSV; enables absolute efficiencies")("debug", "Save individual component graphs as well as combined canvases")("radware-fit", "Fit source data using the RadWare efficiency function")("radware-b", po::value<double>(&fitOptions.fLowSlope), "Fix RadWare low-energy slope B (default: fit B)")("radware-g", po::value<double>(&fitOptions.fTurnover), "Fix RadWare turnover parameter G (default 20)")("radware-scale", po::value<double>(&fitOptions.fScale), "Fixed positive efficiency multiplier (default automatic per fit)")("output,o", po::value<std::string>(&outputName)->default_value("efficiency_plots.root"), "Output ROOT file; a filename alone is placed in efficiency_curves/");
        po::variables_map arguments;
        po::store(po::parse_command_line(argc, argv, options), arguments);
        po::notify(arguments);
        if (arguments.count("help"))
        {
            std::cout << "Usage: " << argv[0] << " [options]\n" << options << std::endl;
            return 0;
        }
        if (raw.empty() || simulation.empty() || outputName.empty()) throw std::runtime_error("Input and output paths must not be empty");
        if (arguments.count("source-info") && sourceInfo.empty()) throw std::runtime_error("Source exposure path must not be empty");
        if (!arguments.count("radware-fit") && (arguments.count("radware-b") || arguments.count("radware-g") ||
                                                arguments.count("radware-scale")))
            throw std::runtime_error("RadWare parameters require --radware-fit");
        fitOptions.fFixLowSlope = arguments.count("radware-b") != 0;
        if (!std::isfinite(fitOptions.fLowSlope) || !std::isfinite(fitOptions.fTurnover) || fitOptions.fTurnover <= 0 ||
            !std::isfinite(fitOptions.fScale) || (arguments.count("radware-scale") && fitOptions.fScale <= 0))
            throw std::runtime_error("RadWare B must be finite and G and scale must be finite and positive");
        if (!sourceInfo.empty() && !fs::exists(sourceInfo)) throw std::runtime_error("Cannot open source exposure CSV: " + sourceInfo);

        fs::path output(outputName);
        if (!output.has_parent_path()) output = fs::path("efficiency_curves") / output;
        if (output.filename().empty()) throw std::runtime_error("Output path must name a file");
        const auto status = fs::symlink_status(output);
        const bool exists = fs::exists(status);
        if (exists && fs::is_directory(output)) throw std::runtime_error("Output path is a directory: " + output.string());
        if (exists && (fs::equivalent(output, raw) || fs::equivalent(output, simulation) ||
                       (!sourceInfo.empty() && fs::equivalent(output, sourceInfo)))) throw std::runtime_error("Output must not overwrite an input file");
        if (exists)
        {
            std::cout << "Output " << output.string() << " already exists. Overwrite? [y/N] " << std::flush;
            std::string answer;
            if (!std::getline(std::cin, answer))
            {
                std::cerr << "No confirmation received; output was not changed." << std::endl;
                return 1;
            }
            answer = EfficiencyPlots::TrimWhitespace(answer);
            std::transform(answer.begin(), answer.end(), answer.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (answer != "y" && answer != "yes")
            {
                std::cout << "Output was not changed." << std::endl;
                return 1;
            }
        }
        fs::create_directories(output.parent_path());
        std::string       tempTemplate = (output.parent_path() / ("." + output.filename().string() + ".XXXXXX")).string();
        std::vector<char> tempName(tempTemplate.begin(), tempTemplate.end());
        tempName.push_back('\0');
        const int fd = mkstemp(tempName.data());
        if (fd == -1) throw std::runtime_error("Cannot create temporary output in " + output.parent_path().string());
        const fs::path temporary(tempName.data());
        if (close(fd) != 0)
        {
            fs::remove(temporary);
            throw std::runtime_error("Cannot close temporary output: " + temporary.string());
        }
        try
        {
            gROOT->SetBatch(kTRUE);
            EfficiencyPlots::WriteEfficiencyPlots(raw, simulation, temporary.string(), sourceInfo,
                                                   arguments.count("radware-fit") != 0, fitOptions, arguments.count("debug") != 0);
            if (exists)
                fs::rename(temporary, output);
            else
            {
                fs::create_hard_link(temporary, output);
                fs::remove(temporary);
            }
        } catch (...)
        {
            std::error_code cleanupError;
            fs::remove(temporary, cleanupError);
            if (cleanupError) std::cerr << "Cannot remove temporary output " << temporary << ": " << cleanupError.message() << std::endl;
            throw;
        }
        std::cout << "Wrote " << output.string() << std::endl;
        return 0;
    } catch (const std::exception& error)
    {
        std::cerr << "Efficiency plots: " << error.what() << std::endl;
        return 1;
    }
}