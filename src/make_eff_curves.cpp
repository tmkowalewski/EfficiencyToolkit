#include "efficiency_plots.hpp"

#include <TROOT.h>
#include <boost/program_options.hpp>
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>

////////////////////////////////////////////////////////////////////////////////
/// Parse options and safely write the experimental efficiency plots.
int main(int argc, char **argv)
{
   namespace po = boost::program_options;
   namespace fs = std::filesystem;
   try {
      std::string raw, simulation, outputName;
      po::options_description options("Options");
      options.add_options()("help,h", "Show this help")(
         "raw,r", po::value<std::string>(&raw)->default_value("data/156Eu-56Co-60Co_RawCounts.csv"),
         "Experimental raw-count CSV")(
         "simulation,s", po::value<std::string>(&simulation)->default_value("data/70Ge_cloverarray.eff.csv"),
         "Simulated efficiency CSV")("output,o",
                                     po::value<std::string>(&outputName)->default_value("efficiency_plots.root"),
                                     "Output ROOT file; a filename alone is placed in efficiency_curves/");
      po::variables_map arguments;
      po::store(po::parse_command_line(argc, argv, options), arguments);
      po::notify(arguments);
      if (arguments.count("help")) {
         std::cout << "Usage: " << argv[0] << " [options]\n" << options << std::endl;
         return 0;
      }
      if (raw.empty() || simulation.empty() || outputName.empty())
         throw std::runtime_error("Input and output paths must not be empty");

      fs::path output(outputName);
      if (!output.has_parent_path())
         output = fs::path("efficiency_curves") / output;
      if (output.filename().empty())
         throw std::runtime_error("Output path must name a file");
      const auto status = fs::symlink_status(output);
      const bool exists = fs::exists(status);
      if (exists && fs::is_directory(output))
         throw std::runtime_error("Output path is a directory: " + output.string());
      if (exists && (fs::equivalent(output, raw) || fs::equivalent(output, simulation)))
         throw std::runtime_error("Output must not overwrite an input file");
      if (exists) {
         std::cout << "Output " << output.string() << " already exists. Overwrite? [y/N] " << std::flush;
         std::string answer;
         if (!std::getline(std::cin, answer)) {
            std::cerr << "No confirmation received; output was not changed." << std::endl;
            return 1;
         }
         answer = EfficiencyPlots::TrimWhitespace(answer);
         std::transform(answer.begin(), answer.end(), answer.begin(),
                        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
         if (answer != "y" && answer != "yes") {
            std::cout << "Output was not changed." << std::endl;
            return 1;
         }
      }
      fs::create_directories(output.parent_path());
      std::string tempTemplate = (output.parent_path() / ("." + output.filename().string() + ".XXXXXX")).string();
      std::vector<char> tempName(tempTemplate.begin(), tempTemplate.end());
      tempName.push_back('\0');
      const int fd = mkstemp(tempName.data());
      if (fd == -1)
         throw std::runtime_error("Cannot create temporary output in " + output.parent_path().string());
      const fs::path temporary(tempName.data());
      if (close(fd) != 0) {
         fs::remove(temporary);
         throw std::runtime_error("Cannot close temporary output: " + temporary.string());
      }
      try {
         gROOT->SetBatch(kTRUE);
         EfficiencyPlots::WriteEfficiencyPlots(raw, simulation, temporary.string());
         if (exists)
            fs::rename(temporary, output);
         else {
            fs::create_hard_link(temporary, output);
            fs::remove(temporary);
         }
      } catch (...) {
         std::error_code cleanupError;
         fs::remove(temporary, cleanupError);
         if (cleanupError)
            std::cerr << "Cannot remove temporary output " << temporary << ": " << cleanupError.message() << std::endl;
         throw;
      }
      std::cout << "Wrote " << output.string() << std::endl;
      return 0;
   } catch (const std::exception &error) {
      std::cerr << "Efficiency plots: " << error.what() << std::endl;
      return 1;
   }
}