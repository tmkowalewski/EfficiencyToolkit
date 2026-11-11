#ifndef EFFICIENCY_TOOLKIT_EFFICIENCY_PLOTS_HPP
#define EFFICIENCY_TOOLKIT_EFFICIENCY_PLOTS_HPP

#include <RtypesCore.h>
#include <map>
#include <string>
#include <vector>

namespace EfficiencyPlots {

/// One energy and value, with their respective uncertainties.
struct Point_t {
   Double_t fEnergy, fEnergyError, fValue, fValueError;
};

using Curve_t = std::vector<Point_t>;
using Sources_t = std::map<std::string, Curve_t>;
using Simulations_t = std::map<std::string, std::map<std::string, Curve_t>>;

/// Experimental detector curves before and after intensity correction.
struct Experiment_t {
   std::map<std::string, Sources_t> fRaw;
   std::map<std::string, Sources_t> fCorrected;
   std::vector<std::string> fDetectors;
};

/// Source curves aligned to a common reference and their scale factors.
struct Alignment_t {
   Sources_t fCurves;
   std::map<std::string, Double_t> fFactors;
   std::string fReference;
};

/// Remove surrounding whitespace from a string.
std::string TrimWhitespace(const std::string &value);

/// Read the input data and write the detector plots to a ROOT file.
void WriteEfficiencyPlots(const std::string &rawFile, const std::string &simulationFile, const std::string &outputFile);
} // namespace EfficiencyPlots

#endif
