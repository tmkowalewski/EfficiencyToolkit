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

/// Run-averaged activity and exposure for an experimental source.
struct SourceExposure_t {
   Double_t fActivityBq = 0, fActivityErrorBq = 0, fElapsedSeconds = 0, fLiveFraction = 0, fLiveFractionError = 0;
   std::string fStatus;
};

using SourceExposures_t = std::map<std::string, SourceExposure_t>;

/// RadWare shape settings and an optional efficiency scale.
struct RadwareFitOptions_t {
   Double_t fLowSlope = 1;
   Double_t fTurnover = 20;
   Double_t fScale = 0;
   bool fFixLowSlope = false;
};

/// Absolute curves and fitted scales for sources without measured exposures.
struct AbsoluteEfficiencies_t {
   Sources_t fCurves;
   std::map<std::string, Double_t> fNormalizationFactors;
   std::map<std::string, std::vector<std::string>> fNormalizationReferences;
};

/// Remove surrounding whitespace from a string.
std::string TrimWhitespace(const std::string &value);

/// Read the input data and write the detector plots to a ROOT file.
void WriteEfficiencyPlots(const std::string &rawFile, const std::string &simulationFile, const std::string &outputFile,
                          const std::string &sourceInfoFile = "", bool radwareFit = false,
                          const RadwareFitOptions_t &fitOptions = {}, bool debug = false);
} // namespace EfficiencyPlots

#endif
