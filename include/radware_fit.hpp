#ifndef EFFICIENCY_TOOLKIT_RADWARE_FIT_HPP
#define EFFICIENCY_TOOLKIT_RADWARE_FIT_HPP

#include "efficiency_plots.hpp"

#include <TF1.h>
#include <TFitResultPtr.h>
#include <memory>

namespace EfficiencyPlots {

/// A RadWare curve and its ROOT minimization result.
struct RadwareFit_t {
   std::unique_ptr<TF1> fFunction;
   TFitResultPtr fResult;
};

/// Fit combined source points with RadWare's efficiency parameterization.
RadwareFit_t FitRadwareEfficiency(const Sources_t& sources, const std::string& name,
                                 const RadwareFitOptions_t& options = {});

} // namespace EfficiencyPlots
#endif
