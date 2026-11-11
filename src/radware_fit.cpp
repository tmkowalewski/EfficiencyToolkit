#include "radware_fit.hpp"

#include <TFitResult.h>
#include <TGraphErrors.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace EfficiencyPlots {

////////////////////////////////////////////////////////////////////////////////
/// Fit the RadWare function on a positive-log efficiency scale.
RadwareFit_t FitRadwareEfficiency(const Sources_t& sources, const std::string& name, const RadwareFitOptions_t& options)
{
    if (!std::isfinite(options.fLowSlope) || !std::isfinite(options.fTurnover) || options.fTurnover <= 0 ||
        !std::isfinite(options.fScale) || options.fScale < 0)
        throw std::runtime_error(name + ": invalid fixed RadWare parameters or scale");
    Curve_t points;
    for (const auto& source : sources)
        points.insert(points.end(), source.second.begin(), source.second.end());
    const size_t freeParameters = 4 + (options.fFixLowSlope ? 0 : 1);
    if (points.size() <= freeParameters)
        throw std::runtime_error(name + ": RadWare fit requires at least " + std::to_string(freeParameters + 1) + " points");
    std::sort(points.begin(), points.end(), [](const Point_t& a, const Point_t& b) { return a.fEnergy < b.fEnergy; });
    const double low = points.front().fEnergy, high = points.back().fEnergy;
    if (low <= 0 || high <= low) throw std::runtime_error(name + ": RadWare fit needs a positive energy range");
    double maxValue = 0;
    for (const auto& point : points)
    {
        if (!std::isfinite(point.fValue) || !std::isfinite(point.fValueError) || point.fValue <= 0 || point.fValueError <= 0)
            throw std::runtime_error(name + ": RadWare fit needs positive finite values and uncertainties");
        maxValue = std::max(maxValue, point.fValue);
    }
    const double scale = options.fScale > 0 ? options.fScale : 10000.0 / maxValue;
    if (!std::isfinite(scale) || scale <= 0) throw std::runtime_error(name + ": invalid RadWare efficiency scale");
    TGraphErrors graph(static_cast<int>(points.size()));
    for (size_t i = 0; i < points.size(); ++i)
    {
        graph.SetPoint(static_cast<int>(i), points[i].fEnergy, points[i].fValue);
        graph.SetPointError(static_cast<int>(i), 0, points[i].fValueError);
    }
    const std::string x = "log(x/100.)", y = "log(x/1000.)";
    const std::string a = "([0]+[1]*" + x + "+[2]*" + x + "*" + x + ")";
    const std::string b = "([3]+[4]*" + y + "+[5]*" + y + "*" + y + ")";
    const std::string formula = "exp(pow(pow(" + a + ",-[6])+pow(" + b + ",-[6]),-1./[6]))/[7]";
    const double xl = std::log(low / 100.), xh = std::log(high / 100.);
    const double yl = std::log(low / 1000.), yh = std::log(high / 1000.);
    RadwareFit_t best{nullptr, TFitResultPtr(-1)};
    const std::vector<double> slopes = options.fFixLowSlope ? std::vector<double>{options.fLowSlope}
                                                           : std::vector<double>{0.5, 1., 2.5, 5.};
    for (double slope : slopes)
    for (bool staged : options.fFixLowSlope ? std::vector<bool>{true} : std::vector<bool>{false, true})
    {
        auto function = std::make_unique<TF1>(name.c_str(), formula.c_str(), low, high);
        function->SetParNames("A", "B", "C", "D", "E", "F", "G", "scale");
        function->SetParameters(std::log(points.front().fValue * scale) - slope * xl + 1,
                                slope, 0., std::log(maxValue * scale), -0.8, -0.05, options.fTurnover, scale);
        if (staged) function->FixParameter(1, slope);
        function->FixParameter(2, 0.);
        function->FixParameter(6, options.fTurnover);
        function->FixParameter(7, scale);
        function->SetLineColor(kMagenta + 2);
        function->SetLineWidth(2);
        function->SetNpx(500);
        TFitResultPtr result = graph.Fit(function.get(), "QSNEX0");
        if (staged && !options.fFixLowSlope)
        {
            function->ReleaseParameter(1);
            result = graph.Fit(function.get(), "QSNEX0");
        }
        if (int(result) != 0 || !result.Get() || !result->IsValid() || result->CovMatrixStatus() < 2 ||
            !std::isfinite(result->Chi2())) continue;
        const double* p = function->GetParameters();
        const auto highBranch = [p](double v) { return p[3] + p[4] * v + p[5] * v * v; };
        double highMinimum = std::min(highBranch(yl), highBranch(yh));
        if (p[5] > 0)
        {
            const double vertex = -p[4] / (2 * p[5]);
            if (vertex > yl && vertex < yh) highMinimum = std::min(highMinimum, highBranch(vertex));
        }
        if (!std::isfinite(highMinimum) || highMinimum <= 0 ||
            std::min(p[0] + p[1] * xl, p[0] + p[1] * xh) <= 0) continue;
        if (!best.fFunction || result->Chi2() < best.fResult->Chi2())
            best = {std::move(function), result};
    }
    if (!best.fFunction)
        throw std::runtime_error(name + ": no converged RadWare fit with reliable covariance and positive branches; "
                                 "try fixing B, changing G, or providing more data");
    return best;
}

} // namespace EfficiencyPlots
