#include "radware_fit.hpp"

#include <TFitResult.h>
#include <TGraphErrors.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace EfficiencyPlots {

namespace {
/// Calculate the one-sigma covariance band for the fitted RadWare parameters.
std::unique_ptr<TGraphErrors> MakeRadwareUncertainty(const TF1& function, const TFitResult& result)
{
    constexpr int n = 501;
    auto band = std::make_unique<TGraphErrors>(n);
    const std::string name = std::string(function.GetName()) + "_uncertainty";
    band->SetName(name.c_str());
    band->SetTitle(name.c_str());
    band->SetFillColor(kMagenta + 2);
    band->SetFillStyle(3001);
    const double* p = function.GetParameters();
    for (int i = 0; i < n; ++i)
    {
        const double energy = function.GetXmin() + (function.GetXmax() - function.GetXmin()) * i / (n - 1);
        const double x = std::log(energy / 100.), y = std::log(energy / 1000.);
        const double a = p[0] + p[1] * x + p[2] * x * x;
        const double b = p[3] + p[4] * y + p[5] * y * y;
        if (!(a > 0) || !(b > 0)) throw std::runtime_error(name + ": nonpositive RadWare branch at " + std::to_string(energy) + " keV");
        const double value = function.Eval(energy);
        const double h = std::log(value * p[7]);
        const double logRatio = -p[6] * (std::log(a) - std::log(b));
        const double weightA = logRatio >= 0 ? 1 / (1 + std::exp(-logRatio))
                                              : std::exp(logRatio) / (1 + std::exp(logRatio));
        const double da = value * h * weightA / a;
        const double db = value * h * (1 - weightA) / b;
        const double gradient[] = {da, da * x, da * x * x, db, db * y, db * y * y};
        double variance = 0, magnitude = 0;
        for (int j = 0; j < 6; ++j)
            for (int k = 0; k < 6; ++k)
            {
                const double term = gradient[j] * result.CovMatrix(j, k) * gradient[k];
                variance += term;
                magnitude += std::abs(term);
            }
        if (!std::isfinite(value) || !std::isfinite(variance) || !std::isfinite(magnitude) ||
            variance < -100 * std::numeric_limits<double>::epsilon() * magnitude)
            throw std::runtime_error(name + ": invalid propagated covariance at " + std::to_string(energy) + " keV");
        band->SetPoint(i, energy, value);
        band->SetPointError(i, 0, std::sqrt(std::max(0., variance)));
    }
    return band;
}
} // namespace

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
    RadwareFit_t best{nullptr, TFitResultPtr(-1), nullptr, low, high};
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
            best = {std::move(function), result, nullptr, low, high};
    }
    if (!best.fFunction)
        throw std::runtime_error(name + ": no converged RadWare fit with reliable covariance and positive branches; "
                                 "try fixing B, changing G, or providing more data");
    const double displayHigh = std::max(high, 8000.);
    const double yhDisplay = std::log(displayHigh / 1000.);
    const double xhDisplay = std::log(displayHigh / 100.);
    const double* p = best.fFunction->GetParameters();
    const auto highBranch = [p](double v) { return p[3] + p[4] * v + p[5] * v * v; };
    double minimum = std::min(highBranch(yl), highBranch(yhDisplay));
    if (p[5] > 0)
    {
        const double vertex = -p[4] / (2 * p[5]);
        if (vertex > yl && vertex < yhDisplay) minimum = std::min(minimum, highBranch(vertex));
    }
    if (!std::isfinite(minimum) || minimum <= 0 ||
        std::min(p[0] + p[1] * xl, p[0] + p[1] * xhDisplay) <= 0)
        throw std::runtime_error(name + ": fitted RadWare function is invalid between the data and 8000 keV");
    best.fFunction->SetRange(low, displayHigh);
    best.fUncertainty = MakeRadwareUncertainty(*best.fFunction, *best.fResult);
    return best;
}

} // namespace EfficiencyPlots
