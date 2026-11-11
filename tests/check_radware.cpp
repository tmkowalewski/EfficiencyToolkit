#include "radware_fit.hpp"

#include <TCanvas.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

/// Evaluate the documented RadWare expression independently of the ROOT formula.
double EvaluateRadware(double energy, const double* p)
{
    const double x = std::log(energy / 100), y = std::log(energy / 1000);
    const double a = p[0] + p[1] * x + p[2] * x * x;
    const double b = p[3] + p[4] * y + p[5] * y * y;
    return std::exp(std::pow(std::pow(a, -p[6]) + std::pow(b, -p[6]), -1 / p[6])) / p[7];
}

/// Verify that the function reproduces a synthetic, known RadWare curve.
void CheckSyntheticFit(double slope, double turnover)
{
    double p[] = {8.5, slope, 0, 9.5, -0.8, -0.04, turnover, 1};
    double lower = 8, upper = 12;
    for (int iteration = 0; iteration < 80; ++iteration)
    {
        p[3] = (lower + upper) / 2;
        double maximum = 0;
        for (int i = 0; i < 50; ++i)
            maximum = std::max(maximum, EvaluateRadware(60 * std::pow(50., i / 49.), p));
        if (maximum < 10000) lower = p[3]; else upper = p[3];
    }
    p[7] = 10000;
    EfficiencyPlots::Sources_t data;
    for (int i = 0; i < 50; ++i)
    {
        const double energy = 60 * std::pow(50., i / 49.);
        const double value = EvaluateRadware(energy, p);
        data[i % 2 ? "one" : "two"].push_back({energy, 0.1, value, value * 0.01});
    }
    const EfficiencyPlots::RadwareFitOptions_t options{slope, turnover, 10000, true};
    auto fit = EfficiencyPlots::FitRadwareEfficiency(data, "synthetic_radware_fit", options);
    if (fit.fResult->Chi2() > 1e-4) throw std::runtime_error("Synthetic RadWare fit did not recover the curve");
    for (int i : {0, 3, 4, 5})
        if (std::abs(fit.fFunction->GetParameter(i) - p[i]) > 0.01)
            throw std::runtime_error("Synthetic RadWare parameter mismatch");
    for (double energy : {75., 244., 1332., 2500.})
        if (std::abs(fit.fFunction->Eval(energy) / EvaluateRadware(energy, p) - 1) > 1e-4)
            throw std::runtime_error("Synthetic RadWare prediction mismatch");
    const EfficiencyPlots::RadwareFitOptions_t freeSlope{1., turnover, 10000, false};
    auto freeFit = EfficiencyPlots::FitRadwareEfficiency(data, "synthetic_free_slope", freeSlope);
    if (freeFit.fResult->Chi2() > 1e-4 || std::abs(freeFit.fFunction->GetParameter(1) - slope) > 0.01)
        throw std::runtime_error("Free-slope fit failed to recover the synthetic turnover");
}

/// Check that freeing B improves the weighted residual without changing input data.
void CompareFits(TFile& fitted, TFile& fixed)
{
    for (const char* detector : {"B2", "B3", "C1+C5", "C3+C7", "B1+B5"})
        for (const char* type : {"relative", "absolute"})
        {
            const std::string path = std::string(detector) + "/" + type + "_radware_fit_result";
            auto* freeResult = fitted.Get<TFitResult>(path.c_str());
            auto* fixedResult = fixed.Get<TFitResult>(path.c_str());
            if (!freeResult || !fixedResult || freeResult->Chi2() > fixedResult->Chi2() * 1.001)
                throw std::runtime_error("Free-slope fit worsened chi-square for " + path);
        }
    auto* freeResult = fitted.Get<TFitResult>("B2/relative_radware_fit_result");
    auto* fixedResult = fixed.Get<TFitResult>("B2/relative_radware_fit_result");
    if (freeResult->Chi2() > fixedResult->Chi2() * 0.5)
        throw std::runtime_error("Expected substantial B2 relative fit improvement");
}

/// Check saved functions, fit diagnostics, canvas overlay, and exact objective.
void CheckStoredFit(TFile& file, const char* type)
{
    const std::string dir = "C1+C5/", prefix = std::string(type);
    auto* fit = file.Get<TF1>((dir + prefix + "_radware_fit").c_str());
    auto* result = file.Get<TFitResult>((dir + prefix + "_radware_fit_result").c_str());
    auto* canvas = file.Get<TCanvas>((dir + prefix).c_str());
    if (!fit || !result || !canvas || !canvas->GetListOfPrimitives()->FindObject(fit->GetName()))
        throw std::runtime_error("Missing saved RadWare function, result or overlay");
    if (!result->IsValid() || result->CovMatrixStatus() < 2 || result->Ndf() <= 0 ||
        fit->GetParameter(2) != 0 || fit->GetParameter(6) != 20)
        throw std::runtime_error("Invalid stored fit result or fixed parameters");
    double chi2 = 0, minimum = 1e10, maximum = 0;
    unsigned int count = 0;
    for (const char* source : {"152Eu", "56Co", "60Co"})
    {
        auto* graph = file.Get<TGraphErrors>((dir + prefix + "_" + source).c_str());
        if (!graph) throw std::runtime_error("Missing source graph");
        for (int i = 0; i < graph->GetN(); ++i)
        {
            const double x = graph->GetPointX(i), y = graph->GetPointY(i);
            const double expected = EvaluateRadware(x, fit->GetParameters());
            if (std::abs(fit->Eval(x) / expected - 1) > 1e-12)
                throw std::runtime_error("Stored formula differs from RadWare function");
            chi2 += std::pow((y - expected) / graph->GetErrorY(i), 2);
            minimum = std::min(minimum, x);
            maximum = std::max(maximum, x);
            ++count;
        }
    }
    if (std::abs(chi2 - result->Chi2()) > 1e-6 * std::max(chi2, 1.) ||
        result->Ndf() != count - 5 || fit->GetXmin() != minimum || fit->GetXmax() != maximum)
        throw std::runtime_error("Stored fit objective, degrees of freedom or range differs");
}

/// Run synthetic and persisted fit checks.
int main(int argc, char** argv)
{
    try
    {
        if (argc != 3) throw std::runtime_error("Usage: check_radware output.root fixed.root");
        CheckSyntheticFit(1, 20);
        CheckSyntheticFit(0.8, 15);
        TFile file(argv[1], "READ");
        if (file.IsZombie()) throw std::runtime_error("Cannot open RadWare fit output");
        CheckStoredFit(file, "relative");
        CheckStoredFit(file, "absolute");
        TFile fixed(argv[2], "READ");
        if (fixed.IsZombie()) throw std::runtime_error("Cannot open fixed-slope output");
        CompareFits(file, fixed);
        std::cout << "RadWare fit checks passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
