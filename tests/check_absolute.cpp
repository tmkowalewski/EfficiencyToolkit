#include <TFile.h>
#include <TGraphErrors.h>
#include <TNamed.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

/// Check a value against the independently calculated expectation.
void CheckNear(double actual, double expected, const std::string& label)
{
    if (!std::isfinite(actual) || std::abs(actual - expected) > 1e-10 * std::max(std::abs(expected), 1e-15))
        throw std::runtime_error(label + " differs from the expected value");
}

/// Retrieve a required graph from a ROOT file.
TGraphErrors& RequireGraph(TFile& file, const char* path)
{
    auto* graph = file.Get<TGraphErrors>(path);
    if (!graph) throw std::runtime_error(std::string("Missing graph: ") + path);
    return *graph;
}

/// Compare existing relative data and verify the absolute calculation.
int main(int argc, char** argv)
{
    try
    {
        if (argc != 4) throw std::runtime_error("Usage: check_absolute relative.root absolute.root normalized.root");
        TFile relative(argv[1], "READ"), absolute(argv[2], "READ"), normalized(argv[3], "READ");
        if (relative.IsZombie() || absolute.IsZombie() || normalized.IsZombie()) throw std::runtime_error("Cannot open test ROOT files");
        auto& relativeGraph = RequireGraph(relative, "C1+C5/relative_152Eu");
        auto& unchangedGraph = RequireGraph(absolute, "C1+C5/relative_152Eu");
        if (relativeGraph.GetN() != unchangedGraph.GetN()) throw std::runtime_error("Relative graph count changed");
        for (int i = 0; i < relativeGraph.GetN(); ++i)
        {
            CheckNear(unchangedGraph.GetPointX(i), relativeGraph.GetPointX(i), "Relative energy");
            CheckNear(unchangedGraph.GetPointY(i), relativeGraph.GetPointY(i), "Relative value");
            CheckNear(unchangedGraph.GetErrorY(i), relativeGraph.GetErrorY(i), "Relative uncertainty");
        }
        if (relative.Get("C1+C5/absolute")) throw std::runtime_error("Absolute output unexpectedly enabled");
        auto& eu = RequireGraph(absolute, "C1+C5/absolute_152Eu");
        const double count = 1340356, countError = 2582.6, intensity = 0.0755, intensityError = 0.00004;
        const double expected = count / (intensity * 1000000 * 40000 * 0.8);
        const double expectedError = expected * std::sqrt(std::pow(countError / count, 2) +
            std::pow(intensityError / intensity, 2) + std::pow(50000.0 / 1000000, 2) + std::pow(0.04 / 0.8, 2));
        CheckNear(eu.GetPointX(0), 244.6974, "Eu energy");
        CheckNear(eu.GetErrorX(0), 0.001, "Eu energy uncertainty");
        CheckNear(eu.GetPointY(0), expected, "Eu efficiency");
        CheckNear(eu.GetErrorY(0), expectedError, "Eu uncertainty");
        auto& co = RequireGraph(absolute, "C1+C5/absolute_56Co");
        CheckNear(co.GetPointY(0), 18659.0 / (0.00191 * 2000000 * 30000 * 0.9), "Co efficiency");
        RequireGraph(absolute, "C1+C5/absolute_60Co");
        auto& sim = RequireGraph(absolute, "C1+C5/simulation_addback");
        auto& overlay = RequireGraph(absolute, "C1+C5/absolute_simulation_addback");
        if (sim.GetN() != overlay.GetN()) throw std::runtime_error("Simulation point count differs");
        for (int i = 0; i < sim.GetN(); ++i)
        {
            CheckNear(overlay.GetPointY(i), sim.GetPointY(i), "Unscaled simulation");
            CheckNear(overlay.GetErrorY(i), sim.GetErrorY(i), "Simulation uncertainty");
        }
        RequireGraph(absolute, "C1+C5/absolute_simulation_sum");
        if (absolute.Get("C1/absolute")) throw std::runtime_error("Unexpected absolute plot without data");
        auto* metadata = absolute.Get<TNamed>("C1+C5/absolute_exposure");
        if (!metadata || std::string(metadata->GetTitle()).find("status=measured") == std::string::npos)
            throw std::runtime_error("Missing exposure provenance");
        auto& normalizedCo = RequireGraph(normalized, "C1+C5/absolute_60Co");
        const double counts[] = {12443610, 11791180};
        const double countErrors[] = {7921.1, 7717.9};
        const double probabilities[] = {0.9985, 0.99983};
        const double probabilityErrors[] = {0.0003, 0.000006};
        double weightedLog = 0, weights = 0;
        int referencesUsed = 0;
        for (auto* reference : {&eu, &co})
        {
            bool used = false;
            for (int i = 0; i < normalizedCo.GetN(); ++i)
            {
                const double x = normalizedCo.GetPointX(i);
                for (int j = 0; j + 1 < reference->GetN(); ++j)
                {
                    if (x < reference->GetPointX(j) || x > reference->GetPointX(j + 1)) continue;
                    const double t = (x - reference->GetPointX(j)) / (reference->GetPointX(j + 1) - reference->GetPointX(j));
                    const double y = (1 - t) * reference->GetPointY(j) + t * reference->GetPointY(j + 1);
                    const double error = std::hypot((1 - t) * reference->GetErrorY(j), t * reference->GetErrorY(j + 1));
                    const double value = counts[i] / probabilities[i];
                    const double relError = std::hypot(countErrors[i] / counts[i], probabilityErrors[i] / probabilities[i]);
                    const double weight = 1 / (relError * relError + std::pow(error / y, 2));
                    weightedLog += weight * std::log(y / value);
                    weights += weight;
                    used = true;
                    break;
                }
            }
            if (used) ++referencesUsed;
        }
        if (normalizedCo.GetN() != 2 || referencesUsed != 2) throw std::runtime_error("Expected two points and references");
        const double factor = std::exp(weightedLog / weights);
        for (int i = 0; i < normalizedCo.GetN(); ++i)
        {
            const double value = counts[i] / probabilities[i] * factor;
            CheckNear(normalizedCo.GetPointY(i), value, "Normalized efficiency");
            CheckNear(normalizedCo.GetErrorY(i), value * std::hypot(countErrors[i] / counts[i],
                      probabilityErrors[i] / probabilities[i]), "Normalized uncertainty");
        }
        auto& unchangedEu = RequireGraph(normalized, "C1+C5/absolute_152Eu");
        CheckNear(unchangedEu.GetPointY(0), eu.GetPointY(0), "Measured efficiency after normalization");
        auto* normalizedMetadata = normalized.Get<TNamed>("C1+C5/absolute_exposure");
        if (!normalizedMetadata || std::string(normalizedMetadata->GetTitle()).find(
                "60Co: status=normalize, normalization_factor=") == std::string::npos ||
            std::string(normalizedMetadata->GetTitle()).find("references=152Eu+56Co") == std::string::npos)
            throw std::runtime_error("Missing normalization provenance");
        std::cout << "Absolute-efficiency integration checks passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
