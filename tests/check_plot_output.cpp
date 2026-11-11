#include <TCanvas.h>
#include <TColor.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TKey.h>
#include <TList.h>
#include <TMultiGraph.h>
#include <TROOT.h>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

/// Retrieve the multigraph embedded in a saved canvas.
TMultiGraph& GetMultiGraph(TDirectory& dir, const std::string& name)
{
    auto* canvas = dir.Get<TCanvas>(name.c_str());
    if (!canvas) throw std::runtime_error("Missing canvas " + name);
    for (auto* object : *canvas->GetListOfPrimitives())
        if (auto* multi = dynamic_cast<TMultiGraph*>(object)) return *multi;
    throw std::runtime_error("Missing multigraph in " + name);
}

/// Compare embedded data with separately saved debug graphs and their drawing style.
void CheckCanvas(TDirectory& normal, TDirectory& debug, const std::string& canvasName)
{
    auto& normalMulti = GetMultiGraph(normal, canvasName);
    auto& debugMulti = GetMultiGraph(debug, canvasName);
    if (normalMulti.GetListOfGraphs()->GetSize() != debugMulti.GetListOfGraphs()->GetSize())
        throw std::runtime_error("Graph count differs in " + canvasName);
    for (auto* link = normalMulti.GetListOfGraphs()->FirstLink(); link; link = link->Next())
    {
        auto* graph = dynamic_cast<TGraphErrors*>(link->GetObject());
        if (!graph) throw std::runtime_error("Unexpected graph type in " + canvasName);
        const std::string name = graph->GetName();
        if (normal.GetKey(name.c_str())) throw std::runtime_error("Unexpected normal-mode component " + name);
        auto* saved = debug.Get<TGraphErrors>(name.c_str());
        auto* embedded = dynamic_cast<TGraphErrors*>(debugMulti.GetListOfGraphs()->FindObject(name.c_str()));
        if (!saved || !embedded || graph->GetN() != saved->GetN() || graph->GetN() != embedded->GetN())
            throw std::runtime_error("Missing debug component " + name);
        for (int i = 0; i < graph->GetN(); ++i)
        {
            for (const auto* other : {saved, embedded})
                if (graph->GetPointX(i) != other->GetPointX(i) || graph->GetPointY(i) != other->GetPointY(i) ||
                    graph->GetErrorX(i) != other->GetErrorX(i) || graph->GetErrorY(i) != other->GetErrorY(i))
                    throw std::runtime_error("Graph data differ for " + name);
        }
        const bool simulation = canvasName == "simulation" || name.find(canvasName + "_simulation_") == 0;
        const std::string option = link->GetOption();
        if (simulation)
        {
            auto* color = gROOT->GetColor(graph->GetFillColor());
            if (option != "3L" || (graph->GetFillStyle() != 3004 && graph->GetFillStyle() != 3005) || !color)
                throw std::runtime_error("Missing line or shaded uncertainty band for " + name);
        }
        else if (option != "P") throw std::runtime_error("Experimental point style changed for " + name);
    }
}

/// Check both output modes across all detector directories.
int main(int argc, char** argv)
{
    try
    {
        if (argc != 3) throw std::runtime_error("Usage: check_plot_output normal.root debug.root");
        TFile normal(argv[1], "READ"), debug(argv[2], "READ");
        if (normal.IsZombie() || debug.IsZombie()) throw std::runtime_error("Cannot open plot outputs");
        int detectorCount = 0, experimentalCount = 0;
        for (auto* object : *normal.GetListOfKeys())
        {
            const std::string name = object->GetName();
            auto* normalDir = normal.GetDirectory(name.c_str());
            auto* debugDir = debug.GetDirectory(name.c_str());
            if (!normalDir || !debugDir) throw std::runtime_error("Missing detector " + name);
            ++detectorCount;
            CheckCanvas(*normalDir, *debugDir, "simulation");
            if (normalDir->GetKey("raw_counts"))
            {
                ++experimentalCount;
                for (const char* canvas : {"raw_counts", "relative", "absolute"})
                    CheckCanvas(*normalDir, *debugDir, canvas);
                for (const char* item : {"normalization", "absolute_exposure", "relative_radware_fit",
                                         "absolute_radware_fit", "relative_radware_fit_result", "absolute_radware_fit_result"})
                    if (!normalDir->GetKey(item) || !debugDir->GetKey(item))
                        throw std::runtime_error("Missing retained metadata or fit " + std::string(item));
            }
        }
        if (detectorCount != 23 || experimentalCount != 5)
            throw std::runtime_error("Unexpected detector counts");
        std::cout << "Verified embedded graphs and debug output in " << detectorCount << " detector directories\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
