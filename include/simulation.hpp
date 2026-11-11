#ifndef EFFICIENCY_TOOLKIT_SIMULATION_HPP
#define EFFICIENCY_TOOLKIT_SIMULATION_HPP

#include <RtypesCore.h>
#include <atomic>
#include <mutex>
#include <string>

/// Thread-safe progress information for a batch of simulation trees.
struct SortProgress_t {
   std::atomic<Long64_t> fFileDone{0}, fFileTotal{0};
   std::atomic<Long64_t> fAllDone{0}, fAllTotal{0};
   std::atomic<long long> fFileStartNs{0}, fAllStartNs{0};
   std::atomic<bool> fStop{false};
   std::mutex fNameMutex;
   std::string fFileName;
   Int_t fFileIndex = 0, fFileCount = 0;
};

/// Input, output, and entry count for one simulation sort.
struct SortJob_t {
   std::string fInputPath, fOutputPath, fName;
   Long64_t fEntries;
};

/// Sort one simulation tree into detector histograms.
void SortSimulationTree(const char *fileName, const char *treeName, const char *outputFile,
                        SortProgress_t *progress = nullptr);
/// Sort every matching simulation tree in a directory.
void SortSimulationTrees(const char *directoryName, const char *outputDirectory, const char *treeName);
/// Extract simulated efficiency curves and a CSV from detector histograms.
void CreateSimulationEfficiencyCurves(const char *directoryName, const char *outputFile);

#endif
