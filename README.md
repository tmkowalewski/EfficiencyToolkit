# Efficiency Toolkit

Build with a C++ compiler, ROOT (including `root-config`), and the Boost
Program Options development library:

```sh
make -C EfficiencyToolkit
```

Run the standalone executables from the parent efficiency calibration
directory, where `data/`, `trees/`, `hists/`, and `efficiency_curves/`
are located:

```sh
./EfficiencyToolkit/make_eff_curves
```

By default this reads `data/156Eu-56Co-60Co_RawCounts.csv` and
`data/70Ge_cloverarray.eff.csv`, and creates
`efficiency_curves/efficiency_plots.root`. To choose different files:

```sh
./EfficiencyToolkit/make_eff_curves --raw raw.csv --simulation simulation.csv --output plots.root
```

A bare output filename is placed in `efficiency_curves/`; an output path with
a directory is used as specified. Paths, including the defaults, are relative
to the current working directory. Use `--help` for all options. If the output
already exists, the program asks for `y` or `yes` before replacing it; a
negative answer or end-of-input leaves the existing file unchanged. Output is
first written to a temporary file in the destination directory and moved into
place after successful generation.

The output contains one top-level directory per individual detector; clover
pairs appearing in the experimental CSV have additional top-level directories,
such as `C1+C5`. Previously generated files retain their original grouped
layout until regenerated. Clover names match `C` or `B` followed by digits;
single- or two-letter simulation
names are classified as CeBr. Longer names and individual clover crystals are
not included. Check the names if the detector configuration changes.

Each detector has a `simulation` canvas and separately saved
`simulation_addback` and `simulation_sum` graphs (or `simulation_single` for
CeBr). Pair curves are the pointwise sums of the individual efficiencies.
Their displayed uncertainties are added in quadrature, assuming independence;
correlated simulation uncertainties are not accounted for. The sum and
addback curves refer to each clover's crystal sum and addback, respectively,
not an event-level addback between separate clovers.

Where experimental data exist, `raw_counts` shows the original peak counts and
their supplied energy and count uncertainties, grouped by source. `relative`
shows the aligned data with both simulation curves; its component graphs are
also saved individually. Individual detectors without experimental columns
have simulation plots only. Blank measurements are omitted.

## Relative normalization

These are **relative**, not absolute, efficiencies:

1. Divide each peak count by its gamma intensity as a fraction per source
   decay. The intensity column is treated as the total photon emission
   probability per decay; the separate `BR` column is not multiplied again.
   Count and intensity uncertainties are propagated independently.
2. Select the source with the most measured points for that detector as the
   reference. Linearly interpolate its corrected counts at the other sources'
   energies, without extrapolation. For each other source, calculate a single
   multiplicative factor from the inverse-variance weighted mean of the log
   ratios over overlapping energies. Each source must overlap the reference.
3. Fit one scale for the simulated addback curve to the aligned experimental
   points, using experimental point uncertainties as weights and linear
   interpolation within the simulated energy range. For CeBr, use the single
   curve. Apply the same scale to both simulated clover curves.

The reference source and factors are saved in each detector's `normalization`
object. Relative point errors include the input count and intensity errors,
but not the uncertainty or correlations in the fitted source factors,
interpolation, or simulation scale. Simulation error bars are saved in the
graphs but are not used in the overlay scale fit.

The source alignment empirically absorbs differences in run time, source
activity, and other source-dependent normalization factors. It is not a
measured live-time correction and does not establish absolute efficiency.
No live times, activities, attenuation, coincidence-summing corrections, or
decay corrections are assumed. The alignment depends on the interpolated
reference data, so sparse overlap or source-specific systematic effects can
bias it. The raw experimental clover peaks are assumed to be addback counts;
the crystal-sum curve is shown for comparison.

## Simulation histograms

The separate `get_sim_efficiencies` executable requires explicit actions:

```sh
# Sort trees, then generate the ROOT graph file and CSV.
./EfficiencyToolkit/get_sim_efficiencies --sort --curves --threads 4

# Or perform just one step.
./EfficiencyToolkit/get_sim_efficiencies --sort
./EfficiencyToolkit/get_sim_efficiencies --curves
```

By default, sorting reads `trees/gamma.NNNNN.root` (tree `edep`) and writes
histograms to `hists/`. Curve generation reads all matching histogram files in
that directory and writes `data/70Ge_cloverarray.eff.root` and the adjacent
`data/70Ge_cloverarray.eff.csv`, the simulation input to `make_eff_curves`.
Use `--trees`, `--histograms`, `--tree-name`, `--output`, and `--threads` to
choose other paths or settings; see `--help`. Paths are relative to the
current working directory. When sorting only a subset, use a separate
histogram directory: curve generation includes *all* matching files present
there, including any from earlier runs.

Existing matching outputs require confirmation before either action starts.
Histogram and simulation output files are written in place, not transactionally;
if a run fails, some outputs may be incomplete. Use a fresh output directory
to preserve previous results.

For addback, each crystal contributes only if its energy is finite and
strictly greater than 150 keV. Earlier versions did not exclude NaNs or apply
the threshold correctly, so their addback histograms and efficiencies are
not equivalent to the corrected results.

## Layout regression check

When a grouped file from the earlier output layout is available, compare its
detectors, graph values, uncertainties, and normalization metadata with a
newly generated flat file:

```sh
make -C EfficiencyToolkit compare_layout
./EfficiencyToolkit/compare_layout grouped.root flat.root
```
