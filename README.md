# Efficiency Toolkit

Build with a C++ compiler, ROOT (including `root-config`), and the Boost
Program Options development library:

```sh
make -C EfficiencyToolkit
```

The executables are built in `EfficiencyToolkit/bin/`. Run
`make -C EfficiencyToolkit clean` to remove them.

Run the standalone executables from the parent efficiency calibration
directory, where `data/`, `trees/`, `hists/`, and `efficiency_curves/`
are located:

```sh
./EfficiencyToolkit/bin/make_eff_curves
```

By default this reads `data/156Eu-56Co-60Co_RawCounts.csv` and
`data/70Ge_cloverarray.eff.csv`, and creates
`efficiency_curves/efficiency_plots.root`. To choose different files:

```sh
./EfficiencyToolkit/bin/make_eff_curves --raw raw.csv --simulation simulation.csv --output plots.root
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

Each detector has a `simulation` canvas. Simulations are drawn as lines with
shaded one-standard-deviation uncertainty bands; the same style is used for
simulation overlays on the relative and absolute plots. Pair curves are the
pointwise sums of the individual efficiencies.
Their displayed uncertainties are added in quadrature, assuming independence;
correlated simulation uncertainties are not accounted for. The sum and
addback curves refer to each clover's crystal sum and addback, respectively,
not an event-level addback between separate clovers.

Where experimental data exist, `raw_counts` shows the original peak counts and
their supplied energy and count uncertainties, grouped by source. `relative`
shows the aligned data with both simulation curves. Individual detectors without experimental columns
have simulation plots only. Blank measurements are omitted.

Use `--debug` to additionally save individual component graphs such as
`simulation_addback`, `simulation_sum`, `raw_counts_152Eu` and
`relative_simulation_addback` as ROOT objects in each detector directory.
Without `--debug`, these graphs remain embedded in their combined canvases
but are not separate directory entries. Fit functions, fit results, and
normalization metadata remain separate objects in either mode.

## Absolute efficiency

Absolute efficiency is optional and is enabled by providing a source-exposure
CSV. The example file may contain placeholder values; verify every entry:

```sh
./EfficiencyToolkit/bin/make_eff_curves --source-info data/source_exposure_example.csv --output absolute_example.root
```

Replace any example numbers before using results scientifically. Each source
name must exactly match the `Source` column in the raw-count CSV; one row per
source is required. The exposure columns are:

| Column | Meaning |
|---|---|
| `source` | Source name, e.g. `152Eu` |
| `activity_bq` | Average source activity during the elapsed run, in decays/s |
| `activity_uncertainty_bq` | One-standard-deviation uncertainty in that activity |
| `elapsed_time_s` | Total elapsed (wall-clock) run time, in seconds, not live time |
| `live_fraction` | Live time divided by elapsed time, between 0 and 1 |
| `live_fraction_uncertainty` | One-standard-deviation uncertainty in the fraction |
| `status` | `measured`, `example`, or `normalize` |

Change `status` to `measured` only after replacing the example values. The
zero uncertainties in an example are placeholders, not evidence that the
quantities are exact. No uncertainty is assigned to elapsed time. Values and
status are saved in each detector's `absolute_exposure` object.

Use `normalize` when a source has no reliable exposure, for example
`60Co,,,,,,normalize`. Its exposure fields are ignored, including any values
already present. For each detector, the program interpolates every overlapping
`measured` absolute curve at the source's energies (without extrapolation) and
fits one multiplicative scale to the intensity-corrected counts using an
inverse-variance weighted mean of the log ratios. At least one measured curve
must overlap; `example` and other `normalize` sources are not references.
The scale and contributing reference sources are saved in `absolute_exposure`.
These points are empirically normalized to an absolute scale, not independently
measured absolute efficiencies.

For each measured peak, efficiency is
`counts / (gamma_emission_probability * activity_bq * elapsed_time_s * live_fraction)`.
The intensity is the emission probability per total source decay; the separate
decay-mode branching ratio is not applied again. The activity must already be
the average during the run, including any decay correction. Multiplication by
the live fraction assumes the live-time fraction is representative throughout
the run; it is not valid for arbitrary time-dependent dead time correlated
with source decay. No decay integration or coincidence-summing correction is
applied. Source-capsule attenuation is only included when explicitly enabled
as described below. These assumptions and the correctness of inputs determine
whether the result is physically absolute.

Point errors propagate the independent supplied count, gamma-intensity,
activity, and live-fraction uncertainties in quadrature. Activity and live
fraction errors are shared by all points from the same source; ROOT graph
error bars show marginal uncertainties only and do not encode correlations.
For `normalize` sources, error bars include only the scaled count and intensity
errors: fitted-scale uncertainty, interpolation correlations, and correlations
with the reference data are not included. Common exposure uncertainties and
shared interpolation points also make the fit weights approximate.
The `absolute` canvas contains one graph per source and overlays simulation
curves fitted to the absolute experimental points. A single scale is fitted
to the addback curve (or the single curve for CeBr) by inverse-variance
weighted least squares, interpolating the simulation at experimental energies
within its range. The same scale is applied to both clover simulation curves
and their error bands. The factor and fitted curve type are saved in
`absolute_simulation_normalization`; the simulation-only canvas remains
unscaled. The fit ignores simulation uncertainties and correlations between
experimental points, and the band does not include fitted-scale uncertainty.
This empirical scale is not a measurement of absolute simulation efficiency.
With `--debug`, separately saved objects are named
`absolute_<source>` and `absolute_simulation_<type>`. Clover experimental
counts are assumed to be addback; the sum simulation curve is also shown for
comparison. Existing raw and relative plots remain unchanged.

Run the absolute-efficiency integration and input-validation checks from the
project directory with `make -C EfficiencyToolkit test`.

## Relative normalization

These are **relative**, not absolute, efficiencies:

1. Divide each peak count by its gamma intensity as a fraction per source
   decay. The intensity column is treated as the total photon emission
   probability per decay; the separate `BR` column is not multiplied again.
   Count and intensity uncertainties are propagated independently. If capsule
   attenuation is enabled, correct the intensity-adjusted counts before alignment.
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
No live times, activities, coincidence-summing corrections, or decay
corrections are assumed in the relative normalization. Capsule attenuation
is used only when explicitly enabled. The alignment depends on the interpolated
reference data, so sparse overlap or source-specific systematic effects can
bias it. The raw experimental clover peaks are assumed to be addback counts;
the crystal-sum curve is shown for comparison.

## Source-capsule attenuation

The three template files in `data/` are intentionally incomplete. Fill in
`source_capsules_template.csv`, `detector_positions_template.csv`, and
`attenuation_coefficients_template.csv` with actual measured geometry and
material data before enabling this correction. No material coefficients or
unknown detector angles are assumed. All three inputs are required together:

```sh
./EfficiencyToolkit/bin/make_eff_curves --source-info data/source_exposure_example.csv \
  --capsules data/source_capsules_template.csv \
  --detector-positions data/detector_positions_template.csv \
  --attenuation-coefficients data/attenuation_coefficients_template.csv \
  --output attenuation_corrected.root
```

Each capsule row specifies the source name, homogeneous `material`,
`density_g_cm3`, cylinder `radius_cm`, full axial `thickness_cm`, and
`axis_theta_deg` and `axis_phi_deg`. The emission point is at the center of
the solid cylinder. These dimensions describe material from the source point
to the *outer surface*, not just the wall of an empty shell. The model is not
suitable for multilayer capsules, off-center sources, or distributed activity.

Detector rows specify the direction from the origin to each individual
detector's center. Angles are in degrees: theta is measured from the +z beam
axis, and phi from +x toward +y in a right-handed coordinate system. C1 is
at (theta, phi) = (90, 0); the cross detectors have theta = 90 degrees.
Other template positions are unspecified. Fill all positions used by the
experimental columns, including individual components of a pair. Fully blank
position rows are ignored, but incomplete or missing required positions fail.
Finite detector solid angle and its variation in material path are not modeled.

The coefficient CSV contains `material,energy_keV,mass_attenuation_cm2_g`,
with positive mass attenuation coefficients for each material. Supply verified
total photon attenuation coefficients covering every source gamma energy;
the program interpolates log-log and refuses extrapolation. Include enough
tabulated points to resolve any absorption edges. Units are keV,
cm²/g, g/cm³, and cm; use total attenuation rather than energy-absorption
coefficients. The source-to-detector center ray intersects the first cylinder
side or endcap at distance `L`; its uncollided transmission is
`exp(-mass_attenuation * density * L)`. Intensity-corrected counts and their
statistical uncertainties are divided by the transmission before relative
normalization and absolute-efficiency calculation. Raw-count plots and the
simulation curves are unchanged.

For combined counts, an exact correction requires the unknown individual
contributions. Instead, the program estimates transmission as the mean of the
individual transmissions weighted by interpolated per-detector simulated
addback efficiencies at that energy. This assumes the simulation describes
the relative detector contributions and does not already contain the source
capsule; it is an approximation, not an event-level reconstruction. Source
attenuation, geometry, material, and simulation-weight uncertainties and their
correlations are not included in point errors. Review these assumptions before
scientific use.

The ROOT file stores copies of the three CSVs as top-level
`attenuation_capsules_csv`, `attenuation_positions_csv`, and
`attenuation_coefficients_csv` objects, and
per-source, per-component path lengths in each experimental detector's
`attenuation_geometry` object. With no attenuation options, the previous
calculation remains unchanged. Run the focused checks with
`make -C EfficiencyToolkit test-attenuation`.

## RadWare efficiency fit

Use `--radware-fit` to fit all experimental source points together for each
detector, separately for the relative plot and, when `--source-info` is given,
the absolute plot:

```sh
./EfficiencyToolkit/bin/make_eff_curves --source-info data/source_exposure_example.csv \
  --radware-fit --output radware_efficiencies.root
```

The overlays and saved `relative_radware_fit` and `absolute_radware_fit` TF1
objects use the RadWare EFFIT function:

```
efficiency = exp(((A + B*x + C*x*x)^(-G) + (D + E*y + F*y*y)^(-G))^(-1/G)) / scale
x = ln(energy_keV / 100)
y = ln(energy_keV / 1000)
```

The RadWare expression requires positive logarithmic branches. The fixed
`scale` multiplies the input efficiency before applying the function and is
saved as its eighth parameter; the output is divided by that same value.
By default it is chosen separately for each fit so the largest input value is
10,000 after scaling. This choice affects the parameterization and fitted
shape, not just the units. Use `--radware-scale N` to fix a different positive
scale for every fit; it must work for both relative and absolute values when
both are fitted. The stored parameters are not directly usable as a RadWare
`.eff` file without accounting for this external scale.

By default `C=0` and `G=20` are fixed, while `A`, `B`, `D`, `E`, and `F`
are fitted. The fit tries multiple starting values for B and retains the
converged, positive-branch solution with the lowest chi-square. Set a different
fixed G with `--radware-g`, or fix B with `--radware-b` (for example,
`--radware-b 1` reproduces the earlier fixed-slope approach). Fixing B can be
useful when there are too few low-energy data to constrain the turnover.
At least six points are required when B is free, or five when it is fixed;
the fit also requires a usable covariance matrix. The fit itself uses only the
observed energy range, but the plotted and saved function extends from the
lowest data energy to 8000 keV (or the highest data energy if greater). Values
above the data are extrapolations and are not validated by the measurements.
The original data range is saved in `*_radware_fit_range`.

The ROOT `*_radware_fit_result` objects contain fit status, chi-square, degrees
of freedom, parameter errors and covariance. Fits minimize residuals in
efficiency divided by the graph's vertical point uncertainty; energy errors
and correlations are ignored. In particular, shared exposure errors and
uncertainties of source-normalization factors are not represented in the fit
covariance. A warning is printed if chi-square per degree of freedom exceeds 5
or the covariance had to be forced positive definite. Convergence alone does
not imply an adequate calibration; inspect residuals, fit range, sensitivity
to fixed parameters and scale, and the resulting curve before scientific use.
The `*_radware_fit_uncertainty` graph and shaded band show the pointwise
one-standard-deviation parameter uncertainty, calculated by linear propagation
of the full fitted-parameter covariance. Fixed parameters contribute no
uncertainty. This is not a prediction interval and does not include model
uncertainty, omitted correlations, source-normalization uncertainty, or
uncertainty from extrapolation beyond the data; poor fit quality can make the
band misleading.

## Simulation histograms

The separate `get_sim_efficiencies` executable requires explicit actions:

```sh
# Sort trees, then generate the ROOT graph file and CSV.
./EfficiencyToolkit/bin/get_sim_efficiencies --sort --curves --threads 4

# Or perform just one step.
./EfficiencyToolkit/bin/get_sim_efficiencies --sort
./EfficiencyToolkit/bin/get_sim_efficiencies --curves
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
newly generated flat file. Generate the flat file with `--debug` so its
individual graphs are available to the comparison tool:

```sh
make -C EfficiencyToolkit compare_layout
./EfficiencyToolkit/bin/compare_layout grouped.root flat.root
```
