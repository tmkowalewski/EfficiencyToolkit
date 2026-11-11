#!/usr/bin/env bash
set -euo pipefail

toolkit=$(cd "$(dirname "$0")/.." && pwd)
cd "$toolkit/.."
test_dir=$(mktemp -d)
trap 'rm -f "$test_dir/relative.root" "$test_dir/absolute.root" "$test_dir/fit.root" "$test_dir/fixed.root" "$test_dir/normal.root" "$test_dir/normalized.root" "$test_dir/multiple.root" "$test_dir/example.root" "$test_dir/invalid.csv" "$test_dir/no-overlap.csv" "$test_dir/few.csv" "$test_dir/protected.root" "$test_dir/log"; rmdir "$test_dir"' EXIT

"$toolkit/bin/make_eff_curves" --debug --output "$test_dir/relative.root"
"$toolkit/bin/make_eff_curves" --debug --source-info "$toolkit/tests/source_exposure_test.csv" --output "$test_dir/absolute.root"
"$toolkit/bin/make_eff_curves" --debug --radware-fit --source-info "$toolkit/tests/source_exposure_test.csv" --output "$test_dir/fit.root" >"$test_dir/log" 2>&1
grep -q 'Warning:.*RadWare fit has chi2/ndf' "$test_dir/log"
"$toolkit/bin/make_eff_curves" --radware-fit --radware-b 1 --source-info "$toolkit/tests/source_exposure_test.csv" --output "$test_dir/fixed.root" >"$test_dir/log" 2>&1
"$toolkit/bin/check_radware" "$test_dir/fit.root" "$test_dir/fixed.root"
"$toolkit/bin/make_eff_curves" --radware-fit --source-info "$toolkit/tests/source_exposure_test.csv" --output "$test_dir/normal.root" >"$test_dir/log" 2>&1
"$toolkit/bin/check_plot_output" "$test_dir/normal.root" "$test_dir/fit.root"
sed 's/^60Co,.*/60Co,ignored,,not-a-number,-1,,normalize/' "$toolkit/tests/source_exposure_test.csv" >"$test_dir/invalid.csv"
"$toolkit/bin/make_eff_curves" --debug --source-info "$test_dir/invalid.csv" --output "$test_dir/normalized.root"
"$toolkit/bin/check_absolute" "$test_dir/relative.root" "$test_dir/absolute.root" "$test_dir/normalized.root"
sed 's/,measured$/,example/' "$toolkit/tests/source_exposure_test.csv" >"$test_dir/invalid.csv"
"$toolkit/bin/make_eff_curves" --source-info "$test_dir/invalid.csv" --output "$test_dir/example.root" >"$test_dir/log" 2>&1
grep -q 'Warning: example exposure values used for 152Eu' "$test_dir/log"

header='source,activity_bq,activity_uncertainty_bq,elapsed_time_s,live_fraction,live_fraction_uncertainty,status'
eu='152Eu,1000000,50000,40000,0.8,0.04,measured'
co56='56Co,2000000,20000,30000,0.9,0.01,measured'
co60='60Co,500000,25000,10000,0.95,0.01,measured'
printf '%s\n%s\n%s\n%s\n' "$header" '152Eu,,,,,,normalize' "$co56" '60Co,,,,,,normalize' >"$test_dir/invalid.csv"
"$toolkit/bin/make_eff_curves" --source-info "$test_dir/invalid.csv" --output "$test_dir/multiple.root"
printf 'original output\n' >"$test_dir/protected.root"
for option in '--radware-g 0' '--radware-scale -1' '--radware-b nan'; do
    if printf 'y\n' | "$toolkit/bin/make_eff_curves" --radware-fit $option --output "$test_dir/protected.root" >"$test_dir/log" 2>&1; then
        echo "Unexpected success for invalid option $option" >&2
        exit 1
    fi
    grep -qx 'original output' "$test_dir/protected.root"
done
if "$toolkit/bin/make_eff_curves" --radware-b 1 --output "$test_dir/protected.root" >"$test_dir/log" 2>&1; then
    echo "Unexpected success for RadWare parameter without fit" >&2
    exit 1
fi
grep -q 'RadWare parameters require --radware-fit' "$test_dir/log"
head -n 5 data/156Eu-56Co-60Co_RawCounts.csv >"$test_dir/few.csv"
if printf 'y\n' | "$toolkit/bin/make_eff_curves" --raw "$test_dir/few.csv" --radware-fit --output "$test_dir/protected.root" >"$test_dir/log" 2>&1; then
    echo "Unexpected success for too few RadWare points" >&2
    exit 1
fi
grep -q 'RadWare fit requires at least 6 points' "$test_dir/log"
grep -qx 'original output' "$test_dir/protected.root"

check_failure() {
    local label=$1 expected=$2
    if printf 'y\n' | "$toolkit/bin/make_eff_curves" --source-info "$test_dir/invalid.csv" --output "$test_dir/protected.root" >"$test_dir/log" 2>&1; then
        echo "Unexpected success for $label" >&2
        exit 1
    fi
    if ! grep -q "$expected" "$test_dir/log"; then
        echo "Unexpected error for $label:" >&2
        cat "$test_dir/log" >&2
        exit 1
    fi
    if [ "$(cat "$test_dir/protected.root")" != 'original output' ]; then
        echo "Existing output changed after $label" >&2
        exit 1
    fi
}

printf '%s\n%s\n%s\n' "$header" "$eu" "$co56" >"$test_dir/invalid.csv"
check_failure "missing source" 'Missing source exposure'
printf '%s\n%s\n%s\n%s\n%s\n' "$header" "$eu" "$co56" "$co60" "$eu" >"$test_dir/invalid.csv"
check_failure "duplicate source" 'Duplicate source exposure'
printf '%s\n%s\n%s\n%s\n' "$header" '152Eu,,,,,,normalize' '56Co,,,,,,normalize' '60Co,,,,,,normalize' >"$test_dir/invalid.csv"
check_failure "no measured reference" 'No usable overlap with measured sources'
printf '%s\n%s\n%s\n%s\n' "$header" "${eu%,measured},example" "${co56%,measured},example" '60Co,,,,,,normalize' >"$test_dir/invalid.csv"
check_failure "example is not a measured reference" 'No usable overlap with measured sources'
printf '%s\n%s\n%s\n%s\n' "$header" "$eu" "$co56" '60Co,,,,,,normalize' >"$test_dir/invalid.csv"
awk -F, 'BEGIN { OFS="," } $7 == "60Co" { $1 += 50000 } { print }' data/156Eu-56Co-60Co_RawCounts.csv >"$test_dir/no-overlap.csv"
if printf 'y\n' | "$toolkit/bin/make_eff_curves" --raw "$test_dir/no-overlap.csv" --source-info "$test_dir/invalid.csv" --output "$test_dir/protected.root" >"$test_dir/log" 2>&1; then
    echo "Unexpected success for non-overlapping measured references" >&2
    exit 1
fi
grep -q 'No usable overlap with measured sources' "$test_dir/log"
grep -qx 'original output' "$test_dir/protected.root"
for invalid in \
    '152Eu,not-a-number,50000,40000,0.8,0.04,measured' \
    '152Eu,0,50000,40000,0.8,0.04,measured' \
    '152Eu,-1,50000,40000,0.8,0.04,measured' \
    '152Eu,1000000,-1,40000,0.8,0.04,measured' \
    '152Eu,1000000,50000,0,0.8,0.04,measured' \
    '152Eu,1000000,50000,40000,0,0.04,measured' \
    '152Eu,1000000,50000,40000,1.1,0.04,measured' \
    '152Eu,1000000,50000,40000,0.8,-1,measured' \
    '152Eu,1000000,50000,40000,0.8,0.04,unknown'; do
    printf '%s\n%s\n%s\n%s\n' "$header" "$invalid" "$co56" "$co60" >"$test_dir/invalid.csv"
    check_failure "$invalid" 'Efficiency plots:'
done
printf '%s\n%s\n' "$header" '152Eu,1000000' >"$test_dir/invalid.csv"
check_failure "malformed row" 'Wrong number of fields'
rm -f "$test_dir/invalid.csv"
check_failure "missing file" 'Cannot open source exposure CSV'

printf '%s\n' "$header" "$eu" "$co56" "$co60" >"$test_dir/invalid.csv"
if "$toolkit/bin/make_eff_curves" --source-info "$test_dir/invalid.csv" --output "$test_dir/invalid.csv" >"$test_dir/log" 2>&1; then
    echo "Unexpectedly allowed input overwrite" >&2
    exit 1
fi
grep -q 'Output must not overwrite an input file' "$test_dir/log"
grep -qx "$eu" "$test_dir/invalid.csv"
echo 'Exposure validation checks passed'
