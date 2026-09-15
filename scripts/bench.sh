#!/bin/sh
# Run the engine benchmark and write a timestamped record.
#
# POSIX shell. The script fails when cargo is not on PATH.
# The script writes one file to benchmarks/results/ and does not touch any
# other tree.

set -eu

if ! command -v cargo >/dev/null 2>&1; then
    echo "error: cargo not found on PATH." >&2
    echo "Install the Rust toolchain from https://rustup.rs and retry." >&2
    exit 127
fi

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
results_dir="$repo_root/benchmarks/results"
mkdir -p "$results_dir"

timestamp=$(date -u +%Y%m%dT%H%M%SZ)
iso_date=$(date -u +%Y-%m-%dT%H:%M:%SZ)
out_file="$results_dir/${timestamp}-engine-timings.txt"
raw_file="$results_dir/.${timestamp}.raw"

bench_args="--warm-up-time 0.5 --measurement-time 1.0 --sample-size 20"
command_line="cargo bench -p nanocad-engine --bench engine_bench -- $bench_args"

# Collect the host facts without failing when a probe is absent.
host_os=$(uname -s 2>/dev/null || echo unknown)
host_arch=$(uname -m 2>/dev/null || echo unknown)
host_name=$(hostname 2>/dev/null || echo unknown)
cpu_brand=unknown
cpu_cores=unknown
if [ "$host_os" = "Darwin" ]; then
    cpu_brand=$(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo unknown)
    cpu_cores=$(sysctl -n hw.logicalcpu 2>/dev/null || echo unknown)
elif [ -r /proc/cpuinfo ]; then
    cpu_brand=$(sed -n 's/^model name[[:space:]]*: //p' /proc/cpuinfo | head -n 1)
    [ -n "$cpu_brand" ] || cpu_brand=unknown
    cpu_cores=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo unknown)
else
    cpu_cores=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo unknown)
fi

cargo_version=$(cargo --version 2>/dev/null || echo unknown)
rustc_version=$(rustc --version 2>/dev/null || echo unknown)

# Run the benchmark. Keep the full output in a temporary raw file.
if ! (
    cd "$repo_root"
    cargo bench -p nanocad-engine --bench engine_bench -- $bench_args
) > "$raw_file" 2>&1; then
    cat "$raw_file" >&2
    rm -f "$raw_file"
    echo "error: cargo bench failed." >&2
    exit 1
fi
cat "$raw_file"

# Extract the "time:" summary lines. Criterion prints the benchmark name
# either on the same line or on the line before.
summary=$(
    awk '
        /Benchmarking/ { next }
        /time:/ {
            line = $0
            sub(/^[ \t]+/, "", line)
            if (line ~ /^time:/) {
                name = last_name
            } else {
                name = line
                sub(/ time:.*/, "", name)
            }
            t = line
            sub(/.*time:[ \t]*/, "", t)
            printf "| `%s` | %s |\n", name, t
            next
        }
        /^[^ \t]/ && !/change:/ && !/Found/ && !/Gnuplot/ && !/Running/ && !/Finished/ && !/warning/ && !/^#/ {
            last_name = $0
        }
    ' "$raw_file"
)
rm -f "$raw_file"

{
    echo "# Engine timing record (reproduced)"
    echo
    echo "Date (UTC): $iso_date"
    echo "Host name: $host_name"
    echo "Host: $cpu_brand, $cpu_cores logical cores"
    echo "OS: $host_os"
    echo "Architecture: $host_arch"
    echo "Toolchain: $cargo_version; $rustc_version"
    echo "Profile: \`bench\` (release, LTO thin, codegen-units 1)"
    echo "Command:"
    echo "\`$command_line\`"
    echo "Bench crate: \`crates/engine/benches/engine_bench.rs\`"
    echo
    echo "## Results"
    echo
    if [ -n "$summary" ]; then
        echo "| Benchmark | Criterion time |"
        echo "|---|---|"
        echo "$summary"
    else
        echo "No summary lines found in the benchmark output."
    fi
    echo
    echo "Criterion writes its raw report to \`target/criterion/\`."
} > "$out_file"

echo "wrote $out_file" >&2
