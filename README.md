# simple_draft_asm

Rust prototype for plant organelle dominant-form draft graph assembly.

The assembler keeps repeats and ambiguous branches in the graph. Plastid uses a
single dominant-form pass by default. Mitochondrial assembly defaults to two
rounds: build a conservative skeleton first, then remap reads to the skeleton to
rescue supported closing links.

## Build

```bash
/Users/zouyinstein-m4max/.cargo/bin/cargo build --release
```

If you want plain `cargo` to work in a new terminal, add Rust to `PATH` first:

```bash
export PATH="$HOME/.cargo/bin:$PATH"
```

## Quick runs

Plastid:

```bash
./target/release/simple_draft_asm --organelle plastid -i data/plastid.fastq.gz -o result_plastid_profile_default -t 8
```

Mitochondrion:

```bash
./target/release/simple_draft_asm --organelle mito -i data/mito.fastq.gz -o result_mito_profile_default -t 8
```

Use `--rounds 1` if you want only the first mitochondrial skeleton pass. Use
`--numt-interference low|high` to switch the mitochondrial profile; the mito
default is currently `high`.

## Compact Datasets

Use `--data-mode compact` or its alias `--small-dataset` for small corrected-read
inputs where standard high-depth profiles can drop low-support links before the
graph is complete. This mode keeps the standard rice plastid logic unchanged and
only applies when explicitly requested.

For the corrected MECAT plastid input:

```bash
./target/release/simple_draft_asm --organelle plastid \
  --small-dataset \
  -i data/mecat_corrected_plastid.fasta.gz \
  -o result_mecat_plastid_compact \
  -t 8
```

For the corrected MECAT mitochondrial input:

```bash
./target/release/simple_draft_asm --organelle mito \
  --small-dataset \
  -i data/mecat_corrected_mito.fasta.gz \
  -o result_mecat_mito_compact \
  -t 8
```

In mitochondrial compact mode only, the second round performs an additional
local bridge-completion step after read remapping. It identifies open skeleton
ends and small disconnected components, extracts reads touching just those
local regions, tests candidate bridges, and keeps only links that improve the
main topology without creating high-degree secondary branching. Unsupported
small components are omitted from the final main graph rather than force-linked.
This compact-mito completion path does not change plastid compact mode or the
standard large-data mito/plastid profiles.

## Link sensitivity

For large plastid datasets, weak secondary links can be retained even when the
dominant graph is already clear. Use `--min-link-ratio` to require each GFA link
to have support close to the best competing link at both of its endpoints:

```bash
./target/release/simple_draft_asm --organelle plastid \
  -i data/rice_plastid.fastq.gz \
  -o result_rice_plastid_clean \
  --min-link-ratio 0.30 \
  -t 8
```

The default is `0`, which preserves the previous fixed `--min-link-support`
behavior. A value such as `0.30` removes low-proportion secondary links while
keeping the primary endpoint-supported links.

For the current rice plastid dataset, use the cleaned 25% read subset as the
working parameter combination:

```bash
./target/release/simple_draft_asm --organelle plastid \
  -i data/rice_plastid.fastq.gz \
  -o result_rice_plastid_best \
  --min-link-ratio 0.30 \
  --subsets=25 \
  -t 8
```

For the current rice mitochondrial dataset, use the same 25% read subset and
relative link filtering with the default two-round mitochondrial workflow:

```bash
./target/release/simple_draft_asm --organelle mito \
  -i data/rice_mito.fastq.gz \
  -o result_rice_mito_best \
  --min-link-ratio 0.30 \
  --subsets=25 \
  -t 8
```

## Read subsampling

Use `--read-subsets` to run deterministic read-level subsampling experiments in
one command. `--subsets` is the same option with a shorter name. Reads are
selected before syncmer/minimizer discovery, and retained reads follow the
normal assembly path. For formal data-volume checks, use a halving series:

```bash
./target/release/simple_draft_asm --organelle plastid \
  -i data/rice_plastid.fastq.gz \
  -o result_rice_plastid_read_subsets \
  --subsets=12.5,25,50,100 \
  -t 8
```

Each subset is written under the output directory as `read_subset_25/`,
`read_subset_50/`, and so on. Decimal subsets use underscores in directory
names, for example `read_subset_12_5/`. The top-level `read_subsets.tsv`
records the elapsed time and read ID file for each subset. Non-100% subsets also
write `read_ids.txt` in their subset directory for downstream read extraction;
IDs are the first whitespace-delimited token in each FASTQ/FASTA header. The
default behavior is unchanged when neither `--read-subsets` nor `--subsets` is
supplied.

For workflows that need to read the data again, such as the default two-round
mitochondrial mode, non-100% subsets also materialize `reads.fasta` in the
subset directory. Round 1, skeleton remapping, and rescue all use that same
subset FASTA, so the subset is applied consistently across the full workflow.

## Main outputs

- `graph.gfa`: final draft graph.
- `unitigs.fasta`: final unitig sequences.
- `depth.tsv`: segment depth estimated by read remapping when available.
- `links.tsv`: junction/link support table.
- `report.txt`: run parameters and summary.

For two-round mitochondrial runs:

- `round1_skeleton/`: conservative first-pass graph.
- `round1_readlinks/`: first-pass graph with read-walk links for rescue evidence.
- `round2_skeleton/`: skeleton remapping and rescued final graph.
- top-level `graph.gfa`, `depth.tsv`, and `links.tsv` are copied from the final
  second-round result.

For mitochondrial compact runs, `round2_skeleton/` also includes local bridge
diagnostics:

- `mito_compact_bridge.report.txt`: bridge completion summary.
- `mito_compact_bridge.links.tsv`: accepted local bridge candidates.
- `mito_compact_bridge.pruned.tsv`: secondary links and small components
  removed from the final main topology.
- `mito_compact_bridge.read_ids.txt` and `mito_compact_bridge.reads.fasta`:
  reads selected for local bridge inspection.

## Current validation

These were regenerated with the quick and compact commands above:

| result | graph | S | L | bases | min | max |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `result_plastid_profile_default` | `graph.gfa` | 3 | 8 | 129,833 | 18,628 | 85,219 |
| `result_mito_profile_default` | `graph.gfa` | 19 | 46 | 376,128 | 741 | 61,105 |
| `result_mito_profile_default/round1_skeleton` | `graph.gfa` | 19 | 40 | 376,128 | 741 | 61,105 |
| `result_mito_profile_default/round1_readlinks` | `graph.gfa` | 19 | 64 | 376,128 | 741 | 61,105 |
| `result_mito_profile_default/round2_skeleton` | `skeleton.linked.gfa` | 19 | 46 | 376,128 | 741 | 61,105 |
| `result_mecat_plastid_compact` | `graph.gfa` | 3 | 8 | 115,020 | 13,345 | 81,135 |
| `result_mecat_mito_compact` | `graph.gfa` | 15 | 34 | 314,441 | 2,072 | 41,260 |
| `codex_mecat_mito_AT007_compact` | `graph.gfa` | 18 | 40 | 365,716 | 917 | 76,599 |
| `codex_mecat_mito_AT008_compact` | `graph.gfa` | 17 | 38 | 363,552 | 699 | 62,428 |

Compact mitochondrial topology note: the MECAT 601-read input does not contain
validated bridge evidence for every small component/end. The compact-mito
completion step produces a single main component with endpoint degree bounded at
2, matching the large rice mitochondrial graph shape at the S/L topology level
without requiring one-to-one segment correspondence. The Arabidopsis compact
checks above close as single main components with no open endpoints; AT007 uses
one focused local PAF bridge, while AT008 closes from skeleton plus full-read PAF
evidence without an extra bridge.

## Benchmark

Benchmark date: 2026-06-11/12. Commands were run on this machine:

- Model: Mac Studio `Mac16,9`
- Chip: Apple M4 Max, 14 cores (10 performance + 4 efficiency)
- Memory: 36 GB
- OS: macOS 26.5.1 (25F80)
- Threads: `-t 8`
- Inputs: `data/mito.fastq.gz` 73 MB, `data/plastid.fastq.gz` 235 MB
- Tool versions: Flye `2.9.6-b1802`, minimap2 `2.30-r1287`,
  OATK/syncasm `1.0`, simple_draft_asm `0.1.0`

Runtime was measured with `/usr/bin/time -p`. OATK was installed from
`c-zhou/oatk` under `external/oatk` and compiled with `make -j8`. The benchmark
uses OATK's core graph assembler `syncasm`, not the full `oatk` wrapper with HMM
annotation/pathfinder.

| tool | dataset | command profile | real | user | sys | output graph | S | L | bases |
| --- | --- | --- | ---: | ---: | ---: | --- | ---: | ---: | ---: |
| OATK/syncasm | mito | `-k 1001 -c 30` | 0.99s | 1.01s | 0.02s | `benchmarks/oatk_mito/oatk_mito.utg.final.gfa` | 9 | 24 | 364,639 |
| OATK/syncasm | plastid | `-k 1001 -c 30` | 2.74s | 3.38s | 0.08s | `benchmarks/oatk_plastid/oatk_plastid.utg.final.gfa` | 3 | 8 | 132,513 |
| simple_draft_asm | mito | `--rounds 1` | 1.65s | 1.25s | 0.16s | `benchmarks/simple_mito_1round/graph.gfa` | 19 | 36 | 376,128 |
| simple_draft_asm | mito | default mito, 2 rounds | 2.18s | 5.89s | 0.31s | `benchmarks/simple_mito/graph.gfa` | 19 | 46 | 376,128 |
| simple_draft_asm | plastid | default plastid, 1 round | 3.08s | 3.60s | 0.19s | `benchmarks/simple_plastid/graph.gfa` | 3 | 8 | 129,833 |
| Flye | mito | `--genome-size 500k` | 226.27s | 700.85s | 2.22s | `benchmarks/flye_mito_full/assembly_graph.gfa` | 12 | 17 | 370,282 |
| Flye | plastid | `--genome-size 160k` | >10m50s | n/a | n/a | aborted during read extension | n/a | n/a | n/a |

Speed summary:

- Mito: OATK/syncasm was about 1.7x faster than simple_draft_asm `--rounds 1`
  and about 2.2x faster than simple_draft_asm's default two-round run.
- Mito: Flye was about 103.8x slower than simple_draft_asm's default two-round
  run and about 228.6x slower than OATK/syncasm.
- Plastid: OATK/syncasm was about 1.1x faster than simple_draft_asm.
- Flye plastid was not completed because the high-coverage run was still in read
  extension after more than 10 minutes.
- The default simple_draft_asm mito two-round path now reuses the first-round
  anchor assembly when writing the readlink rescue graph, avoiding a duplicate
  first-round assembly pass.
- The simple_draft_asm anchor pass avoids repeated k-mer/suffix string
  allocation in hot loops; this keeps graph output unchanged while reducing
  runtime.

Benchmark commands:

```bash
/usr/bin/time -p ./target/release/simple_draft_asm --organelle mito -i data/mito.fastq.gz -o benchmarks/simple_mito -t 8
/usr/bin/time -p ./target/release/simple_draft_asm --organelle mito --rounds 1 -i data/mito.fastq.gz -o benchmarks/simple_mito_1round -t 8
/usr/bin/time -p ./target/release/simple_draft_asm --organelle plastid -i data/plastid.fastq.gz -o benchmarks/simple_plastid -t 8

/usr/bin/time -p external/oatk/syncasm -k 1001 -c 30 -t 8 -o benchmarks/oatk_mito/oatk_mito data/mito.fastq.gz
/usr/bin/time -p external/oatk/syncasm -k 1001 -c 30 -t 8 -o benchmarks/oatk_plastid/oatk_plastid data/plastid.fastq.gz

/usr/bin/time -p /opt/homebrew/bin/flye --pacbio-hifi data/mito.fastq.gz --extra-params output_gfa_before_rr=1 --genome-size 500k -t 8 -o benchmarks/flye_mito_full
/usr/bin/time -p /opt/homebrew/bin/flye --pacbio-hifi data/plastid.fastq.gz --extra-params output_gfa_before_rr=1 --genome-size 160k -t 8 -o benchmarks/flye_plastid_full
```

## Help

Common options are shown with:

```bash
./target/release/simple_draft_asm --help
```

Low-level tuning parameters are hidden from the normal interface but documented
with:

```bash
./target/release/simple_draft_asm --help-advanced
```
