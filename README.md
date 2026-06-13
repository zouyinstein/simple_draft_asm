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
./target/release/simple_draft_asm -p ps -i data/plastid.fastq.gz -o result_plastid_profile_default -t 8
```

Mitochondrion:

```bash
./target/release/simple_draft_asm -p ms -i data/mito.fastq.gz -o result_mito_profile_default -t 8
```

Use `--rounds 1` if you want only the first mitochondrial skeleton pass. Use
`--numt-interference low|high` to switch the mitochondrial profile; the mito
default is currently `high`.

## Presets

Use `--preset` or `-p` for the common organelle/profile combinations:

| preset | aliases | meaning |
| --- | --- | --- |
| `ml` | `mito_low` | mitochondrial low-data mode, replacing the user-facing compact command |
| `ms` | `mito_standard` | mitochondrial standard mode, including the default two-round skeleton workflow |
| `mh` | `mito_high` | mitochondrial high-data mode: standard mode plus `--min-link-ratio 0.30 --subsets=25,50,100` |
| `mx` | `mito_stable`, `mito_complex` | full mitochondrial complex-repeat mode with internal-junction splitting and stable bridge selection |
| `pl` | `plastid_low` | plastid low-data mode, replacing the user-facing compact command |
| `ps` | `plastid_standard` | plastid standard mode, keeping the plastid one-round workflow |
| `ph` | `plastid_high` | plastid high-data mode: standard mode plus `--min-link-ratio 0.30 --subsets=25,50,100` |

Presets are shorthand over the existing options; the older `--organelle`,
`--data-mode compact`, `--small-dataset`, `--min-link-ratio`, and `--subsets`
flags are still supported. Later explicit options can still override a preset.
The preset layer does not merge plastid and mito internals: plastid keeps its
one-round graph logic, while mito keeps its two-round skeleton/remapping logic
unless `--rounds` is explicitly changed.

## Compact Datasets

Use low presets (`-p ml` or `-p pl`) for small corrected-read inputs where
standard high-depth profiles can drop low-support links before the graph is
complete. The older `--data-mode compact` and `--small-dataset` flags are kept
as aliases for compatibility. This mode keeps the standard large-data plastid
and mito profiles unchanged and only applies when explicitly requested.

For the corrected MECAT plastid input:

```bash
./target/release/simple_draft_asm -p pl \
  -i data/mecat_corrected_plastid.fasta.gz \
  -o result_mecat_plastid_compact \
  -t 8
```

For the corrected MECAT mitochondrial input:

```bash
./target/release/simple_draft_asm -p ml \
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

## Operating Logic

Start with standard mode unless the data volume already tells you otherwise.
Plastid standard mode is a one-round dominant graph workflow. Mitochondrial
standard mode is two-round: the first round builds a conservative skeleton and
the second round remaps reads to that skeleton to extend/rescue supported open
graph structure.

Use compact mode only for genuinely small corrected-read inputs. Randomly
sampling 500 or 1000 raw reads from a large mitochondrial dataset is not a
replacement for the full standard run; on `data/Col-0_mito.fastq.gz`, 500-read
and 1000-read compact tests did not match the full standard topology.

Use high-data mode for very large datasets when weak secondary links need to be
controlled by read subsetting and link filtering:

```bash
./target/release/simple_draft_asm -p mh \
  -i data/mito.fastq.gz \
  -o result_mito_high \
  -t 8
```

When a mitochondrial standard graph has suspicious `1-2` or `2-1` nodes, use
`mx` to analyze or repair only the selected repeat neighborhood rather than
continuing to tune global thresholds.

## Experimental Mito Stable Mode

Use `-p mx` for mitochondrial datasets where the normal standard graph may be
closed but repeat topology is unstable. The recommended workflow is manual: run
standard `ms` first, inspect the graph, then rerun `mx` when the standard graph
has suspicious repeat-node structure or sample-to-sample instability.

```bash
./target/release/simple_draft_asm -p ms \
  -i data/CRR958891_mito.fastq.gz \
  -o result_crr_mito_ms \
  -t 8

./target/release/simple_draft_asm -p mx \
  -i data/CRR958891_mito.fastq.gz \
  -o result_crr_mito_mx \
  -t 8
```

The `mx` rerun is not a global threshold sweep. It changes the topology search
itself:

1. Build the normal mitochondrial two-round skeleton, including the first-pass
   read-walk rescue graph.
2. Remap all reads to the skeleton and look for strong internal split points
   inside existing segments, not only at segment ends.
3. Split skeleton segments at supported internal junctions and remap reads when
   new breakpoints are found.
4. Build a candidate graph from existing skeleton links, PAF links, local
   read-supported bridges, and rescue-graph links.
5. For paired `1-2`/`2-1` nodes, test whether their single sides are already
   connected through a short `2-2` bridge-chain. When the bridge-chain support
   strongly dominates single-copy reconnect alternatives, add the direct
   double-copy chord.
6. Score candidate bridge sets globally by connectedness, open ends, endpoint
   overload, branch count, cycle rank, link count, and support.
7. Keep only bridges that improve the main topology without creating excessive
   secondary links.
8. Apply repeat-aware shortcut replacement when a short local chain can explain
   an otherwise ambiguous connection while preserving valid node degrees.
9. Prune closed redundant links only after confirming the graph remains
   connected, closed, and within degree guards.

`mx` has three user-facing modes:

- `--mx-mode auto`: full automatic mito-stable search. This can add bridges,
  perform repeat-aware shortcut replacement, and prune redundant closed links.
- `--mx-mode forbid-selected --selected-nodes LIST`: only analyze the selected
  nodes and try to repair those nodes with local candidate bridges. Automatic
  repeat expansion and pruning are disabled in this mode.
- `--mx-mode allow-selected --selected-nodes LIST`: analyze the selected nodes
  without changing the graph. This is useful when deciding whether selected
  `1-2`/`2-1` nodes should remain as repeat-copy structure.

For example, to test whether `edge_12` and `edge_17` should be repaired:

```bash
./target/release/simple_draft_asm -p mx \
  --skeleton-gfa data/check/all_mito_500K_before_rr.gfa \
  -i data/CRR958891_mito.fastq.gz \
  -o benchmarks/crr_mx_repair_20260613/forbid_edge12_edge17 \
  -t 8 \
  --mx-mode forbid-selected \
  --selected-nodes edge_12,edge_17
```

To inspect whether `edge_4` and `edge_8` are better interpreted as retained
repeat-copy nodes:

```bash
./target/release/simple_draft_asm -p mx \
  --skeleton-gfa data/check/all_mito_500K_before_rr.gfa \
  -i data/CRR958891_mito.fastq.gz \
  -o benchmarks/crr_mx_repair_20260613/allow_edge4_edge8 \
  -t 8 \
  --mx-mode allow-selected \
  --selected-nodes edge_4,edge_8
```

For unstable mitochondrial repeat data, this is meant to separate stable
regions from unstable repeat neighborhoods and then solve those neighborhoods
with local split/bridge/repeat repair. The topology scan treats node names as
graph-local identifiers; compare roles, degrees, neighbors, and sequence
placement between runs rather than assuming names are stable.

The main `mx` diagnostics are:

- `mito_stable_splits.tsv`: internal breakpoints used to expose hidden branch
  endpoints.
- `mito_stable_bridge.candidates.tsv` and `mito_stable_bridge.selected.tsv`:
  tested and accepted bridge links.
- `mito_stable_bridge.report.txt`: topology score summary for the selected
  bridge set.
- `mito_stable_selected_nodes.tsv`: selected node degree, incident base links,
  and incident read-supported candidates.
- `mito_stable_copy_choice.tsv`: local single-copy versus double-copy
  interpretations for three-way node pairs, including added-link support,
  support source such as `read_walk` or `bridge_chain:*`, and topology outcome.
- `mito_stable_node_degrees.tsv`: physical left/right degree class per final
  node.
- `mito_stable_node_repairs.tsv`: before/after degree class and repair events
  for nodes changed by selected bridges, manual edits, pruning, or repeat
  expansion.
- `mito_stable_topology_scan.tsv`: final node-role scan, including accepted
  three-way merge candidates and invalid node shapes.
- `mito_stable_manual_edits.tsv`: advanced forced/dropped link edits, when
  supplied.
- `mito_stable_repeat_expansions.tsv`: repeat-aware shortcut replacements.
- `mito_stable_pruned.tsv`: redundant links removed after topology validation.

If you already have an incomplete or suspicious GFA and want to repair that
specific graph, add `--skeleton-gfa`; this is an input graph to inspect and
repair, not a topology reference:

```bash
./target/release/simple_draft_asm -p mx \
  --skeleton-gfa result_crr_mito_ms/graph.gfa \
  -i data/CRR958891_mito.fastq.gz \
  -o result_crr_mito_mx_from_gfa \
  -t 8
```

Advanced `--force-link` and `--drop-link` options remain available as local
graph-edit controls. They use `from:from_orient:to:to_orient` syntax, for
example `edge_8:-:edge_4:-`, and are recorded in
`mito_stable_manual_edits.tsv`.

The manual two-step repeat audit remains useful when reviewing a suspicious GFA:
first run `forbid-selected` on `edge_12,edge_17`, then run `forbid-selected` on
`edge_4,edge_8` from that repaired graph. The second run only lists
`edge_8:-:edge_4:-` as selected because `edge_17:-:edge_12:-` is already
present in its input graph; the two steps together make both repeat cassettes
`2-2`.

The current direct CRR reruns show the bridge-chain behavior from direct `ms`
graphs. `mx auto` resolves CRR958891 by adding `utg12:+:utg14:+` through the
`utg13` bridge-chain and `utg17:+:utg20:+` through the `utg21` bridge-chain.
The same motif also triggers on CRR958893, adding `utg12:-:utg14:+` through
`utg13` and `utg0:+:utg18:+` through `utg19`. See each run's
`mito_stable_copy_choice.tsv` for support and source.

For CRR958891, use
`CRR958891_mx_auto_from_cleaned_ms/graph.gfa`
as the final graph. The cleaned input drops `utg10` because its `ms` depth is 0,
and drops the manually rejected `utg8--utg26` direct link before `mx auto`.
This keeps sample-level cleanup separate from `mx`: `mx auto` then only adds
the two repeat-cassette chords `utg12:+:utg14:+` and `utg17:+:utg20:+`.

Why the `utg17--utg20` direct link is easy to lose in `ms`: the evidence is
not absent, but the standard graph representation assigns it to the intervening
short repeat bridge `utg21`. In the CRR958891 `ms` graph, `utg21` is 1,520 bp
and the retained links are very strong:

```text
utg17:+ -> utg21:-  skeleton 146, PAF 117
utg21:- -> utg20:+  skeleton 140, PAF 127
```

The direct chord `utg17:+ -> utg20:+` only appears as a weak raw candidate in
`links.tsv` with PAF support 3 and link ratios about 0.02, so it does not pass
the normal `ms` graph filters and is absent from `ms/graph.gfa`. Biologically
and topologically, however, this is the same double-copy cassette evidence:
`utg17:R--utg21:R` plus `utg21:L--utg20:L` says the two single sides are joined
through a short `2-2` bridge. `mx auto` therefore converts that bridge-chain
support into the direct double-copy chord, writing
`utg17:+:utg20:+` with support 140 and `support_source=bridge_chain:utg21`.

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

For experimental mito-stable runs, `round2_skeleton/` and the top-level output
also include:

- `mito_stable_splits.tsv`: high-support internal skeleton breakpoints used to
  expose hidden branch endpoints before final bridge selection.
- `mito_stable_bridge.report.txt`: global bridge-selection topology summary.
- `mito_stable_bridge.candidates.tsv` and `mito_stable_bridge.selected.tsv`:
  candidate and accepted bridge links after internal splitting.
- `mito_stable_node_degrees.tsv`: physical left/right degree class for each
  final node.
- `mito_stable_topology_scan.tsv`: self-check of node roles; highlights
  accepted three-way merge candidates through small `2-2` repeat bridges,
  ordinary `2-2` repeat nodes, linear nodes, and invalid node shapes.
- `mito_stable_repeat_expansions.tsv`: copied repeat segments used to replace
  low-overlap shortcut links with local high-overlap chains; this file can be
  empty when node-degree constraints reject the repair.
- `mito_stable_pruned.tsv`: closed redundant links removed after repeat repair.

## Current validation

The latest validation set was regenerated under
`benchmarks/simple_refresh_20260612/` and
`benchmarks/tool_refresh_20260612/` on 2026-06-12. The simple refresh keeps
Col-0 large datasets in standard mode only, rice large datasets in standard and
`--min-link-ratio 0.30 --subsets=25,100` modes, and the two Arabidopsis MECAT
mitochondrial datasets in compact mode. OATK was rerun for all six inputs. Flye
was rerun only for Col-0 mitochondrial and the two compact mitochondrial inputs;
rice Flye and Col-0 plastid Flye are skipped because those inputs are too large
for the current comparison budget. For future Flye comparison rows, runs that
exceed 10 minutes should be recorded as `>10m` with `n/a` graph statistics.

The Arabidopsis mitochondrial compact datasets are the positive regression
tests for the compact-mito bridge workflow. They should keep the compact-mito
specific local bridge behavior while leaving plastid compact and large-data
profiles unchanged.

| input | mode | output graph | S | L | bases | components | open endpoints | bridge behavior |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| `data/mecat_mito_Arb-0.fasta.gz` | compact mito | `benchmarks/simple_refresh_20260612/results/mecat_mito_Arb-0_compact/graph.gfa` | 18 | 40 | 365,716 | 1 | 0 | 179 focused reads, 1 accepted local PAF bridge |
| `data/mecat_mito_AUZE-A-5.fasta.gz` | compact mito | `benchmarks/simple_refresh_20260612/results/mecat_mito_AUZE-A-5_compact/graph.gfa` | 17 | 38 | 363,552 | 1 | 0 | already closed by skeleton plus full-read PAF |

## Benchmark

Benchmark date: 2026-06-12. Commands were run on this machine:

- Model: Mac Studio `Mac16,9`
- Chip: Apple M4 Max, 14 cores (10 performance + 4 efficiency)
- Memory: 36 GB
- OS: macOS 26.5.1 (25F80)
- Threads: `-t 8`
- Tool versions: Flye `2.9.6-b1802`, minimap2 `2.30-r1287`,
  OATK/syncasm `1.0`, simple_draft_asm `0.1.0`

Standard and compact simple rows use `/usr/bin/time -p` wall time. Subset rows
use the per-subset `elapsed_seconds` values from `read_subsets.tsv`. OATK and
Flye rows use `/usr/bin/time -p` wall time from
`benchmarks/tool_refresh_20260612/logs/`.

| tool | dataset | profile | elapsed | graph | S | L | bases | components |
| --- | --- | --- | ---: | --- | ---: | ---: | ---: | ---: |
| simple_draft_asm | Col-0 mito | standard | 2.33s | `benchmarks/simple_refresh_20260612/results/col0_mito_standard/graph.gfa` | 19 | 46 | 376,128 | 1 |
| OATK/syncasm | Col-0 mito | `-k 1001 -c 30` | 0.98s | `benchmarks/tool_refresh_20260612/oatk/col0_mito/col0_mito.utg.final.gfa` | 9 | 24 | 364,639 | 1 |
| Flye | Col-0 mito | `--genome-size 500k` | 225.64s | `benchmarks/tool_refresh_20260612/flye/col0_mito/assembly_graph.gfa` | 12 | 17 | 370,282 | 1 |
| simple_draft_asm | Col-0 plastid | standard | 3.63s | `benchmarks/simple_refresh_20260612/results/col0_plastid_standard/graph.gfa` | 3 | 8 | 129,833 | 1 |
| OATK/syncasm | Col-0 plastid | `-k 1001 -c 30` | 2.77s | `benchmarks/tool_refresh_20260612/oatk/col0_plastid/col0_plastid.utg.final.gfa` | 3 | 8 | 132,513 | 1 |
| Flye | Col-0 plastid | skipped, too slow | n/a | n/a | n/a | n/a | n/a | n/a |
| simple_draft_asm | rice mito | standard | 17.60s | `benchmarks/simple_refresh_20260612/results/rice_mito_standard/graph.gfa` | 151 | 179 | 1,127,187 | 75 |
| simple_draft_asm | rice mito | `--min-link-ratio 0.30 --subsets=25` | 3.928s | `benchmarks/simple_refresh_20260612/results/rice_mito_subsets_ratio030/read_subset_25/graph.gfa` | 18 | 44 | 364,738 | 1 |
| simple_draft_asm | rice mito | `--min-link-ratio 0.30 --subsets=100` | 27.039s | `benchmarks/simple_refresh_20260612/results/rice_mito_subsets_ratio030/read_subset_100/graph.gfa` | 151 | 153 | 1,127,187 | 84 |
| OATK/syncasm | rice mito | `-k 1001 -c 30` | 4.19s | `benchmarks/tool_refresh_20260612/oatk/rice_mito/rice_mito.utg.final.gfa` | 263 | 442 | 2,877,837 | 110 |
| Flye | rice mito | skipped, dataset too large | n/a | n/a | n/a | n/a | n/a | n/a |
| simple_draft_asm | rice plastid | standard | 11.40s | `benchmarks/simple_refresh_20260612/results/rice_plastid_standard/graph.gfa` | 5 | 29 | 115,873 | 1 |
| simple_draft_asm | rice plastid | `--min-link-ratio 0.30 --subsets=25` | 18.103s | `benchmarks/simple_refresh_20260612/results/rice_plastid_subsets_ratio030/read_subset_25/graph.gfa` | 3 | 8 | 115,266 | 1 |
| simple_draft_asm | rice plastid | `--min-link-ratio 0.30 --subsets=100` | 18.842s | `benchmarks/simple_refresh_20260612/results/rice_plastid_subsets_ratio030/read_subset_100/graph.gfa` | 5 | 14 | 115,873 | 1 |
| OATK/syncasm | rice plastid | `-k 1001 -c 30` | 7.03s | `benchmarks/tool_refresh_20260612/oatk/rice_plastid/rice_plastid.utg.final.gfa` | 263 | 326 | 3,505,084 | 150 |
| Flye | rice plastid | skipped, dataset too large | n/a | n/a | n/a | n/a | n/a | n/a |
| simple_draft_asm | mecat mito Arb-0 | compact mito | 0.65s | `benchmarks/simple_refresh_20260612/results/mecat_mito_Arb-0_compact/graph.gfa` | 18 | 40 | 365,716 | 1 |
| OATK/syncasm | mecat mito Arb-0 | `-k 1001 -c 30` | 0.14s | `benchmarks/tool_refresh_20260612/oatk/mecat_mito_Arb-0/mecat_mito_Arb-0.utg.final.gfa` | 4 | 0 | 286,320 | 4 |
| Flye | mecat mito Arb-0 | `--genome-size 500k` | 35.46s | `benchmarks/tool_refresh_20260612/flye/mecat_mito_Arb-0/assembly_graph.gfa` | 2 | 2 | 368,875 | 2 |
| simple_draft_asm | mecat mito AUZE-A-5 | compact mito | 0.67s | `benchmarks/simple_refresh_20260612/results/mecat_mito_AUZE-A-5_compact/graph.gfa` | 17 | 38 | 363,552 | 1 |
| OATK/syncasm | mecat mito AUZE-A-5 | `-k 1001 -c 30` | 0.14s | `benchmarks/tool_refresh_20260612/oatk/mecat_mito_AUZE-A-5/mecat_mito_AUZE-A-5.utg.final.gfa` | 4 | 0 | 301,896 | 4 |
| Flye | mecat mito AUZE-A-5 | `--genome-size 500k` | 37.08s | `benchmarks/tool_refresh_20260612/flye/mecat_mito_AUZE-A-5/assembly_graph.gfa` | 3 | 4 | 364,696 | 1 |

Benchmark command set:

```bash
bash benchmarks/simple_refresh_20260612/run_simple_benchmarks.sh
bash benchmarks/tool_refresh_20260612/run_oatk_benchmarks.sh
bash benchmarks/tool_refresh_20260612/run_flye_remaining_benchmarks.sh
```

OATK uses `external/oatk/syncasm -k 1001 -c 30 -t 8`. Flye uses
`--pacbio-hifi`, `--extra-params output_gfa_before_rr=1`, and
`--genome-size 500k` for mitochondrial inputs or `--genome-size 160k` for
plastid inputs. Flye rows that are intentionally skipped or still active after
10 minutes are recorded as `n/a`.

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
