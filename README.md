# simple_draft_asm

Rust prototype for plant organelle dominant-form draft graph assembly.

The assembler keeps repeats and ambiguous branches in the graph. Plastid uses a
single dominant-form pass by default. Mitochondrial assembly defaults to two
rounds: build a conservative skeleton first, then remap reads to the skeleton to
rescue supported closing links.

## Build

```bash
cargo build --release
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

## Current validation

These were regenerated with the quick commands above:

| result | graph | S | L | bases | min | max |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `result_plastid_profile_default` | `graph.gfa` | 3 | 8 | 129,833 | 18,628 | 85,219 |
| `result_mito_profile_default` | `graph.gfa` | 19 | 46 | 376,128 | 741 | 61,105 |
| `result_mito_profile_default/round1_skeleton` | `graph.gfa` | 19 | 40 | 376,128 | 741 | 61,105 |
| `result_mito_profile_default/round1_readlinks` | `graph.gfa` | 19 | 64 | 376,128 | 741 | 61,105 |
| `result_mito_profile_default/round2_skeleton` | `skeleton.linked.gfa` | 19 | 46 | 376,128 | 741 | 61,105 |

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
