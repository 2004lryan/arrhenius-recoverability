# Datasets — sources, licenses, and availability

This repository ships **data cards and derived artefacts only** (see `../docs/`),
not the raw source data. Obtain each dataset from its official source below.
SHA-256 checksums are the full 64-hex digests recorded in `../docs/datacard-c5.md`;
they are reproduced by `python code/120datacard_sha256.py`.

No new data were collected for this work. All three real systems are
previously published, openly licensed third-party data.

## External public datasets (download from the official source; not re-hosted)

| Dataset | Role in the paper | Source | License | SHA-256 (16) |
|---|---|---|---|---|
| Hass avocado postharvest ripening (photographic) | Real system 3 — the only system with a genuine sampling distribution (478 fruit, 14,722 records) | [Mendeley Data, doi:10.17632/3xd9n945v8.1](https://doi.org/10.17632/3xd9n945v8.1) | CC BY 4.0 | `f8abeaba6eedf678` |
| Apple firmness temperature gradient (Table A1) | Real system 1 — 4 cultivars × 5 storage temperatures × 5 time points; design-window sweep | [Foods 2023, doi:10.3390/foods12112113](https://doi.org/10.3390/foods12112113) | CC BY 4.0 | `d22131f037a120e7` |
| Long-shelf-life food kinetics, ammonia / TBARS (Table 3) | Real system 2 — 5 foods × 4 storage temperatures, 10 sequences; independent replication | [Foods 2022, doi:10.3390/foods11142004](https://doi.org/10.3390/foods11142004) | CC BY 4.0 | `ee6da3d22db941fe` |

The last two are **published tables inside open-access articles**. This repository
ships the deterministic extraction (`pandas.read_html`) of those tables, checked
value by value against the Europe PMC full text and by hand. The avocado xlsx is
read directly with no numerical transformation and is **not** re-hosted here.

## Derived artefacts (shipped in this repository)

Every number in the manuscript traces to one of these files.

| File | Content | SHA-256 (16) |
|---|---|---|
| `43trichotomy_formal.json` | Semi-synthetic trichotomy: terminal Ψ, collapse slopes, Krug-line match | `3ddee85bf7b57f03` |
| `54c5_realkinetics.json` | Apple system: per-cultivar Eₐ and R², design sweep, over-determined subsets | `d2f7f2b2cd777d7f` |
| `54c5_realkinetics-designsweep.csv` | The 26 apple temperature subsets, one row each (source of Table S2) | `2b1769aa0b095d56` |
| `91c5_realkinetics_killer.json` | Ψ versus inexpensive heuristics; Ψ_het across systems | `64fc189dacbb29c7` |
| `99_realkinetics2.json` | Food system: Eₐ range, homogeneous / heterogeneous Ψ replication | `7b541f1778481b00` |
| `102_avocado_c5_pilot.json` | Avocado system: genuine statistical estimation error, sweep along n | `125ba608e1256276` |
| `119pazman_pronzato_bridge.json` | Coordinate check against the Pázman–Pronzato (2006) worked example | `6d5f827620897f26` |

Semi-synthetic data are not shipped as files — they are regenerated
deterministically from fixed seeds by the scripts.

Each dataset is used under its own license; please cite the original sources.
`code/fetch_datasets.py --list` prints the same pointers.
