# Data cards only — no raw data in this repository

This directory holds the **data card** for the manuscript:

- [`datacard-c5.md`](datacard-c5.md) — all three real systems plus the
  semi-synthetic data, with source & license, full SHA-256 digests, granularity
  & units, sample statistics, split logic, known biases, label reliability,
  preprocessing pointers, and ethics status.

The raw data are **not** here. No new data were collected for this work: every
real system is previously published third-party data under CC BY 4.0, obtainable
from the DOIs listed in [`../data/DATA.md`](../data/DATA.md). The large avocado
xlsx is referenced through its Mendeley DOI and is not re-hosted. There is no
restricted or self-collected data in this project.

The derived result files backing each number in the paper *are* shipped, with
their SHA-256 digests; regenerate the digests with
`python code/120datacard_sha256.py`.
