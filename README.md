# Recoverability trichotomy for Arrhenius inversion

Code and reproducibility materials for the manuscript

> **Recoverability trichotomy for Arrhenius inversion: the kinetic compensation line as a Fisher null direction, and its consequences for experimental design**
> (under review, *Chemometrics and Intelligent Laboratory Systems*)

## What the paper claims

When activation energies are inverted from a narrow temperature window, the Arrhenius
parameters become almost indistinguishable along one direction and extra sampling does not
repair it. Recoverability is governed by one computable scalar

$$\Psi_n = n\,\lambda_{\min}(\hat I_1),$$

the sample size times the smallest eigenvalue of the normalised Fisher information, and it is
strictly trichotomous: as the design window shrinks like a power $\rho$ of the sample size, the
error vanishes at rate $\Psi^{-1/2}$ for $\rho<1/2$, is constant at $\rho=1/2$, and stays bounded
away from zero for $\rho>1/2$. The degenerate direction is the classical kinetic compensation line.

**What we do not claim**, and the code shows why:

- $\Psi$ is, as a criterion, classical **E-optimality**; no novelty is claimed for that.
- The branch $\Psi\to\infty$ is the established least-squares consistency criterion.
- The compensation-line correspondence is an identity.
- The design prescription (widen the window) is not new as a prescription.
- Three claims of ours are **falsified in the paper itself**, by these scripts:
  $\Psi$ has no discriminating power among $\Psi$-equivalent designs; a D-optimality
  comparison we withdrew; and $\Psi$ as a *calibration* diagnostic (`124coverage_study.py`).

## Repository layout

```
code/       analysis scripts, numbered as in the manuscript
outputs/    every result file behind a number in the paper
data/       extracted source tables (CC BY 4.0) + SHA-256 manifest
figures/    the four main-text figures
docs/       the data card
tools/      pre-push privacy guard
```

## Reproduction

```bash
pip install -r requirements.txt

python code/123montecarlo_criterion_study.py   # 1000 independent designs, ~105 s
python code/124coverage_study.py               # CI coverage study
python code/122criterion_bakeoff.py            # E / D / A / c comparison
python code/120datacard_sha256.py              # regenerate the SHA-256 manifest
python code/121audit_en_assurance.py           # manuscript/number consistency audit
```

All runs are single-threaded CPU, seeded, and complete in minutes; no GPU is used.

## Data

No new data were collected. All three real systems are previously published third-party
data under CC BY 4.0 — see [`data/DATA.md`](data/DATA.md) for sources and DOIs, and
[`docs/datacard-c5.md`](docs/datacard-c5.md) for the full data card. The two literature
tables are shipped here as extractions; the avocado dataset is large and is referenced
through its Mendeley DOI rather than re-hosted.

**Not in this repository, by design:** the group's own unpublished spectral datasets, which
belong to other projects and are not used by this paper. `tools/check_no_leak.sh` fails the
push if any of them, any raw-data blob, or any local absolute path appears.

## Citation

See [`CITATION.cff`](CITATION.cff). A BibTeX entry will be added on acceptance.

## License

Code: MIT (see [`LICENSE`](LICENSE)). Each dataset remains under its own licence; please
cite the original sources.

## Contact

Hua Huang (corresponding author) — huanghua@xjau.edu.cn — Xinjiang Agricultural University.
