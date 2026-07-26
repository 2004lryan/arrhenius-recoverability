# Data card — recoverability trichotomy for Arrhenius inversion

Covers every dataset used by the manuscript. Fields follow the project data-card
standard (source & license, SHA-256, granularity & units, sample statistics,
split logic, known biases, label reliability, preprocessing pointer, ethics).

Generated from `code/120datacard_sha256.py`; the full 64-hex digests below are
the authoritative record. **No new data were collected for this work.**

> **What SHA-256 does and does not certify.** The digests certify that the files
> have not changed *since registration*. They say nothing about the period before
> registration, and they are not evidence that the extraction step is correct —
> extraction correctness rests on `pandas.read_html` plus the value-by-value
> re-check described under each system.

---

## System 1 — Apple firmness under a temperature gradient

- **Source & license.** Table A1 of an open-access article, Foods 2023,
  [doi:10.3390/foods12112113](https://doi.org/10.3390/foods12112113)
  (PMC10253207). CC BY 4.0. Collected by another group, not by us.
- **SHA-256.** `d22131f037a120e71a321d358e141ea0b8b2b8b97c4299706f4b0c8fad66665f`
  (`PMC10253207_T16.csv`, 1,770 bytes)
- **Granularity & units.** One row per (storage temperature, elapsed time).
  Temperature in °C; time in hours; firmness in kg·cm⁻², reported as
  mean ± standard deviation.
- **Sample statistics.** 25 data rows (the file carries a second header row):
  5 storage temperatures (2, 4, 6, 8, 10 °C) × 5 time points (0, 12, 24, 36,
  48 h), with one firmness column per cultivar — 4 cultivars (Granny Smith,
  Pink Lady, Red Delicious, Royal Gala). Only published per-temperature means
  are available — **there are no replicate observations.**
- **Split logic.** No train/test split. The design-window sweep enumerates all
  26 temperature subsets of size ≥ 2; these subsets are **nested and therefore
  not independent**, which is why the paper reports rank correlations without
  p-values.
- **Known biases.** (i) Pink Lady is non-monotone — firmness rises first at 2 °C
  through physiological after-ripening — and yields an unstable inversion
  (Eₐ = 858 kJ/mol, R² = 0.55); the paper keeps it precisely as the pathological
  case the diagnostic should flag. (ii) The temperature range is narrow
  (2–10 °C), which is the regime the paper is about, not a defect.
- **Label reliability.** Firmness is an instrumental measurement reported with
  its standard deviation by the original authors; we did not re-measure it.
- **Preprocessing pointer.** `code/55fetch_pmc_tables.py` (extraction via
  `pandas.read_html`, then checked by hand row by row against the published
  table); `code/54c5_realkinetics.py` and `code/91c5_realkinetics_killer.py`
  (inversion, design sweep, heuristic comparison).
- **Ethics.** Plant material only. No human or animal subjects; no personal data.

## System 2 — Long-shelf-life food kinetics (ammonia / TBARS)

- **Source & license.** Table 3 of an open-access article, Foods 2022;11(14):2004,
  [doi:10.3390/foods11142004](https://doi.org/10.3390/foods11142004)
  (PMC9319022). CC BY 4.0. Collected by another group.
- **SHA-256.** `ee6da3d22db941fe2936bba9817e3aa7873e408024ae9ab95d6bf48d536ccfbe`
  (`PMC9319022_T3.csv`, 5,239 bytes)
- **Granularity & units.** One row per (food, storage time, storage temperature).
  Time in months; temperature in °C; ammonia content in mg/kg; TBARS in A538/mg.
  Values are mean ± standard deviation with the original significance letters kept.
  Note: the published table wraps long food names across lines, so the raw
  `Sample` column is fragmented (e.g. "Instant" / "goulash" / "soup" are three
  pieces of one name); the food identity is reassembled downstream, not in the
  extracted CSV.
- **Sample statistics.** 145 data rows (the file carries a second header row);
  5 foods × 4 storage temperatures (−18, 5, 25, 40 °C; the initial time point is
  recorded as "-"), giving 10 usable kinetic sequences.
- **Split logic.** No split; each (food, response) pair is inverted separately.
- **Known biases.** Two chemically distinct responses (ammonia accumulation and
  lipid oxidation) are pooled as independent kinetic sequences; they share the
  same storage batches, so the sequences are not fully independent across
  responses within a food.
- **Label reliability.** Published means with standard deviations; the extraction
  was **independently re-scraped from the Europe PMC `fullTextXML` and compared
  value by value — 132 rows, 0 discrepancies.**
- **Preprocessing pointer.** `code/55fetch_pmc_tables.py` (extraction),
  `code/99_realkinetics2.py` (inversion).
- **Ethics.** Processed food products. No human or animal subjects.

## System 3 — Hass avocado postharvest ripening

- **Source & license.** Mendeley Data,
  [doi:10.17632/3xd9n945v8.1](https://doi.org/10.17632/3xd9n945v8.1). CC BY 4.0.
  Collected by another group. **Referenced read-only and not re-hosted here.**
- **SHA-256.** `f8abeaba6eedf67869907b3a71ee132cd89b42bb3dcce6b25ddb22e6269b51a7`
  (`Avocado Ripening Dataset.xlsx`, 633,351 bytes)
- **Granularity & units.** One row per fruit per photography day. Ripening index
  on a 5-level ordinal scale; storage temperature in °C; time in days.
- **Sample statistics.** 478 fruit × 3 storage conditions; 14,722 records of
  daily two-sided photography. At the two nominally known temperatures (10 °C
  and 20 °C) there are 155 and 127 **fully independent individuals**
  (282 in total), with between-individual standard deviations of 2.7 and 0.6
  days in days-to-ripeness.
- **Split logic.** Resampling is at the **individual** level (cluster bootstrap
  over the 282 independent fruit), never at the record level — records from one
  fruit are correlated.
- **Known biases.** Only two of the three storage conditions have a nominally
  known temperature, so this system covers the sample-size axis of
  Ψ = n·(ΔT)² but **not** the design-window axis. The paper states this
  limitation explicitly (§Scope, L4/L5).
- **Label reliability.** The ripening index is a visual ordinal grading by the
  original authors; the derived rate is the reciprocal of the number of days to
  first reach level 5. Its dispersion is genuine biological variation, not
  injected noise.
- **Preprocessing pointer.** `code/102_avocado_c5_pilot.py` — the original xlsx
  is read directly with **no numerical transformation**.
- **Ethics.** Plant material only; the photographs contain no people.

## Semi-synthetic data

- **Source.** Generated, not measured. Produced deterministically from fixed
  seeds by `code/43trichotomy_formal.py`; not shipped as files.
- **Role.** The only evidence for the phase transition itself. The real systems
  validate the design rule, **not** the transition — stated as such in the paper.
- **Ethics.** Not applicable.

---

## Derived artefacts

The SHA-256 digest of every derived result file backing a number in the
manuscript is listed in `../data/DATA.md` and reproduced by
`code/120datacard_sha256.py`.
