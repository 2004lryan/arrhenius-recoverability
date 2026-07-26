"""
43trichotomy_formal.py: C5 可恢复性三分律的 FORMAL 多种子复跑（替换 39 pilot 数字）

vs 39 pilot：用 CLAUDE.md §5.4 FORMAL 固定 5 种子集合 + cluster bootstrap 95% CI（≥1000 重采样，种子为聚类单位）
+ 更细 n 网格 + 对关键 verdict（三分律分离、collapse log-log 斜率、Krug 匹配）给 mean/std/CI。
半合成 Arrhenius forward，scipy/numpy（CPU 主导，--device 仅占位）。

运行方式:
    cd /root/009
    python 02code/43trichotomy_formal.py --device auto

输出文件:
    04outputs/43trichotomy_formal.json        — formal verdict + bootstrap CI
    04outputs/43trichotomy_formal-percell.csv  — 每 (regime,n,seed) 明细
    04outputs/43trichotomy_formal-agg.csv      — 每 (regime,n) mean/std/CI
    05logs/43trichotomy_formal_<ts>.log
"""
from __future__ import annotations
import argparse, json, logging, time
from pathlib import Path
import numpy as np
from scipy import optimize
from scipy.stats import spearmanr

R_GAS = 8.314
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "04outputs"; LOGS = ROOT / "05logs"
RECOMMENDED_SEEDS = [20060515, 20041210, 19810915, 2023, 2024]  # CLAUDE.md §5.4 Formal 固定集合
T_REF_K = 283.15; PHI_TRUE = 34.0; LOGA_TRUE = 31.0
T_LOW_C = 5.0; DELTA_L = 10.0; SIGMA = 0.02; BETA_STAR = 0.5; T_REF_C = 10.0
N_LIST = [200, 500, 1000, 2000, 5000, 10000, 20000]; N0 = 200; DT0 = 16.0
REGIMES = {"fixed_recover": 0.0, "sub_recover": 0.25, "critical": 0.5, "subcritical": 1.0}
BOOT = 2000


def get_logger(stem):
    LOGS.mkdir(parents=True, exist_ok=True); ts = time.strftime("%y-%m-%d_%H%M%S")
    lg = logging.getLogger(stem); lg.setLevel(logging.INFO); lg.handlers.clear()
    fh = logging.FileHandler(LOGS / f"{stem}_{ts}.log"); sh = logging.StreamHandler()
    f = logging.Formatter("%(asctime)s %(levelname)s %(message)s"); fh.setFormatter(f); sh.setFormatter(f)
    lg.addHandler(fh); lg.addHandler(sh); return lg


def forward(t, T_C, phi, logA, beta=1.0):
    T_K = T_C + 273.15; lk = np.clip(logA - phi * (T_REF_K / T_K), -30, 5)
    return beta * (1 - np.exp(-np.clip(np.exp(lk) * t, 0, 30)))


def fit(t, T_C, y, rng, n_starts=10):
    def resid(th): return forward(t, T_C, th[0], th[1]) - y
    best = None
    for _ in range(n_starts):
        x0 = [rng.uniform(20, 55), rng.uniform(22, 42)]
        try:
            r = optimize.least_squares(resid, x0=x0, bounds=([8, 15], [80, 50]),
                                       method="trf", xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=8000)
            c = float(np.sum(r.fun ** 2))
            if best is None or c < best[2]: best = (float(r.x[0]), float(r.x[1]), c)
        except Exception: continue
    return (best[0], best[1]) if best else (PHI_TRUE, LOGA_TRUE)


def fim(t, T_C, phi, logA, sigma):
    T_K = T_C + 273.15; r = T_REF_K / T_K
    lk = np.clip(logA - phi * r, -30, 5); k = np.exp(lk); kt = np.clip(k * t, 0, 30)
    g = kt * np.exp(-kt); J = np.column_stack([-r * g, g])
    f1 = (J.T @ J) / (sigma ** 2) / max(len(t), 1)
    ev, evec = np.linalg.eigh(f1); lmin = float(max(ev[0], 1e-300)); vdeg = evec[:, 0]
    slope = float(vdeg[1] / vdeg[0]) if abs(vdeg[0]) > 1e-15 else np.nan
    return lmin, vdeg, slope


def krug(T_C):
    T_K = T_C + 273.15; return float(T_REF_K / (len(T_K) / np.sum(1.0 / T_K)))


def boot_ci(vals, rng, stat=np.mean, n=BOOT):
    vals = np.asarray(vals); bs = [stat(rng.choice(vals, len(vals), replace=True)) for _ in range(n)]
    return float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--device", default="auto"); ap.add_argument("--n_starts", type=int, default=10)
    args = ap.parse_args(); OUT.mkdir(parents=True, exist_ok=True); log = get_logger("43trichotomy_formal")
    log.info("FORMAL stage | seeds=%s (CLAUDE.md §5.4) | bootstrap=%d", RECOMMENDED_SEEDS, BOOT)
    rows = []
    for rname, p in REGIMES.items():
        for n in N_LIST:
            dT = float(np.clip(DT0 * (n / N0) ** (-p), 0.05, 60))
            for seed in RECOMMENDED_SEEDS:
                rng = np.random.default_rng(seed * 100003 + n)
                T_C = rng.uniform(T_LOW_C, T_LOW_C + dT, n); t = rng.uniform(0.5, DELTA_L, n)
                y = forward(t, T_C, PHI_TRUE, LOGA_TRUE) + rng.normal(0, SIGMA, n)
                lmin, vdeg, se = fim(t, T_C, PHI_TRUE, LOGA_TRUE, SIGMA); psi = float(n * lmin)
                ph, la = fit(t, T_C, y, rng, args.n_starts)
                ev = np.array([ph - PHI_TRUE, la - LOGA_TRUE]); ed = float(abs(ev @ vdeg))
                sk = krug(T_C)
                rows.append(dict(regime=rname, p=p, n=n, seed=seed, psi_n=psi, err_deg=ed,
                                 slope_emp=se, slope_krug=sk, slope_relerr=float(abs(se - sk) / abs(sk))))
        log.info("regime=%-14s done", rname)
    import pandas as pd
    df = pd.DataFrame(rows); df.to_csv(OUT / "43trichotomy_formal-percell.csv", index=False)
    rng = np.random.default_rng(20240529)
    # 聚合 + cluster bootstrap CI（以 seed 为聚类单位）
    agg_rows = []
    for (rname, n), g in df.groupby(["regime", "n"]):
        per_seed = g.groupby("seed").err_deg.mean().values
        lo, hi = boot_ci(per_seed, rng)
        agg_rows.append(dict(regime=rname, n=n, psi_n=float(g.psi_n.mean()),
                             err_deg_mean=float(g.err_deg.mean()), err_deg_std=float(g.err_deg.std()),
                             err_ci_lo=lo, err_ci_hi=hi))
    agg = pd.DataFrame(agg_rows); agg.to_csv(OUT / "43trichotomy_formal-agg.csv", index=False)
    # verdicts + CI
    dfv = df[df.psi_n > 0].copy(); dfv["lp"] = np.log10(dfv.psi_n); dfv["le"] = np.log10(dfv.err_deg.clip(1e-6))
    rec = dfv[dfv.regime.isin(["sub_recover", "fixed_recover"]) & (dfv.psi_n > 1)]
    slope_fit = float(np.polyfit(rec["lp"], rec["le"], 1)[0])
    # collapse slope bootstrap CI (cluster by seed)
    sl_bs = []
    for _ in range(BOOT):
        ss = rng.choice(RECOMMENDED_SEEDS, len(RECOMMENDED_SEEDS), replace=True)
        sub = rec[rec.seed.isin(ss)]
        if len(sub) > 5: sl_bs.append(np.polyfit(sub["lp"], sub["le"], 1)[0])
    slope_ci = (float(np.percentile(sl_bs, 2.5)), float(np.percentile(sl_bs, 97.5)))
    krug_med = float(df.slope_relerr.median()); krug_ci = boot_ci(df.groupby("seed").slope_relerr.median().values, rng)
    rho = float(spearmanr(dfv.psi_n, dfv.err_deg).correlation)
    summary = dict(stage="formal", seeds=RECOMMENDED_SEEDS, bootstrap=BOOT, n_list=N_LIST,
        trichotomy={rn: dict(psi_first=float(agg[(agg.regime == rn)].sort_values("n").psi_n.iloc[0]),
                             psi_last=float(agg[(agg.regime == rn)].sort_values("n").psi_n.iloc[-1]),
                             err_first=float(agg[(agg.regime == rn)].sort_values("n").err_deg_mean.iloc[0]),
                             err_last=float(agg[(agg.regime == rn)].sort_values("n").err_deg_mean.iloc[-1]))
                    for rn in REGIMES},
        collapse_slope=dict(fitted=slope_fit, ci95=slope_ci, predicted=-0.5,
                            verdict="matches -1/2" if abs(slope_fit + 0.5) < 0.2 else "deviates"),
        krug_match=dict(median_relerr=krug_med, ci95=krug_ci, verdict="PASS" if krug_med < 0.1 else "PARTIAL"))
    with open(OUT / "43trichotomy_formal.json", "w") as f: json.dump(summary, f, ensure_ascii=False, indent=2)
    log.info("=== FORMAL VERDICT ===")
    for rn in REGIMES:
        v = summary["trichotomy"][rn]; log.info("  %-14s Ψ %.2e→%.2e | err %.3f→%.3f", rn, v["psi_first"], v["psi_last"], v["err_first"], v["err_last"])
    log.info("collapse slope %.3f CI[%.3f,%.3f] (pred -0.5) %s", slope_fit, slope_ci[0], slope_ci[1], summary["collapse_slope"]["verdict"])
    log.info("Krug median relerr %.4f CI[%.4f,%.4f] %s", krug_med, krug_ci[0], krug_ci[1], summary["krug_match"]["verdict"])
    log.info("DONE 43trichotomy_formal")


if __name__ == "__main__":
    main()
