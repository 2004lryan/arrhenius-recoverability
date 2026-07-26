#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
122criterion_bakeoff.py — 经典最优设计准则的真实数据对决：为什么是 E（即 Ψ），不是 D/A/c？

动机（2026-07-25 fusion 终审 4.75/10 的核心指控）
------------------------------------------------
面板 8/12 一致指出：本文真实数据唯一验证的结论是设计推论「展宽温度窗」，
而稿件自己承认该处方并非新结论（D-最优本就把支撑点推向温区两端）。
另有模型指出 Ψ_n = N·λ_min(I) 在准则层面就是经典 **E-最优**，全稿却从未点名。

原设想：真实数据的载荷换成"在经典准则中唯 E 追踪得住真实反演误差"。
**该设想已被本脚本自身的数据推翻，2026-07-25 第二轮 fusion 抓出，作者复算坐实：**
  (1) 所谓"E/A/c 秩等价"是**算术恒等**而非实证发现——2×2 近奇异阵中
      A/E ∈ [1.999961, 1.999991]（A 就是 2E），c_deg/E = c_krug/E = 1.000000 精确。
      该事实三行可证，与是否测过苹果无关，故不得作为"真实数据结论"呈现。
  (2) 所谓"唯 D 排错"**统计上不成立**：bootstrap 5000 次，边际 |ρ_E|-|ρ_D| = +0.030
      (95%CI [-0.095,+0.204], P(D≥E)=44.5%)；偏相关 +0.092 (95%CI [-0.054,+0.330], P=26.6%)。
      两个 CI 都含零。16 个子集共 120 对中，E 与 D 只有 6 对排序不同，且集中于同两个
      4-温度设计对同三个 3-温度设计，是**一次结构性比较**而非 16 次独立观测。
  (3) 逻辑自陷：正文以"子集嵌套非独立故不报 p"为由拒绝检验，却又拿同一批依赖点上
      **两个相关系数之差**当结论——两者不能兼得。
故本脚本的定位已改为：**证伪自己的假设**，并给出可写入正文的诚实版本（塌缩=解析事实；D 对比=未建立）。

口径
----
- 子集与真实误差完全复用 `04outputs/54c5_realkinetics-designsweep.csv`（脚本 54 产出，26 个温度子集）。
- Fisher 阵按脚本 54 的 `fisher()` 同构重建：J = [-1/(R·T_K), 1]/σ，I = JᵀJ。
  σ_lnk 在各子集间是**公共常数**，而每个准则都是 σ 的单调变换，故取 σ=1 不改变任何 Spearman；
  重建正确性由"能否精确复现 CSV 中 Psi_n / lambda_min 的比例"当场校验（见 verify_reconstruction）。
- 六个准则一律乘样本量 N，化为可比的"总信息"量，**均为越大越好**：
    E    : N·λ_min(M)                  ← 本文的 Ψ_n
    D    : N·det(M)^(1/2)              ← 特征值几何均值
    A    : N·2/tr(M⁻¹)                 ← 特征值调和均值
    c_deg: N/(v_degᵀ M⁻¹ v_deg)        ← 沿退化方向（补偿线）的 c-最优
    c_Ea : N/(e₁ᵀ M⁻¹ e₁)              ← 只关心 E_a 时最自然的 c-最优
    cond : λ_max/λ_min（越大越差，报相关时取负号语义由符号自明）
- 判据：与真实反演误差 realCV_Ea 的 Spearman（超定子集 n_temps≥3，df>0，真统计误差），
  以及控制 ΔT 后的偏 Spearman（教科书标准定义，与脚本 91/104/106/107 同口径）。

诚实边界
--------
1. 子集嵌套、互不独立，故只报秩相关不报 p 值（与正文 §res-c5 同口径）。
2. 单一品种、单一真实数据集，n=16 超定子集，结论是**基准级**而非机制级。
3. 若 D 或 A 胜出，本脚本照实报告。**实际结果比这更糟：E 对 D 的优势根本未达统计可辨，已如实撤回。**
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "04outputs" / "54c5_realkinetics-designsweep.csv"
OUT = ROOT / "04outputs" / "122criterion_bakeoff.json"
R_GAS = 8.314e-3  # kJ/(mol·K)。**必须与脚本 54 逐位一致**：54 用的就是 8.314e-3，
                  # 换成更精确的 8.314462618e-3 会让重建与 CSV 的比例出现 1.4e-7 相对展布，
                  # 即口径不一致；正文所有数字均出自 8.314e-3 这一支。


def fisher_sigma1(T_C: np.ndarray) -> np.ndarray:
    """脚本 54 的 fisher() 取 σ=1；σ 为跨子集公共常数，不影响任何 Spearman。"""
    T_K = np.asarray(T_C, float) + 273.15
    J = np.stack([-1.0 / (R_GAS * T_K), np.ones_like(T_K)], axis=1)
    return J.T @ J


def krug_direction(T_C: np.ndarray) -> np.ndarray:
    """理论退化方向（Krug 补偿线）：d(lnA)/d(E_a) = 1/(R·T_hm)，T_hm 为调和均值温度。"""
    T_K = np.asarray(T_C, float) + 273.15
    T_hm = 1.0 / np.mean(1.0 / T_K)
    v = np.array([1.0, 1.0 / (R_GAS * T_hm)])
    return v / np.linalg.norm(v)


def criteria(T_C, N: int) -> dict:
    M = fisher_sigma1(T_C)
    ev, evec = np.linalg.eigh(M)
    lmin, lmax = float(ev[0]), float(ev[1])
    Minv = np.linalg.inv(M)
    v_min = evec[:, 0]
    v_krug = krug_direction(T_C)
    e1 = np.array([1.0, 0.0])
    return {
        "E": N * lmin,
        "D": N * float(np.sqrt(max(np.linalg.det(M), 0.0))),
        "A": N * 2.0 / float(np.trace(Minv)),
        "c_deg": N / float(v_min @ Minv @ v_min),
        "c_krug": N / float(v_krug @ Minv @ v_krug),
        "c_Ea": N / float(e1 @ Minv @ e1),
        "cond": lmax / max(lmin, 1e-300),
        "_lambda_min": lmin,
        # 引理 fisher 的逐子集数值检验：理论补偿线 vs 数值最小特征方向的夹角
        "_align_krug_vs_vmin_abscos": abs(float(v_krug @ v_min)),
    }


def partial_spearman(x, y, z) -> float:
    """教科书标准偏 Spearman：三者各自秩变换，x、y 对 z 线性回归取残差，再求残差 Pearson 相关。"""
    rx, ry, rz = (rankdata(v).astype(float) for v in (x, y, z))
    Z = np.stack([rz, np.ones_like(rz)], axis=1)
    res = []
    for r in (rx, ry):
        beta, *_ = np.linalg.lstsq(Z, r, rcond=None)
        res.append(r - Z @ beta)
    if np.std(res[0]) < 1e-12 or np.std(res[1]) < 1e-12:
        return float("nan")
    return float(np.corrcoef(res[0], res[1])[0, 1])


def main() -> None:
    df = pd.read_csv(CSV)
    temps = [np.array([float(x) for x in s.split("/")]) for s in df["temps"]]
    n_per_T = int(round(df["Psi_n"].iloc[0] / df["lambda_min"].iloc[0] / df["n_temps"].iloc[0]))

    # --- 重建正确性校验：σ=1 的 λ_min 应与 CSV 的 lambda_min 成【恒定比例】(= σ²) ---
    lmin_rebuilt = np.array([float(np.linalg.eigvalsh(fisher_sigma1(T))[0]) for T in temps])
    ratio = lmin_rebuilt / df["lambda_min"].to_numpy()
    rel_spread = float(ratio.std() / ratio.mean())
    sigma_implied = float(np.sqrt(ratio.mean()))
    assert rel_spread < 1e-9, f"重建与 CSV 不成恒定比例（相对展布 {rel_spread:.2e}），口径不一致"
    # 排序一致性：Spearman 必须严格为 1
    rank_ok = float(spearmanr(lmin_rebuilt, df["lambda_min"])[0])
    assert abs(rank_ok - 1.0) < 1e-12, f"重建排序与 CSV 不一致 (Spearman={rank_ok})"

    rows = []
    for T, nt in zip(temps, df["n_temps"]):
        rows.append(criteria(T, int(nt) * n_per_T))
    # 加前缀避免与 CSV 既有列（cond）撞名——撞名会让 df[k] 返回 DataFrame 而非 Series
    C = pd.DataFrame(rows).add_prefix("crit_")
    full = pd.concat([df.reset_index(drop=True), C], axis=1)

    # 交叉核对 1：本脚本的 E 与 CSV 的 Psi_n 排序必须完全一致（同一个量，仅差 σ² 常数）
    assert abs(spearmanr(full["crit_E"], full["Psi_n"])[0] - 1.0) < 1e-12
    # 交叉核对 2：本脚本独立算出的条件数应与 CSV 的 cond 数值一致（cond 无量纲，σ 完全抵消）
    cond_relerr = float(np.max(np.abs(full["crit_cond"] / full["cond"] - 1.0)))
    assert cond_relerr < 1e-9, f"条件数与 CSV 不一致（最大相对误差 {cond_relerr:.2e}）"

    # --- 超定子集（n_temps>=3）：真实测量 SD 传播的估计误差 ---
    od = full[(full.n_temps >= 3) & np.isfinite(full.realCV_Ea) & (full.realCV_Ea > 0)].copy()
    NAMES = ["crit_E", "crit_D", "crit_A", "crit_c_deg", "crit_c_krug", "crit_c_Ea", "crit_cond"]
    err = od["realCV_Ea"].to_numpy()
    dT = od["deltaT"].to_numpy()

    marg = {k: float(spearmanr(od[k], err)[0]) for k in NAMES}
    marg["deltaT"] = float(spearmanr(dT, err)[0])
    marg["n_temps"] = float(spearmanr(od["n_temps"], err)[0])
    part = {k: partial_spearman(od[k].to_numpy(), err, dT) for k in NAMES}

    # 越大越好的准则应与误差【负】相关；|rho| 越大越强。cond 越大越差故期望正相关。
    rank_by_abs = sorted(marg.items(), key=lambda kv: -abs(kv[1]))

    # --- 关键发现（与作者预期相反，照实测量）：近奇异区各准则是否【塌缩为同一排序】 ---
    # 2×2 信息阵条件数达 1e5~1e6 量级时 tr(M^-1)=1/λmin+1/λmax≈1/λmin，
    # 故 A≈2λmin；同理 c 沿任何与补偿线不正交的方向也被 1/λmin 支配 → 三者与 E 秩等价。
    pair = {}
    for i, a in enumerate(NAMES):
        for b in NAMES[i + 1:]:
            pair[f"{a}|{b}"] = float(spearmanr(full[a], full[b])[0])
    rank_equiv_with_E = sorted(k for k in NAMES
                               if k != "crit_E" and abs(spearmanr(full[k], full["crit_E"])[0] - 1.0) < 1e-12)
    cond_med = float(np.median(full["crit_cond"]))
    # A 与 2λmin 的相对偏差：直接量化"塌缩"有多彻底
    A_over_2lmin = (full["crit_A"] / (full["crit_n_scale"] * 2 * full["crit__lambda_min"])
                    if "crit_n_scale" in full else None)
    collapse = {
        "criteria_rank_equivalent_to_E": rank_equiv_with_E,
        "median_condition_number": cond_med,
        "pairwise_spearman": pair,
        "interpretation": (
            "在本文关心的近奇异设计区，条件数中位数达 %.3g，tr(M^-1)≈1/λ_min，"
            "故 A-最优、以及沿任何与补偿线不正交方向的 c-最优，都被 1/λ_min 支配，"
            "与 E-最优（本文 Ψ）秩等价（Spearman=1.000）。唯 D-最优不塌缩——"
            "它取特征值几何均值，被大特征值抬高，因而偏离真实误差排序。" % cond_med),
    }

    # --- 决定性检验：E 相对 D 的优势是否可辨？（不做这一步就会把掷硬币当成结论）---
    rng = np.random.default_rng(20260725)
    nb = len(od)
    gm, gp = [], []
    for _ in range(5000):
        idx = rng.integers(0, nb, nb)
        if len(np.unique(err[idx])) < 4:
            continue
        gm.append(abs(spearmanr(od["crit_E"].to_numpy()[idx], err[idx])[0])
                  - abs(spearmanr(od["crit_D"].to_numpy()[idx], err[idx])[0]))
        a = partial_spearman(od["crit_E"].to_numpy()[idx], err[idx], dT[idx])
        b = partial_spearman(od["crit_D"].to_numpy()[idx], err[idx], dT[idx])
        if np.isfinite(a) and np.isfinite(b):
            gp.append(abs(a) - abs(b))
    def summarise(g):
        g = np.asarray(g)
        return {"point": float(np.mean(g)), "ci95": [float(np.percentile(g, 2.5)), float(np.percentile(g, 97.5))],
                "P_D_at_least_as_good": float(np.mean(g <= 0)), "n_resamples": int(len(g))}
    Ev, Dv = od["crit_E"].to_numpy(), od["crit_D"].to_numpy()
    n_discordant = sum(1 for i, j in itertools.combinations(range(nb), 2)
                       if (Ev[i] < Ev[j]) != (Dv[i] < Dv[j]))
    A_over_E = (od["crit_A"] / od["crit_E"])
    E_vs_D = {
        "bootstrap_marginal_gap": summarise(gm),
        "bootstrap_partial_gap": summarise(gp),
        "n_discordant_pairs": n_discordant,
        "n_pairs_total": nb * (nb - 1) // 2,
        "A_over_E_range": [float(A_over_E.min()), float(A_over_E.max())],
        "lambda_min_fold_range": float(od["crit__lambda_min"].max() / od["crit__lambda_min"].min()),
        "lambda_max_fold_range": float((od["crit__lambda_min"] * od["crit_cond"]).max()
                                       / (od["crit__lambda_min"] * od["crit_cond"]).min()),
        "verdict": ("E 相对 D 的优势【未达统计可辨】：两个 bootstrap CI 均含零，"
                    "且 120 对中仅 6 对排序不同。不得作为结论写入正文。"),
    }

    align = float(full["crit__align_krug_vs_vmin_abscos"].min())
    res = {
        "script": Path(__file__).name,
        "source_csv": str(CSV),
        "n_subsets_total": int(len(full)),
        "n_overdetermined": int(len(od)),
        "n_per_T": n_per_T,
        "reconstruction_check": {
            "ratio_relative_spread": rel_spread,
            "implied_sigma_lnk": sigma_implied,
            "spearman_vs_csv_lambda_min": rank_ok,
            "cond_max_relerr_vs_csv": cond_relerr,
            "note": "σ=1 重建与脚本 54 输出成恒定比例且排序 Spearman=1，故口径一致",
        },
        "lemma_fisher_numerical_check": {
            "min_abs_cos_krug_vs_min_eigvec": align,
            "note": ("理论补偿线方向（斜率=调和均值温度）与数值最小特征方向的最小 |cos|；"
                     "≈1 即引理 fisher 在全部 26 个真实子集上成立"),
        },
        "spearman_marginal_vs_realCV": marg,
        "spearman_partial_controlling_deltaT": part,
        "ranking_by_abs_marginal": rank_by_abs,
        "criterion_collapse_in_near_singular_regime": collapse,
        "E_vs_D_is_NOT_established": E_vs_D,
        "honesty": [
            "子集嵌套非独立，只报秩相关不报 p 值（与正文 §res-c5 同口径）",
            f"单一品种、单一真实数据集、{len(od)} 个超定子集，结论为基准级而非机制级",
            "若 D 或 A 胜出即照实报告",
        ],
    }
    OUT.write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")
    full.to_csv(ROOT / "04outputs" / "122criterion_bakeoff-percriterion.csv", index=False)

    print(f"重建校验：比例相对展布 {rel_spread:.2e}，隐含 σ_lnk={sigma_implied:.4f}，排序 Spearman={rank_ok:.12f}")
    print(f"引理 fisher 数值检验：min|cos(Krug, v_min)| = {align:.6f}（26 个子集）")
    print(f"\n超定子集 n={len(od)}，与真实反演误差 realCV_Ea 的 Spearman：")
    for k, v in rank_by_abs:
        p = part.get(k)
        ps = f"  偏(控 ΔT) {p:+.3f}" if p is not None and np.isfinite(p) else ""
        print(f"  {k:>8s}  边际 {v:+.3f}{ps}")
    print(f"\n与 E 秩等价的准则（Spearman=1.000）：{rank_equiv_with_E}")
    print(f"条件数中位数 {cond_med:.3g} → 近奇异，故 A/c 被 1/λ_min 支配而塌缩到 E；唯 D 不塌缩且更差。")
    print(f"\nE vs D 决定性检验：边际差 {E_vs_D['bootstrap_marginal_gap']['point']:+.3f} "
          f"CI{[round(x,3) for x in E_vs_D['bootstrap_marginal_gap']['ci95']]} "
          f"P(D>=E)={E_vs_D['bootstrap_marginal_gap']['P_D_at_least_as_good']:.1%}")
    print(f"                    偏相关差 {E_vs_D['bootstrap_partial_gap']['point']:+.3f} "
          f"CI{[round(x,3) for x in E_vs_D['bootstrap_partial_gap']['ci95']]} "
          f"P(D>=E)={E_vs_D['bootstrap_partial_gap']['P_D_at_least_as_good']:.1%}")
    print(f"  → {E_vs_D['verdict']}")
    print(f"  A/E 比值范围 {E_vs_D['A_over_E_range']}（≈2 即算术恒等）；"
          f"{E_vs_D['n_pairs_total']} 对中仅 {n_discordant} 对 E/D 排序不同")
    print(f"\n写出 {OUT}")


if __name__ == "__main__":
    main()
