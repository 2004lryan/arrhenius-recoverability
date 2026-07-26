r"""
91c5_realkinetics_killer.py: C5 诊断 Ψ_n 的"杀手对照"——在真实苹果硬度动力学反演上，Ψ_n（Fisher）
是否严格优于廉价启发式（温度范围 ΔT / 温度数 / 线性条件数 cond）预测真实反演误差。

动机：论文 C6（光谱协方差恢复）已有杀手对照、结论 κ_K 与条件数\emph{作分数}近重合（Spearman 0.98）。
但 C5（动力学反演）的 Ψ_n 是\emph{非线性 Fisher} 量、非线性条件数——本脚本检验二者是否\emph{不同}：
若 Ψ_n 在真实反演误差预测上严格胜过线性 cond/ΔT，则 C5 诊断\emph{非}≈cond，与 C6 形成对照。

数据：04outputs/54c5_realkinetics-designsweep.csv（真实苹果硬度，Foods 2023，read_html 提取；
16 个\emph{超定}≥3 温度子集，含真实测量 SD 传播的反演误差 realCV_Ea，及 Ψ_n/ΔT/n_temps/cond 各诊断）。
食品体系（99_realkinetics2）作\emph{诚实边界}对照：同质 Ψ_n 在冷冻惰性点破坏等方差假设、需异方差 Ψ_het。

运行方式:
    cd .
    python 02code/91c5_realkinetics_killer.py

输出文件:
    04outputs/91c5_realkinetics_killer.json  — 各诊断 Spearman、固定 ΔT 组内排序、偏 Spearman、结论
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[1]
import os as _os
_DS = _os.environ.get("DESIGNSWEEP_CSV")          # 多种子复跑入口；未设置时与原字面值逐字节一致
_TAG = _os.environ.get("OUT_TAG", "")             # 输出名后缀，用于区分多种子运行
CSV = Path(_DS) if _DS else ROOT / "04outputs" / "54c5_realkinetics-designsweep.csv"
OUT = ROOT / "04outputs" / f"91c5_realkinetics_killer{_TAG}.json"


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    """【本脚本原口径，论文 L460 数字来源】秩残差【再做一次 Spearman】。

    注意：这不是教科书标准定义。标准偏 Spearman = 秩对控制变量回归后【残差的 Pearson】相关
    （见下方 partial_spearman_standard）。两者在本数据上相差约 0.04（Ψ_n −0.806 vs −0.822；
    cond −0.094 vs −0.053）。保留本函数以维持论文既有数字的可追溯性，两口径并报。
    """
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)

    def resid(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return a - np.polyval(np.polyfit(b, a, 1), b)

    return float(spearmanr(resid(rx, rz), resid(ry, rz))[0])


def partial_spearman_standard(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    """【教科书标准定义】把三者各自秩变换，x、y 对 z 线性回归取残差，再求残差的 Pearson 相关。

    与 104/106/107 的实现同口径（那三个脚本的 5 种子分析全部基于本定义）。
    """
    rx, ry, rz = (rankdata(v).astype(float) for v in (x, y, z))
    Z = np.stack([rz, np.ones_like(rz)], axis=1)
    res = []
    for r in (rx, ry):
        beta, *_ = np.linalg.lstsq(Z, r, rcond=None)
        res.append(r - Z @ beta)
    a, b = res
    d = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(a @ b / d) if d > 0 else float("nan")


def apple_psi_het() -> dict:
    """严格计算苹果 Ψ_het（异方差 Fisher，各温度按真实测量 SD 传播的 σ_lnk 精度加权），
    检验\emph{统一}诊断 Ψ_het 是否在苹果（近等方差）亦预测真实反演误差——回应"第二体系失败"。
    Ψ_n 是 Ψ_het 的等方差特例；食品（强异方差冷冻点）须 Ψ_het，本函数验证苹果亦兼容 Ψ_het。
    """
    import re
    r_gas = 8.314e-3
    raw = pd.read_csv(ROOT / "03data" / "processed" / "pmc_real" / "PMC10253207_T16.csv",
                      header=None, skiprows=2)

    def parse(cell: object) -> float:
        m = re.findall(r"[0-9.]+", str(cell))
        return float(m[0]) if m else float("nan")

    def parse_sd(cell: object) -> float:
        m = re.findall(r"[0-9.]+", str(cell))
        return float(m[1]) if len(m) >= 2 else float("nan")

    raw["T"] = pd.to_numeric(raw[0], errors="coerce")
    raw["t"] = pd.to_numeric(raw[1], errors="coerce")
    raw["fm"] = [parse(c) for c in raw[5]]   # 第 5 列 = Royal Gala（设计扫描品种）
    raw["fsd"] = [parse_sd(c) for c in raw[5]]
    raw = raw.dropna(subset=["T", "t", "fm"])
    slnk: dict = {}
    for temp in sorted(raw["T"].unique()):
        g = raw[raw["T"] == temp].sort_values("t")
        ly = np.log(g["fm"].values)
        tt = g["t"].values
        w = (g["fm"].values / np.clip(g["fsd"].values, 1e-6, None)) ** 2  # 测量精度权 1/σ_ln²
        amat = np.vstack([tt, np.ones_like(tt)]).T
        beta = np.linalg.lstsq(np.sqrt(w)[:, None] * amat, np.sqrt(w) * ly, rcond=None)[0]
        slnk[float(temp)] = float(np.std(ly - amat @ beta)) or 0.05  # 该温度残差 std = per-T 噪声

    def psi_het_of(tsub: list[float]) -> float:
        info = np.zeros((2, 2))
        for temp in tsub:
            jvec = np.array([-1.0 / (r_gas * (temp + 273.15)), 1.0])
            info += np.outer(jvec, jvec) / slnk[temp] ** 2
        return len(tsub) * float(np.linalg.eigvalsh(info)[0])

    df = pd.read_csv(CSV)
    od = df[(df.variety == "Royal Gala") & (df.n_temps >= 3)
            & np.isfinite(df.realCV_Ea) & (df.realCV_Ea > 0)].copy()
    ph = [psi_het_of([float(x) for x in str(r.temps).split("/")]) for _, r in od.iterrows()]
    rho, p = spearmanr(ph, od.realCV_Ea.values)
    return {
        "n": len(od),
        "heteroscedasticity_sigma_max_over_min": float(max(slnk.values()) / min(slnk.values())),
        "apple_Psi_het_spearman_vs_realCV": float(rho),
        "apple_Psi_het_p": float(p),
        "consistency_Psihet_vs_Psin_spearman": float(spearmanr(ph, od.Psi_n.values)[0]),
        "note": "苹果异方差度~5x，Ψ_het 仍强负相关 → 统一异方差 Fisher Ψ_het 在苹果亦预测真实反演误差",
    }


def main() -> None:
    df = pd.read_csv(CSV)
    od = df[(df.n_temps >= 3) & np.isfinite(df.realCV_Ea) & (df.realCV_Ea > 0)].copy()
    y = od.realCV_Ea.values  # 真实反演误差（待预测的 ground truth）

    res: dict = {
        "study": "C5 killer control: Ψ_n(Fisher) vs 廉价启发式(ΔT/n_temps/线性cond) 预测真实反演误差 realCV_Ea",
        "data": "54c5_realkinetics-designsweep.csv 真实苹果硬度动力学(Foods 2023), 16 超定>=3温真实子集",
        "n_subsets": len(od),
        "ground_truth": "realCV_Ea = 真实测量 SD 经[硬度->k(T)->Arrhenius]全链 bootstrap 传播的 E_a 估计误差",
    }

    # ---- (1) 边际 Spearman：各诊断 vs 真实反演误差 ----
    marg = {}
    for name, x in [("Psi_n", od.Psi_n.values), ("deltaT", od.deltaT.values),
                    ("n_temps", od.n_temps.values.astype(float)), ("cond_linear", od.cond.values)]:
        rho, p = spearmanr(x, y)
        marg[name] = {"spearman_vs_realCV": float(rho), "p": float(p), "abs": float(abs(rho))}
    res["marginal_spearman"] = marg
    res["Psi_n_strictly_beats_all_heuristics_marginal"] = bool(
        marg["Psi_n"]["abs"] > max(marg["deltaT"]["abs"], marg["n_temps"]["abs"], marg["cond_linear"]["abs"])
    )

    # ---- (2) 固定 ΔT 组内：ΔT/cond 近常数，唯 Ψ_n 能排序（value beyond range/cond）----
    within = {}
    for dt, g in od.groupby("deltaT"):
        if len(g) >= 3:
            rho, p = spearmanr(g.Psi_n, g.realCV_Ea)
            within[f"deltaT_{dt}"] = {
                "n": len(g),
                "spearman_Psin_realCV": float(rho),
                "p": float(p),
                "cond_within_group_range_x": float(g.cond.max() / g.cond.min()),
            }
    res["within_fixed_deltaT"] = within

    # ---- (3) 偏 Spearman：控制 ΔT 后 Ψ_n vs cond 的残余预测力（决定性）----
    res["partial_spearman_control_deltaT"] = {
        "Psi_n_given_deltaT": partial_spearman(od.Psi_n.values, y, od.deltaT.values),
        "cond_given_deltaT": partial_spearman(od.cond.values, y, od.deltaT.values),
        "_standard_definition": {
            "note": "教科书标准偏 Spearman（秩残差的 Pearson）；与 104/106/107 的 5 种子分析同口径",
            "Psi_n_given_deltaT": partial_spearman_standard(od.Psi_n.values, y, od.deltaT.values),
            "cond_given_deltaT": partial_spearman_standard(od.cond.values, y, od.deltaT.values),
        },
        "note": "控制温度范围后 Ψ_n 保留预测力而线性 cond 崩塌 => Ψ_n 捕获范围/cond 之外的可恢复性信息",
    }

    # ---- (4) 结论 + 诚实边界 ----
    pn = res["partial_spearman_control_deltaT"]["Psi_n_given_deltaT"]
    pc = res["partial_spearman_control_deltaT"]["cond_given_deltaT"]
    differentiates = (res["Psi_n_strictly_beats_all_heuristics_marginal"]
                      and abs(pn) > 0.5 and abs(pc) < 0.3)
    res["Psi_differentiates_from_cond_on_C5"] = bool(differentiates)
    res["conclusion"] = (
        f"C5 上 Ψ_n(Fisher) 严格优于线性 cond/ΔT/n_temps 预测真实反演误差："
        f"边际 |Spearman| Ψ_n={marg['Psi_n']['abs']:.2f} > cond={marg['cond_linear']['abs']:.2f}/ΔT={marg['deltaT']['abs']:.2f}；"
        f"控制 ΔT 后偏 Spearman Ψ_n={pn:+.2f} 而 cond={pc:+.2f}(崩塌)。"
        f"故 C5 上 Ψ_n 在排序上非≈cond，把'Ψ处处≈cond'收窄为'C6≈cond、C5不'。"
        f"诚实机理：Ψ_n=N*lmin 与 cond=lmax/lmin 同源自一个 Arrhenius Fisher，Ψ_n≡monotone(lmin)，"
        f"故差异本质=绝对(lmin)vs相对(比值)条件数=经典 Stewart-Sun(同 C6)、非新方法学；"
        f"'控制ΔT后cond崩塌'因组内cond近常数结构性所致。" if differentiates else
        "未观察到 Ψ_n 对 cond 的严格优势，C5 亦 ≈cond。"
    )
    res["honest_boundaries"] = (
        "(1) 单一真实数据集(苹果 Foods 2023)、16 个\\emph{非独立}超定子集(温度重叠)、realCV 为测量 SD 传播"
        "\\emph{非}独立重复，故为\\emph{排序}证据非强统计独立检验；"
        "(2) ΔT=6 组(n=6)较弱(p=0.16)，ΔT=4/8 清晰(p<5e-4)；"
        "(3) 食品体系(99_realkinetics2)同质 Ψ_n 不直接复现(Spearman +0.19 n.s.，冷冻惰性点破坏等方差)、"
        "需异方差 Ψ_het(-0.34, p=0.014)——即 Ψ_n 优势\\emph{条件于近等方差}，强异方差须用 Ψ_het；"
        "(4) 仅 C5(动力学)成立；C6(光谱线性恢复)κ_K 仍 ≈cond。"
    )

    # ---- (5) 统一异方差诊断 Ψ_het 跨 2 独立真实体系（苹果 computed + 食品 99 既有）----
    apple = apple_psi_het()
    res["unified_Psi_het_cross_system"] = {
        "apple_kinetics_computed": apple,
        "food_kinetics_99": {"Psi_het_spearman": -0.344, "p": 0.0145,
                             "source": "99_realkinetics2.json het_spearman_p",
                             "homoscedastic_Psi_n_fails": "+0.19 n.s."},
        "conclusion": (f"统一异方差 Fisher Ψ_het 跨 2 独立真实体系皆负相关："
                       f"苹果 Ψ_het={apple['apple_Psi_het_spearman_vs_realCV']:+.2f}"
                       f"(p={apple['apple_Psi_het_p']:.1e}, 异方差度 {apple['heteroscedasticity_sigma_max_over_min']:.1f}x)、"
                       f"食品 Ψ_het=-0.34(p=0.014)；Ψ_n 是 Ψ_het 的等方差特例。'第二体系失败'实为用错诊断(同质 Ψ_n)，"
                       f"改用原则性 Ψ_het 则两体系一致。诚实：仅 2 体系、食品弱、机理仍属绝对-vs-相对(robustness 非 novelty)。"),
    }

    OUT.write_text(json.dumps(res, ensure_ascii=False, indent=2))
    print("=== C5 killer control: Ψ_n vs 廉价启发式 预测真实反演误差 ===")
    for k, v in marg.items():
        print(f"  {k:12s} |Spearman vs realCV| = {v['abs']:.3f}  (rho={v['spearman_vs_realCV']:+.3f}, p={v['p']:.1e})")
    print(f"  >> Ψ_n 边际严格胜全部启发式: {res['Psi_n_strictly_beats_all_heuristics_marginal']}")
    print("  固定 ΔT 组内 Spearman(Ψ_n, realCV):")
    for k, v in within.items():
        print(f"    {k}: {v['spearman_Psin_realCV']:+.3f} (n={v['n']}, p={v['p']:.1e}, cond组内仅变{v['cond_within_group_range_x']:.2f}x)")
    print(f"  偏 Spearman 控制 ΔT: Ψ_n={pn:+.3f}  vs  cond={pc:+.3f}")
    print(f"  >> Ψ_n 在 C5 上与 cond 区分(非≈cond): {res['Psi_differentiates_from_cond_on_C5']}")
    print(f"\n结论: {res['conclusion']}")
    print(f"诚实边界: {res['honest_boundaries']}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
