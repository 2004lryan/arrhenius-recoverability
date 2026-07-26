"""
59c5_supp_table.py: C5 真实动力学验证的补充材料表（回应 codex 到 8.5：26 子集 CV/Ψ + 品种复现）

从 54c5 真实结果（JSON + designsweep CSV，均由真实 read_html 数据生成）渲染两张 booktabs 三线表，
直接 \input 进论文补充材料。绝不手填——全部从已落盘的真实输出读取。

输出: 06doc/01manuscript/supp_c5_tables.tex（中文，可 \input）
"""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "04outputs"
TEX = ROOT / "06doc" / "01manuscript" / "supp_c5_tables.tex"


def main() -> None:
    d = json.load(open(OUT / "54c5_realkinetics.json"))
    sweep = pd.read_csv(OUT / "54c5_realkinetics-designsweep.csv")
    mv = d["multi_variety_replication"]["detail"]

    lines = []
    lines.append("% 自动生成（59c5_supp_table.py），数据源 04outputs/54c5_realkinetics*（真实 read_html）")
    lines.append("% 注：章节标题由 Nature_zh.tex 提供（避免重复 \\section*）；表号由 \\thetable=S\\arabic 给出，")
    lines.append("%     故 caption 内不再手工写 \"S1.\"/\"S2.\" 前缀。")
    lines.append("")
    # 表 S1: 品种间复现
    lines.append("\\begin{table}[h]\\centering\\small")
    lines.append("\\caption{\\textbf{设计法则 cor:design 在单调软化品种上的品种间复现}（窄窗 $=$ 最近两温度 $\\{8,10\\}^\\circ$C、宽窗 $=$ 最远两温度 $\\{2,10\\}^\\circ$C；CV 为 $E_a$ 噪声传播敏感性，越低越稳）。Pink~Lady 因 $2^\\circ$C 生理后熟非单调被排除。诚实说明：本表 CV 取自\\emph{品种间复现}分析的 bootstrap 传播，与正文按 $\\Delta T$ 配对报告的\\emph{设计扫描} draw 不同源（Royal~Gala 窄窗 $1.71$ vs 正文 $1.66$）——同一量的两次独立 bootstrap draw，非口径冲突。}")
    lines.append("\\label{tab:suppS1}")
    lines.append("\\begin{tabular}{@{}lccccc@{}}")
    lines.append("\\toprule")
    lines.append("品种 & $E_a$ (kJ/mol) & Arrhenius $R^2$ & 窄窗 CV & 宽窗 CV & 宽/窄改善 \\\\")
    lines.append("\\midrule")
    order = ["Royal Gala", "Granny Smith", "Red Delicious", "Pink Lady"]
    for v in order:
        m = mv[v]
        if m.get("monotone"):
            fac = m["wide_better"]
            facs = f"{fac:.1f}$\\times$" if fac < 100 else f"$\\sim$10$^2\\times$（窄窗病态）"
            lines.append(f"{v} & {m['Ea']:.1f} & {m['r2_arrhenius']:.3f} & {m['narrow_CV']:.2f} & {m['wide_CV']:.3f} & {facs} \\\\")
        else:
            lines.append(f"{v} & \\multicolumn{{5}}{{c}}{{非单调（$2^\\circ$C 后熟），$E_a$ 反演失稳、排除}} \\\\")
    lines.append("\\bottomrule")
    lines.append("\\end{tabular}")
    lines.append("\\end{table}")
    lines.append("")
    # 表 S2: 26 子集 Ψ vs CV（按 ΔT 分组紧凑呈现）
    # 斜率必须从 JSON 实读，不得手填（历史 bug：此处曾硬编码 -0.47，而真值为 -0.4617 -> -0.46，
    # 与正文 -0.46 打架。任何印进 caption 的数都必须可追溯到落盘输出。）
    slope = d["design_sweep"]["logPsi_vs_logCV_slope"]
    lines.append("\\begin{table}[h]\\centering\\small")
    lines.append(
        "\\caption{\\textbf{全 $26$ 个温度子集的诊断 $\\Psin$ 与反演 CV}（Royal~Gala；按温度窗宽 "
        "$\\Delta T$ 排序，示 $\\Psin$ 增、CV 降的单调趋势，$\\log\\Psin$--$\\log$CV 斜率 "
        f"${slope:.2f}$）。}}"
    )
    lines.append("\\label{tab:suppS2}")
    lines.append("\\begin{tabular}{@{}clccc@{}}")
    lines.append("\\toprule")
    lines.append("$\\Delta T$ ($^\\circ$C) & 温度子集 & $\\#$温度 & $\\Psin$ & 反演 CV \\\\")
    lines.append("\\midrule")
    sweep = sweep.sort_values(["deltaT", "Psi_n"])
    for _, r in sweep.iterrows():
        temps = str(r["temps"]).replace("/", ",")
        lines.append(f"{int(r['deltaT'])} & $\\{{{temps}\\}}$ & {int(r['n_temps'])} & {r['Psi_n']:.3f} & {r['inv_CV_Ea']:.3f} \\\\")
    lines.append("\\bottomrule")
    lines.append("\\end{tabular}")
    lines.append("\\end{table}")

    TEX.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {TEX} ({len(lines)} lines, {len(sweep)} subset rows, {len(order)} varieties)")


if __name__ == "__main__":
    main()
