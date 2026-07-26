#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
120 — C5 稿件数据卡与 SHA-256 登记（补齐 01CLAUDE.md §4.6 / §202 的强制登记）

背景
----
2026-07-25 核查发现：本项目**从未**对 C5 稿件所用数据做过 SHA-256 登记——
`.datalock-manifest.tsv` 不存在，`08github/data/DATA.md` 仍是含 {{占位符}} 的模板，
`08github/docs/` 无任何数据卡，而 `03data/processed/*-meta.json` 不含 SHA 字段。
但稿件表 tab:data 的表注声称各数据集"经确定性 L1 预处理并附 SHA-256 校验"。
本脚本把该声明做实：对稿件实际依赖的每个文件计算 SHA-256 并落盘登记。

诚实注记（2026-07-26 修正）
--------------------------------
本清单最初写在 `03data/.datalock-manifest.tsv`。那是错的：该路径由 `datalock.sh`
独占，其 schema 为 `sha256⇥size⇥ts⇥相对 03data 的路径`，而本脚本写的是
`kind⇥id⇥sha256⇥bytes⇥绝对路径`。撞名的后果是 `datalock.sh verify` 把本清单
15 条连同表头全判为 MISSING，`checkpoint.sh` 因「数据完整性校验失败」而永久拒绝提交。
本清单的覆盖范围本来就与 datalock 不同（跨 03data / 04outputs / 只读的跨项目共享库），
无法用相对 DATA_DIR 的路径表达，故改放到数据卡旁边，两套登记各管各的。

覆盖范围（只登记 C5 稿件真正用到的文件，不含 C6 遗留产物）
----------------------------------------------------------
A. 源数据（他组已发表表格的提取件 + 第三体系原始数据）
B. 派生产物（正文与补充材料每个数字的来源 JSON/CSV）

输出
----
- `06doc/0NDATACARD_C5_manifest.tsv`  逐文件 SHA-256 清单（TSV）
- `06doc/0NDATACARD_C5.md`         中文数据卡（§1 项目内部工作区中文强制）
- `04outputs/120datacard_sha256.json`

诚实边界
--------
SHA-256 仅证明**自登记之时起**文件字节未变；它不追认登记之前的历史，也不构成
对提取步骤本身正确性的证明（提取正确性另由脚本 55 的 read_html + 人工逐条核验承担）。
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW_ROOT = Path("<EXTERNAL_DATA_ROOT>")

# --- A. 源数据 -------------------------------------------------------------
SOURCES = [
    {
        "id": "apple_firmness",
        "name": "苹果硬度温度梯度（4 品种 × 5 贮温 × 5 时间点）",
        "role": "C5 真实动力学 / 体系 1",
        "path": ROOT / "03data/processed/pmc_real/PMC10253207_T16.csv",
        "origin": "Foods 2023, doi:10.3390/foods12112113 (PMC10253207) 表 A1（开放获取）",
        "extraction": "pandas.read_html 提取（脚本 55）+ 人工逐条核验；非本组采集",
        "license": "CC BY 4.0（Foods 为全开放获取）",
    },
    {
        "id": "food_ammonia_tbars",
        "name": "长货架期食品 氨/TBARS 累积动力学（5 食品 × 4 贮温，10 序列）",
        "role": "C5 真实动力学 / 体系 2",
        "path": ROOT / "03data/processed/pmc_real/PMC9319022_T3.csv",
        "origin": "Foods 2022 11(14):2004, doi:10.3390/foods11142004 (PMC9319022) 表 3",
        "extraction": "pandas.read_html 提取（脚本 55）+ Europe PMC fullTextXML 独立重抓逐值核对（132 行 0 处不一致）",
        "license": "CC BY 4.0",
    },
    {
        "id": "hass_avocado",
        "name": "Hass 牛油果采后成熟（478 果 × 3 组贮藏，14,722 条记录）",
        "role": "C5 真实动力学 / 体系 3（真统计估计误差）",
        "path": RAW_ROOT / "004_hass_avocado_rgb_ripening/Hass Avocado Ripening Photographic Dataset/Avocado Ripening Dataset.xlsx",
        "origin": "Mendeley Data, doi:10.17632/3xd9n945v8.1",
        "extraction": "原始 xlsx 直接读取（未做 L1 变换）；按 4.5 只读引用，未复制进项目",
        "license": "CC BY 4.0（Mendeley Data 默认）",
    },
]

# --- B. 派生产物（稿件每个数字的来源） -------------------------------------
DERIVED = [
    ("04outputs/43trichotomy_formal.json", "半合成三分律：Ψ 终值、塌缩斜率、Krug 匹配"),
    ("04outputs/54c5_realkinetics.json", "苹果体系：品种 Ea/R²、设计扫描、超定子集"),
    ("04outputs/54c5_realkinetics-designsweep.csv", "苹果 26 温度子集逐条（补充表 S2 数据源）"),
    ("04outputs/91c5_realkinetics_killer.json", "Ψ vs 廉价启发式判别力、Ψ_het 跨体系"),
    ("04outputs/99_realkinetics2.json", "食品体系：Ea 范围、同质/异质 Ψ 复现"),
    ("04outputs/102_avocado_c5_pilot.json", "牛油果体系：真统计估计误差、n 轴扫描"),
    ("04outputs/119pazman_pronzato_bridge.json", "与 Pázman–Pronzato 2006 算例的坐标对照"),
    ("04outputs/122criterion_bakeoff.json", "准则对决（含对自身假设的证伪：塌缩=代数、D 对比未建立）"),
    ("04outputs/123montecarlo_criterion_study.json", "1000 独立设计 Monte-Carlo：两轴速率、饱和、决策 regret"),
    ("04outputs/123montecarlo_criterion_study-designs.csv", "逐设计表：准则值 + 真值 RMSE（正文 §sec:mc 数据源）"),
    ("04outputs/124coverage_study.json", "置信区间覆盖率研究：Ψ 作为标定诊断被证伪 + 区间口径警告"),
    ("04outputs/124coverage_study-designs.csv", "逐设计覆盖率表"),
]


def sha256(p: Path) -> tuple[str, int]:
    h = hashlib.sha256()
    n = 0
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
            n += len(chunk)
    return h.hexdigest(), n


def main() -> None:
    stamp = datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%d %H:%M:%S %z")
    rows, missing = [], []

    for s in SOURCES:
        p = Path(s["path"])
        if not p.is_file():
            missing.append(str(p))
            continue
        d, n = sha256(p)
        rows.append({"kind": "source", "id": s["id"], "path": str(p), "sha256": d, "bytes": n, **{k: s[k] for k in ("name", "role", "origin", "extraction", "license")}})

    for rel, desc in DERIVED:
        p = ROOT / rel
        if not p.is_file():
            missing.append(str(p))
            continue
        d, n = sha256(p)
        rows.append({"kind": "derived", "id": Path(rel).stem, "path": str(p), "sha256": d, "bytes": n, "name": desc})

    # 1) manifest TSV
    man = ROOT / "06doc/0NDATACARD_C5_manifest.tsv"
    with man.open("w", encoding="utf-8") as f:
        f.write("kind\tid\tsha256\tbytes\tpath\n")
        for r in rows:
            f.write(f"{r['kind']}\t{r['id']}\t{r['sha256']}\t{r['bytes']}\t{r['path']}\n")

    # 2) 中文数据卡（§1 项目内部工作区中文强制）
    card = ROOT / "06doc/0NDATACARD_C5.md"
    L = [
        "# 0NDATACARD_C5.md — C5 稿件数据卡（SHA-256 登记）",
        f"\n生成时间：{stamp} · 生成脚本：`02code/120datacard_sha256.py`",
        "\n> 依据 `01CLAUDE.md` §4.6 / §202。本卡只登记 **C5 稿件（`Elsevier_zh.tex` / `Elsevier_en.tex`）实际依赖**的文件，",
        "> 不含 C6/光谱方向的历史产物（`03data/processed/` 下 35 项与本稿无关）。",
        "\n> **诚实边界**：SHA-256 仅证明**自本次登记起**文件字节未变，不追认登记之前的历史，",
        "> 也不构成对提取步骤正确性的证明——提取正确性由脚本 55 的 `read_html` 与人工/Europe PMC 逐条核验承担。",
        "\n---\n\n## 一、源数据",
    ]
    for r in [x for x in rows if x["kind"] == "source"]:
        L += [
            f"\n### {r['id']} — {r['name']}",
            f"- **角色**：{r['role']}",
            f"- **来源**：{r['origin']}",
            f"- **获取/提取方式**：{r['extraction']}",
            f"- **许可**：{r['license']}",
            f"- **在盘路径**：`{r['path']}`",
            f"- **大小**：{r['bytes']:,} bytes",
            f"- **SHA-256**：`{r['sha256']}`",
        ]
    L += ["\n---\n\n## 二、派生产物（稿件每个数字的来源）", "",
          "| 产物 | 内容 | bytes | SHA-256 (前 16) |", "| --- | --- | ---: | --- |"]
    for r in [x for x in rows if x["kind"] == "derived"]:
        L.append(f"| `{Path(r['path']).name}` | {r['name']} | {r['bytes']:,} | `{r['sha256'][:16]}` |")
    L += [
        "\n完整 64 位摘要见 `06doc/0NDATACARD_C5_manifest.tsv`。",
        "\n---\n\n## 三、与稿件的对应",
        "- 稿件表 `tab:data`（数据卡表）列出的三个真实体系即本卡第一节三项。",
        "- 稿件正文与补充材料的每个数字，其来源产物见本卡第二节，并在 `06doc/02sub/claim_evidence_map.md` 中逐条绑定。",
        f"\n## 四、缺失项\n\n{'无。' if not missing else chr(10).join('- 未找到：`' + m + '`' for m in missing)}",
    ]
    card.write_text("\n".join(L) + "\n", encoding="utf-8")

    out = ROOT / "04outputs/120datacard_sha256.json"
    out.write_text(json.dumps({
        "script": Path(__file__).name, "generated": stamp,
        "n_source": sum(1 for r in rows if r["kind"] == "source"),
        "n_derived": sum(1 for r in rows if r["kind"] == "derived"),
        "missing": missing, "rows": rows,
    }, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"登记 {len(rows)} 个文件（源 {sum(1 for r in rows if r['kind']=='source')} / 派生 {sum(1 for r in rows if r['kind']=='derived')}）")
    for r in rows:
        print(f"  [{r['kind']:7s}] {r['sha256'][:16]}  {r['bytes']:>10,}  {Path(r['path']).name}")
    print(f"缺失: {missing or '无'}")
    print(f"写出 {man}\n写出 {card}\n写出 {out}")


if __name__ == "__main__":
    main()
