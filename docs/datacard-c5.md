# 0NDATACARD_C5.md — C5 稿件数据卡（SHA-256 登记）

生成时间：2026-07-26 04:17:49 -0400 · 生成脚本：`02code/120datacard_sha256.py`

> 依据 `01CLAUDE.md` §4.6 / §202。本卡只登记 **C5 稿件（`Elsevier_zh.tex` / `Elsevier_en.tex`）实际依赖**的文件，
> 不含 C6/光谱方向的历史产物（`03data/processed/` 下 35 项与本稿无关）。

> **诚实边界**：SHA-256 仅证明**自本次登记起**文件字节未变，不追认登记之前的历史，
> 也不构成对提取步骤正确性的证明——提取正确性由脚本 55 的 `read_html` 与人工/Europe PMC 逐条核验承担。

---

## 一、源数据

### apple_firmness — 苹果硬度温度梯度（4 品种 × 5 贮温 × 5 时间点）
- **角色**：C5 真实动力学 / 体系 1
- **来源**：Foods 2023, doi:10.3390/foods12112113 (PMC10253207) 表 A1（开放获取）
- **获取/提取方式**：pandas.read_html 提取（脚本 55）+ 人工逐条核验；非本组采集
- **许可**：CC BY 4.0（Foods 为全开放获取）
- **在盘路径**：`03data/processed/pmc_real/PMC10253207_T16.csv`
- **大小**：1,770 bytes
- **SHA-256**：`d22131f037a120e71a321d358e141ea0b8b2b8b97c4299706f4b0c8fad66665f`

### food_ammonia_tbars — 长货架期食品 氨/TBARS 累积动力学（5 食品 × 4 贮温，10 序列）
- **角色**：C5 真实动力学 / 体系 2
- **来源**：Foods 2022 11(14):2004, doi:10.3390/foods11142004 (PMC9319022) 表 3
- **获取/提取方式**：pandas.read_html 提取（脚本 55）+ Europe PMC fullTextXML 独立重抓逐值核对（132 行 0 处不一致）
- **许可**：CC BY 4.0
- **在盘路径**：`03data/processed/pmc_real/PMC9319022_T3.csv`
- **大小**：5,239 bytes
- **SHA-256**：`ee6da3d22db941fe2936bba9817e3aa7873e408024ae9ab95d6bf48d536ccfbe`

### hass_avocado — Hass 牛油果采后成熟（478 果 × 3 组贮藏，14,722 条记录）
- **角色**：C5 真实动力学 / 体系 3（真统计估计误差）
- **来源**：Mendeley Data, doi:10.17632/3xd9n945v8.1
- **获取/提取方式**：原始 xlsx 直接读取（未做 L1 变换）；按 4.5 只读引用，未复制进项目
- **许可**：CC BY 4.0（Mendeley Data 默认）
- **在盘路径**：`<EXTERNAL_DATA_ROOT>/004_hass_avocado_rgb_ripening/Hass Avocado Ripening Photographic Dataset/Avocado Ripening Dataset.xlsx`
- **大小**：633,351 bytes
- **SHA-256**：`f8abeaba6eedf67869907b3a71ee132cd89b42bb3dcce6b25ddb22e6269b51a7`

---

## 二、派生产物（稿件每个数字的来源）

| 产物 | 内容 | bytes | SHA-256 (前 16) |
| --- | --- | ---: | --- |
| `43trichotomy_formal.json` | 半合成三分律：Ψ 终值、塌缩斜率、Krug 匹配 | 1,317 | `3ddee85bf7b57f03` |
| `54c5_realkinetics.json` | 苹果体系：品种 Ea/R²、设计扫描、超定子集 | 5,987 | `9ba2ab597b83d50a` |
| `54c5_realkinetics-designsweep.csv` | 苹果 26 温度子集逐条（补充表 S2 数据源） | 4,115 | `2b1769aa0b095d56` |
| `91c5_realkinetics_killer.json` | Ψ vs 廉价启发式判别力、Ψ_het 跨体系 | 4,103 | `64fc189dacbb29c7` |
| `99_realkinetics2.json` | 食品体系：Ea 范围、同质/异质 Ψ 复现 | 14,092 | `7b541f1778481b00` |
| `102_avocado_c5_pilot.json` | 牛油果体系：真统计估计误差、n 轴扫描 | 2,519 | `125ba608e1256276` |
| `119pazman_pronzato_bridge.json` | 与 Pázman–Pronzato 2006 算例的坐标对照 | 2,236 | `6d5f827620897f26` |
| `122criterion_bakeoff.json` | 准则对决（含对自身假设的证伪：塌缩=代数、D 对比未建立） | 4,854 | `273a0263ec9b06b8` |
| `123montecarlo_criterion_study.json` | 1000 独立设计 Monte-Carlo：两轴速率、饱和、决策 regret | 5,924 | `ebaea2757742c77a` |
| `123montecarlo_criterion_study-designs.csv` | 逐设计表：准则值 + 真值 RMSE（正文 §sec:mc 数据源） | 335,166 | `a93b4a8da05dd96a` |
| `124coverage_study.json` | 置信区间覆盖率研究：Ψ 作为标定诊断被证伪 + 区间口径警告 | 7,492 | `ec75a190effb964a` |
| `124coverage_study-designs.csv` | 逐设计覆盖率表 | 264,231 | `2b5afdae281cd949` |

完整 64 位摘要见 `06doc/0NDATACARD_C5_manifest.tsv`。

---

## 三、与稿件的对应
- 稿件表 `tab:data`（数据卡表）列出的三个真实体系即本卡第一节三项。
- 稿件正文与补充材料的每个数字，其来源产物见本卡第二节，并在 `06doc/02sub/claim_evidence_map.md` 中逐条绑定。

## 四、缺失项

无。
