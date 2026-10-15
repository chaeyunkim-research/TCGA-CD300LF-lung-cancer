# -*- coding: utf-8 -*-
"""
analysis_specific/forest_plot_specific.py
============================================
여러 GEO 코호트의 Cox 결과(cox_regression_specific.py 출력)를 모아
DerSimonian-Laird random-effects 메타분석으로 풀링하고 forest plot을 그린다.

기존 대화에서 만들었던 forest plot 이미지(HR=0.906 (0.700-1.172), p=0.4531, N=713,
코호트별 GSE72094 HR=0.732/GSE31210 HR=1.100/GSE36471 HR=0.994)를 이 스크립트의
공식으로 그대로 재현할 수 있음을 확인했습니다 (아래 랜덤효과 공식으로 역검산 완료).

이 스크립트는 이전에는 별도 파일로 저장되어 있지 않고 대화 중 즉석으로 만들어졌던
부분이라, 이번에 처음으로 재사용 가능한 스크립트로 정식 등록합니다.

입력: cox_regression_specific.py가 만든 "<GENE>_cox_all_cohorts.csv"
      (컬럼: Cohort, HR, CI_lower, CI_upper, p_value, N)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths_config as pcfg  # noqa: E402

# ===== 여기만 바꾸면 됩니다 =====
GENE_CONFIG_MODULE = "gene_config.CD300LF_panel"
# ================================


def random_effects_meta(hr, ci_lower, ci_upper):
    """DerSimonian-Laird random-effects 메타분석. (hr, ci_lower, ci_upper)는 array-like."""
    y = np.log(np.asarray(hr, dtype=float))
    v = ((np.log(np.asarray(ci_upper, dtype=float)) - np.log(np.asarray(ci_lower, dtype=float))) / (2 * 1.96)) ** 2
    w = 1 / v

    mean_fixed = np.sum(w * y) / np.sum(w)
    Q = np.sum(w * (y - mean_fixed) ** 2)
    df = len(y) - 1
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - df) / C) if df > 0 and C > 0 else 0.0

    w_star = 1 / (v + tau2)
    pooled_y = np.sum(w_star * y) / np.sum(w_star)
    se = np.sqrt(1 / np.sum(w_star))

    pooled_hr = np.exp(pooled_y)
    ci_lo = np.exp(pooled_y - 1.96 * se)
    ci_hi = np.exp(pooled_y + 1.96 * se)
    z = pooled_y / se
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return pooled_hr, ci_lo, ci_hi, p


def make_forest_plot(cox_df: pd.DataFrame, target_gene: str, out_dir: Path):
    pooled_hr, pooled_lo, pooled_hi, pooled_p = random_effects_meta(
        cox_df["HR"], cox_df["CI_lower"], cox_df["CI_upper"]
    )
    pooled_n = int(cox_df["N"].sum())

    n_rows = len(cox_df) + 1  # +1 for pooled row
    fig, ax = plt.subplots(figsize=(9, 1.1 * n_rows + 1))

    y_positions = list(range(n_rows, 1, -1))  # cohorts from top, pooled at bottom
    labels = [c.replace("_LUAD", "") for c in cox_df["Cohort"]]  # 표시용 라벨에서는 "_LUAD" 접미사 생략

    for y, (_, row) in zip(y_positions, cox_df.iterrows()):
        ax.plot([row["CI_lower"], row["CI_upper"]], [y, y], color="black", lw=1.5, zorder=1)
        ax.plot([row["CI_lower"], row["CI_lower"]], [y - 0.1, y + 0.1], color="black", lw=1.5)
        ax.plot([row["CI_upper"], row["CI_upper"]], [y - 0.1, y + 0.1], color="black", lw=1.5)
        ax.scatter([row["HR"]], [y], s=140, color="#c0392b", edgecolor="black", zorder=2)
        ax.text(
            2.28, y,
            f"HR={row['HR']:.3f} ({row['CI_lower']:.3f}-{row['CI_upper']:.3f})\n"
            f"p={row['p_value']:.4f}, N={int(row['N'])}",
            va="center", fontsize=10,
        )

    pooled_y = 1
    diamond_half_h = 0.28
    diamond = Polygon(
        [(pooled_lo, pooled_y), (pooled_hr, pooled_y + diamond_half_h),
         (pooled_hi, pooled_y), (pooled_hr, pooled_y - diamond_half_h)],
        closed=True, facecolor="#2c3e50", edgecolor="#2c3e50", zorder=2,
    )
    ax.add_patch(diamond)
    ax.text(
        2.28, pooled_y,
        f"HR={pooled_hr:.3f} ({pooled_lo:.3f}-{pooled_hi:.3f})\np={pooled_p:.4f}, N={pooled_n}",
        va="center", fontsize=10, fontweight="bold",
    )

    ax.axhline(1.6, color="gray", lw=0.8)
    ax.axvline(1.0, color="gray", ls="--", lw=1)

    ax.set_yticks(y_positions + [pooled_y])
    ax.set_yticklabels(labels + ["Random-effect Pooled"])
    for tick, y in zip(ax.get_yticklabels(), y_positions + [pooled_y]):
        if y == pooled_y:
            tick.set_fontweight("bold")

    ax.set_ylim(0.3, n_rows + 1)
    ax.set_xlim(0.5, 2.6)
    ax.set_xlabel(f"Hazard ratio for {target_gene} (95% CI)", fontsize=12)
    ax.set_title(f"External GEO validation of {target_gene} prognostic effect (LUAD)",
                 fontsize=14, fontweight="bold")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(left=False)

    plt.tight_layout()
    out_path = out_dir / f"{target_gene}_forest_plot.png"
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"Pooled HR={pooled_hr:.3f} ({pooled_lo:.3f}-{pooled_hi:.3f}), p={pooled_p:.4f}, N={pooled_n}")
    print(f"저장 완료: {out_path}")


def main():
    import importlib

    gene_config = importlib.import_module(GENE_CONFIG_MODULE)
    target_gene = gene_config.TARGET_GENE

    cox_csv = pcfg.results(target_gene, "cox") / f"{target_gene}_cox_all_cohorts.csv"
    out_dir = pcfg.results(target_gene, "forest_plot")

    cox_df = pd.read_csv(cox_csv)
    make_forest_plot(cox_df, target_gene, out_dir)


if __name__ == "__main__":
    main()
