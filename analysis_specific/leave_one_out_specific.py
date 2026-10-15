# -*- coding: utf-8 -*-
"""
analysis_specific/leave_one_out_specific.py
============================================
forest_plot_specific.py가 만든 "<GENE>_cox_all_cohorts.csv"를 입력으로 받아,
코호트를 한 번에 하나씩 빼고 나머지로 다시 random-effects 메타분석을 돌리는
leave-one-out sensitivity analysis.

목적: pooled 유의성이 특정 코호트 하나에 과도하게 의존하고 있는지 확인.
(예: GSE72094 하나를 빼면 유의성이 사라지는지 등)

입력: cox_regression_specific.py가 만든 "<GENE>_cox_all_cohorts.csv"
출력:
  - results/<GENE>/leave_one_out/<GENE>_leave_one_out.csv (코호트별 결과 표)
  - results/<GENE>/leave_one_out/<GENE>_leave_one_out_plot.png (sensitivity forest plot)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths_config as pcfg  # noqa: E402
from forest_plot_specific import random_effects_meta  # noqa: E402

# ===== 여기만 바꾸면 됩니다 =====
GENE_CONFIG_MODULE = "gene_config.CD300LF_panel"
# ================================


def run_leave_one_out(cox_df: pd.DataFrame, target_gene: str, out_dir: Path):
    all_hr, all_lo, all_hi, all_p = random_effects_meta(
        cox_df["HR"], cox_df["CI_lower"], cox_df["CI_upper"]
    )
    all_n = int(cox_df["N"].sum())

    rows = []
    for excl_cohort in cox_df["Cohort"]:
        subset = cox_df[cox_df["Cohort"] != excl_cohort]
        hr, lo, hi, p = random_effects_meta(subset["HR"], subset["CI_lower"], subset["CI_upper"])
        rows.append({
            "Excluded_Cohort": excl_cohort,
            "Pooled_HR": round(hr, 3),
            "CI_lower": round(lo, 3),
            "CI_upper": round(hi, 3),
            "p_value": round(p, 4),
            "N_remaining": int(subset["N"].sum()),
            "Significant_at_0.05": p < 0.05,
        })

    loo_df = pd.DataFrame(rows)
    csv_path = out_dir / f"{target_gene}_leave_one_out.csv"
    loo_df.to_csv(csv_path, index=False)
    print(f"=== {target_gene} Leave-one-out sensitivity analysis ===")
    print(f"전체 {len(cox_df)}개 코호트 pooled: HR={all_hr:.3f} ({all_lo:.3f}-{all_hi:.3f}), "
          f"p={all_p:.4f}, N={all_n}\n")
    print(loo_df.to_string(index=False))
    print(f"\n저장 완료: {csv_path}")

    # --- sensitivity forest plot ---
    n_rows = len(loo_df) + 1  # +1 for "전체(All cohorts)" 기준선
    fig, ax = plt.subplots(figsize=(9, 0.9 * n_rows + 1.2))

    y_positions = list(range(n_rows, 1, -1))
    labels = [f"Excluding {c.replace('_LUAD', '')}" for c in loo_df["Excluded_Cohort"]]

    # 표준 forest plot 관례: 개별 행(코호트 제외 시나리오)에는 유의성 표시(*, n.s.)나
    # 색 구분을 넣지 않습니다. 신뢰구간이 HR=1(세로 점선) 을 가로지르는지가 유의성을
    # 보여주는 시각적 단서이고, 색/마커 구분은 pooled(전체) 추정치의 다이아몬드에만
    # 씁니다 (metafor, RevMan, Stata metaninf 등에서 공통적으로 쓰는 방식).
    row_color = "#34495e"
    for y, (_, row) in zip(y_positions, loo_df.iterrows()):
        ax.plot([row["CI_lower"], row["CI_upper"]], [y, y], color=row_color, lw=1.5, zorder=1)
        ax.plot([row["CI_lower"], row["CI_lower"]], [y - 0.1, y + 0.1], color=row_color, lw=1.5)
        ax.plot([row["CI_upper"], row["CI_upper"]], [y - 0.1, y + 0.1], color=row_color, lw=1.5)
        ax.scatter([row["Pooled_HR"]], [y], s=120, color=row_color, edgecolor="black", zorder=2)
        ax.text(
            2.05, y,
            f"HR={row['Pooled_HR']:.3f} ({row['CI_lower']:.3f}-{row['CI_upper']:.3f})\n"
            f"p={row['p_value']:.4f}, N={row['N_remaining']}",
            va="center", fontsize=9.5,
        )

    all_y = 1
    ax.plot([all_lo, all_hi], [all_y, all_y], color="#2c3e50", lw=2.2, zorder=1)
    ax.plot([all_lo, all_lo], [all_y - 0.12, all_y + 0.12], color="#2c3e50", lw=2.2)
    ax.plot([all_hi, all_hi], [all_y - 0.12, all_y + 0.12], color="#2c3e50", lw=2.2)
    ax.scatter([all_hr], [all_y], s=160, marker="D", color="#2c3e50", edgecolor="black", zorder=2)
    ax.text(
        2.05, all_y,
        f"HR={all_hr:.3f} ({all_lo:.3f}-{all_hi:.3f})\np={all_p:.4f}, N={all_n}",
        va="center", fontsize=9.5, fontweight="bold",
    )

    ax.axhline(1.6, color="gray", lw=0.8)
    ax.axvline(1.0, color="gray", ls="--", lw=1)

    ax.set_yticks(y_positions + [all_y])
    ax.set_yticklabels(labels + ["All 5 cohorts (reference)"])
    for tick, y in zip(ax.get_yticklabels(), y_positions + [all_y]):
        if y == all_y:
            tick.set_fontweight("bold")

    ax.set_ylim(0.3, n_rows + 1)
    ax.set_xlim(0.5, 2.3)
    ax.set_xlabel(f"Pooled hazard ratio for {target_gene} (95% CI)", fontsize=12)
    ax.set_title(f"Leave-one-out sensitivity analysis ({target_gene}, LUAD)",
                 fontsize=13, fontweight="bold")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(left=False)

    plt.tight_layout()
    plot_path = out_dir / f"{target_gene}_leave_one_out_plot.png"
    plt.savefig(plot_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"저장 완료: {plot_path}")

    return loo_df


def main():
    import importlib

    gene_config = importlib.import_module(GENE_CONFIG_MODULE)
    target_gene = gene_config.TARGET_GENE

    cox_csv = pcfg.results(target_gene, "cox") / f"{target_gene}_cox_all_cohorts.csv"
    out_dir = pcfg.results(target_gene, "leave_one_out")

    cox_df = pd.read_csv(cox_csv)
    run_leave_one_out(cox_df, target_gene, out_dir)


if __name__ == "__main__":
    main()
