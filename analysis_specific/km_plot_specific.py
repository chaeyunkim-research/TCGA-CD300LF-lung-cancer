# -*- coding: utf-8 -*-
"""
analysis_specific/km_plot_specific.py
========================================
TARGET_GENE 중앙값 분리(median split) Kaplan-Meier plot + number-at-risk 테이블.

기존 파일과의 대응: km_plot.py, km_plot_GSE31210.py, km_plot_GSE72094.py,
km_plot_styled.py, km_plot_styled_v2.py
-> 전부 GSE31210 버전(number-at-risk 테이블 포함, x축 20개월 단위)이 가장 완성된
형태였으므로 그걸 표준으로 삼아 COHORT_NAME/TARGET_GENE만 바꾸면 되게 통합.
(styled.py / styled_v2.py는 GSE31210 버전으로 가는 중간 시행착오였습니다.)

여러 코호트를 한 번에 뽑으려면 COHORT_RUNS 리스트에 추가하세요.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths_config as pcfg  # noqa: E402

# ===== 여기만 바꾸면 됩니다 =====
GENE_CONFIG_MODULE = "gene_config.CD300LF_panel"

COHORT_RUNS = [
    "GSE31210",
    "GSE72094",
    "GSE36471",
    "GSE13213",
    "GSE42127",
    # TCGA를 추가하려면 tcga/TCGA_runner_com.py 실행 후 아래 주석 해제:
    # "TCGA_LUAD",
]

COLOR_HIGH = "#F8766D"
COLOR_LOW = "#00BFC4"
# ================================


def make_km_plot(cohort_name: str, merged_valid_path: str, target_gene: str, out_dir: Path):
    merged = pd.read_csv(merged_valid_path)
    df = merged[[target_gene, "survival_months", "survival_status"]].dropna()

    median_val = df[target_gene].median()
    df = df.copy()
    df["group"] = df[target_gene].apply(lambda x: "High" if x >= median_val else "Low")

    high = df[df["group"] == "High"]
    low = df[df["group"] == "Low"]

    result = logrank_test(
        high["survival_months"], low["survival_months"],
        event_observed_A=high["survival_status"], event_observed_B=low["survival_status"],
    )
    p_value = result.p_value

    kmf_high, kmf_low = KaplanMeierFitter(), KaplanMeierFitter()
    kmf_high.fit(high["survival_months"], event_observed=high["survival_status"], label="High")
    kmf_low.fit(low["survival_months"], event_observed=low["survival_status"], label="Low")

    fig, (ax_km, ax_risk) = plt.subplots(
        2, 1, figsize=(7.5, 8), gridspec_kw={"height_ratios": [4, 1], "hspace": 0.05}
    )

    kmf_high.plot_survival_function(ax=ax_km, color=COLOR_HIGH, ci_show=False, show_censors=True)
    kmf_low.plot_survival_function(ax=ax_km, color=COLOR_LOW, ci_show=False, show_censors=True)

    ax_km.get_legend().remove()
    ax_km.set_xlabel("")
    ax_km.set_ylabel("Survival probability", fontsize=11)
    ax_km.set_ylim(0, 1.02)
    ax_km.text(0.05, 0.12, f"p = {p_value:.3f}", transform=ax_km.transAxes, fontsize=13)
    ax_km.spines[["top", "right"]].set_visible(False)

    xmax = df["survival_months"].max()
    time_points = np.arange(0, xmax + 20, 20)
    ax_km.set_xlim(-5, xmax + 10)
    ax_km.set_xticks(time_points)

    fig.suptitle(f"{cohort_name} LUAD | {target_gene} median split", fontsize=14, fontweight="bold", y=0.98)
    fig.text(0.5, 0.925, target_gene, ha="center", fontsize=10, fontweight="bold")
    fig.text(0.44, 0.905, "— High", ha="center", fontsize=10, color=COLOR_HIGH)
    fig.text(0.58, 0.905, "— Low", ha="center", fontsize=10, color=COLOR_LOW)

    def n_at_risk(fitted_data, tps):
        return [(fitted_data["survival_months"] >= t).sum() for t in tps]

    risk_high = n_at_risk(high, time_points)
    risk_low = n_at_risk(low, time_points)

    ax_risk.set_xlim(ax_km.get_xlim())
    ax_risk.set_ylim(0, 1)
    ax_risk.axis("off")
    ax_risk.text(-0.20, 0.95, "Number at risk", transform=ax_risk.transAxes,
                 fontsize=10, fontweight="bold", ha="left", va="top")
    ax_risk.text(-0.02, 0.62, "High", transform=ax_risk.transAxes,
                 fontsize=10, color=COLOR_HIGH, ha="right", va="center")
    ax_risk.text(-0.02, 0.22, "Low", transform=ax_risk.transAxes,
                 fontsize=10, color=COLOR_LOW, ha="right", va="center")

    for t, r_high, r_low in zip(time_points, risk_high, risk_low):
        ax_risk.text(t, 0.62, str(r_high), fontsize=9, ha="center", va="center",
                     color=COLOR_HIGH, transform=ax_risk.transData)
        ax_risk.text(t, 0.22, str(r_low), fontsize=9, ha="center", va="center",
                     color=COLOR_LOW, transform=ax_risk.transData)

    fig.text(0.55, 0.045, "Time (months)", ha="center", fontsize=11)
    plt.subplots_adjust(left=0.20, right=0.95, top=0.88, bottom=0.10)

    out_path = out_dir / f"{cohort_name}_{target_gene}_KMplot.png"
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"[{cohort_name}] High n={len(high)}, Low n={len(low)}, log-rank p={p_value:.4f}")
    print(f"저장 완료: {out_path}")


def main():
    import importlib

    gene_config = importlib.import_module(GENE_CONFIG_MODULE)
    target_gene = gene_config.TARGET_GENE

    out_dir = pcfg.results(target_gene, "km_plots")

    for cohort_name in COHORT_RUNS:
        path = pcfg.processed(cohort_name) / f"{cohort_name}_merged_valid.csv"
        make_km_plot(cohort_name, path, target_gene, out_dir)


if __name__ == "__main__":
    main()
