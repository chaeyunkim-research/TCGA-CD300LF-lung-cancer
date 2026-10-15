# -*- coding: utf-8 -*-
"""
analysis_specific/cox_regression_specific.py  ("_specific" = 타겟 유전자에 종속된 분석)
========================================================================
merged_valid.csv 하나에 대해 TARGET_GENE(z-score 표준화) 기준 단변량 Cox 회귀를 돌리고,
forest_plot_specific.py에서 바로 합칠 수 있는 표준 포맷(csv)으로 저장한다.

기존 파일과의 대응: cox_GSE31210.py, cox_GSE72094.py, cox_regression.py (GSE36471용)
-> 셋이 사실상 동일한 로직이었으므로 COHORT_NAME + 입력파일만 바꾸는 이 스크립트 하나로 통합.

여러 코호트를 처리하려면 COHORT_RUNS 리스트에 코호트명만 추가하면 됩니다
(파일 경로는 paths_config를 통해 자동으로 찾습니다 - 실행 위치와 무관하게 항상
같은 data_processed/<코호트>/ 를 봅니다). 결과는 results/<GENE>/cox/ 밑에
코호트별 csv + 하나로 합친 "<GENE>_cox_all_cohorts.csv"로 저장됩니다
(forest_plot_specific.py의 입력).
"""

import sys
from pathlib import Path

import pandas as pd
from lifelines import CoxPHFitter

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths_config as pcfg  # noqa: E402

# ===== 여기만 바꾸면 됩니다 =====
GENE_CONFIG_MODULE = "gene_config.CD300LF_panel"

# 코호트명 목록. 코호트를 추가/제거하려면 여기만 수정.
COHORT_RUNS = [
    "GSE31210",
    "GSE72094",
    "GSE36471",
    "GSE13213",
    "GSE42127",
    # TCGA를 추가하려면 tcga/TCGA_runner_com.py 실행 후 아래 주석 해제:
    # "TCGA_LUAD",
]
# ================================


def run_cox_for_cohort(cohort_name: str, merged_valid_path: str, target_gene: str) -> dict:
    merged = pd.read_csv(merged_valid_path)
    df = merged[[target_gene, "survival_months", "survival_status"]].dropna().copy()

    gene_z_col = f"{target_gene}_z"
    df[gene_z_col] = (df[target_gene] - df[target_gene].mean()) / df[target_gene].std()

    cph = CoxPHFitter()
    cph.fit(
        df[[gene_z_col, "survival_months", "survival_status"]],
        duration_col="survival_months",
        event_col="survival_status",
    )

    hr = cph.summary.loc[gene_z_col, "exp(coef)"]
    ci_low = cph.summary.loc[gene_z_col, "exp(coef) lower 95%"]
    ci_high = cph.summary.loc[gene_z_col, "exp(coef) upper 95%"]
    p_val = cph.summary.loc[gene_z_col, "p"]

    print(f"=== {cohort_name} | {target_gene} Cox regression 결과 ===")
    print(f"HR = {hr:.3f} (95% CI: {ci_low:.3f}-{ci_high:.3f}), p = {p_val:.4f}, N = {len(df)}")

    return {
        "Cohort": f"{cohort_name}_LUAD",
        "HR": round(hr, 3),
        "CI_lower": round(ci_low, 3),
        "CI_upper": round(ci_high, 3),
        "p_value": round(p_val, 4),
        "N": len(df),
    }


def main():
    import importlib

    gene_config = importlib.import_module(GENE_CONFIG_MODULE)
    target_gene = gene_config.TARGET_GENE

    out_dir = pcfg.results(target_gene, "cox")

    rows = []
    for cohort_name in COHORT_RUNS:
        path = pcfg.processed(cohort_name) / f"{cohort_name}_merged_valid.csv"
        row = run_cox_for_cohort(cohort_name, path, target_gene)
        rows.append(row)
        pd.DataFrame([row]).to_csv(out_dir / f"{cohort_name}_{target_gene}_cox_result.csv", index=False)

    all_df = pd.DataFrame(rows)
    all_path = out_dir / f"{target_gene}_cox_all_cohorts.csv"
    all_df.to_csv(all_path, index=False)
    print(f"\n저장 완료: {all_path} (forest_plot_specific.py 입력으로 사용)")


if __name__ == "__main__":
    main()
