# -*- coding: utf-8 -*-
"""
analysis_specific/correlation_table_specific.py
==================================================
TARGET_GENE과 SIGNATURE_GENES 각각의 Pearson 상관분석 표 (R, p-value, N).
THEME_GROUPS(gene_config)에 정의된 테마별로 표를 나눠서도 출력한다.

기존 파일과의 대응: correlation_analysis.py의 3번 섹션(GSE36471용)
-> 코호트/유전자 파라미터화해서 일반화.
"""

import sys
from pathlib import Path

import pandas as pd
from scipy import stats

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
# ================================


def run_correlation(cohort_name: str, merged_path: str, target_gene: str, signature_genes: list, out_dir: Path):
    merged = pd.read_csv(merged_path)

    rows = []
    for gene in signature_genes:
        valid = merged[[target_gene, gene]].dropna()
        r, p = stats.pearsonr(valid[target_gene], valid[gene])
        rows.append({"Gene": gene, "R": round(r, 3), "P_value": round(p, 4), "N": len(valid)})

    corr_df = pd.DataFrame(rows)
    out_path = out_dir / f"{cohort_name}_{target_gene}_correlation_results.csv"
    corr_df.to_csv(out_path, index=False)
    print(f"=== {cohort_name} | {target_gene} 상관분석 ===")
    print(corr_df)
    print(f"저장 완료: {out_path}\n")
    return corr_df


def main():
    import importlib

    gene_config = importlib.import_module(GENE_CONFIG_MODULE)
    target_gene = gene_config.TARGET_GENE
    signature_genes = gene_config.SIGNATURE_GENES

    out_dir = pcfg.results(target_gene, "correlation")

    for cohort_name in COHORT_RUNS:
        path = pcfg.processed(cohort_name) / f"{cohort_name}_merged_final.csv"
        corr_df = run_correlation(cohort_name, path, target_gene, signature_genes, out_dir)

        for theme, genes in getattr(gene_config, "THEME_GROUPS", {}).items():
            sub = corr_df[corr_df["Gene"].isin(genes)]
            print(f"--- {cohort_name} | 테마: {theme} ---")
            print(sub)


if __name__ == "__main__":
    main()
