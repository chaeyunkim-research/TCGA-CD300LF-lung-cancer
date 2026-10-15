# -*- coding: utf-8 -*-
"""
analysis_specific/coexpression_clustermap_specific.py
========================================================
TARGET_GENE + SIGNATURE_GENES 패널의 Spearman 상관계수 clustermap.

기존 파일과의 대응: clustermap_GSE31210.py, clustermap_GSE72094.py,
make_clustermap.py (GSE36471용) -> 코호트명/유전자만 다르고 로직은 100% 동일했으므로 통합.
"""

import sys
from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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


def make_clustermap(cohort_name: str, expr_path: str, target_gene: str, out_dir: Path):
    expr = pd.read_csv(expr_path, index_col=0)
    corr_matrix = expr.T.corr(method="spearman")

    g = sns.clustermap(
        corr_matrix, annot=True, fmt=".2f",
        cmap="RdBu_r", vmin=-1, vmax=1, center=0,
        figsize=(11, 10), linewidths=0.5, linecolor="white",
        dendrogram_ratio=0.15, cbar_pos=(1.02, 0.3, 0.03, 0.4),
        cbar_kws={"label": "Spearman correlation"}, method="average",
    )
    g.ax_heatmap.set_xlabel("GENE", fontsize=11)
    g.ax_heatmap.set_ylabel("GENE", fontsize=11)
    g.fig.suptitle(f"{cohort_name} LUAD | Co-expression with {target_gene} (clustered)",
                   fontsize=13, fontweight="bold", y=1.02)
    for text in g.ax_heatmap.texts:
        text.set_fontsize(8)

    out_path = out_dir / f"{cohort_name}_{target_gene}_coexpression_clustermap.png"
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(g.fig)
    print(f"저장 완료: {out_path}")


def main():
    import importlib

    gene_config = importlib.import_module(GENE_CONFIG_MODULE)
    target_gene = gene_config.TARGET_GENE

    out_dir = pcfg.results(target_gene, "clustermaps")

    for cohort_name in COHORT_RUNS:
        path = pcfg.processed(cohort_name) / f"{cohort_name}_{target_gene}_panel_expression.csv"
        make_clustermap(cohort_name, path, target_gene, out_dir)


if __name__ == "__main__":
    main()
