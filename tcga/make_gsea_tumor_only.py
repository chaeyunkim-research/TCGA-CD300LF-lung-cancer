# -*- coding: utf-8 -*-
"""
Build tumor-only GSEA input files (.gct + continuous-phenotype .cls) for CD300LF
from TCGA Firehose RNASeqV2 RSEM normalized data (run 2016-01-28).

Keeps primary tumor samples only (TCGA barcode sample type "01"), one sample per
patient, and drops genes without a symbol ("?") and duplicated symbols.
Expression values are left on the original RSEM normalized scale, as in the
earlier GSEA run.

Usage (from the Downloads folder):
    python make_gsea_tumor_only.py
Outputs are written to ./CD300LF_GSEA_tumor_only/
"""

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
TARGET_GENE = "CD300LF"
RUN = "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data"
SOURCES = {
    "LUAD": HERE / f"gdac.broadinstitute.org_LUAD.Merge_{RUN}.Level_3.2016012800.0.0" / f"LUAD.{RUN}.data.txt",
    "LUSC": HERE / f"gdac.broadinstitute.org_LUSC.Merge_{RUN}.Level_3.2016012800.0.0" / f"LUSC.{RUN}.data.txt",
}
OUT_DIR = HERE / "CD300LF_GSEA_tumor_only"


def build(cancer: str, path: Path):
    if not path.exists():
        print(f"[{cancer}] file not found: {path}")
        return
    expr = pd.read_csv(path, sep="\t", index_col=0, skiprows=[1], low_memory=False)
    n_all = expr.shape[1]

    # primary tumors only, one sample per patient
    tumor_cols = sorted(c for c in expr.columns if c[13:15] == "01")
    seen, keep = set(), []
    for c in tumor_cols:
        if c[:12] not in seen:
            seen.add(c[:12])
            keep.append(c)
    expr = expr[keep]

    # gene symbols: "TP53|7157" -> "TP53"; drop unnamed and duplicated symbols
    expr.index = [i.split("|")[0] for i in expr.index]
    expr = expr[(expr.index != "?") & (expr.index != "")]
    expr = expr[~expr.index.duplicated(keep="first")]
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0)
    expr.columns = [c[:12] for c in expr.columns]

    if TARGET_GENE not in expr.index:
        print(f"[{cancer}] {TARGET_GENE} not found")
        return

    OUT_DIR.mkdir(exist_ok=True)
    gct = OUT_DIR / f"{cancer}_tumor_only.gct"
    with open(gct, "w") as f:
        f.write("#1.2\n")
        f.write(f"{expr.shape[0]}\t{expr.shape[1]}\n")
        f.write("NAME\tDESCRIPTION\t" + "\t".join(expr.columns) + "\n")
        for gene, row in expr.iterrows():
            f.write(gene + "\tna\t" + "\t".join(f"{v:.4f}" for v in row.values) + "\n")

    # GSEA treats each distinct value of a continuous phenotype as one class, so
    # tied values (e.g. two LUSC samples at exactly 133.3333) make it fail.
    # Break ties with a negligible offset (1e-6 per repeat); this does not
    # change the Pearson ranking in any meaningful way.
    values = expr.loc[TARGET_GENE].astype(float).tolist()
    counts, n_ties = {}, 0
    for i, v in enumerate(values):
        k = round(v, 4)
        if k in counts:
            counts[k] += 1
            values[i] = v + 1e-6 * counts[k]
            n_ties += 1
        else:
            counts[k] = 0
    if n_ties:
        print(f"[{cancer}] {n_ties} tied {TARGET_GENE} value(s) offset by 1e-6")

    cls = OUT_DIR / f"{cancer}_{TARGET_GENE}_tumor_only.cls"
    with open(cls, "w") as f:
        f.write("#numeric\n")
        f.write(f"#{TARGET_GENE}\n")
        f.write(" ".join(f"{v:.6f}" for v in values) + "\n")

    print(f"[{cancer}] {n_all} samples in source -> {expr.shape[1]} primary tumors, "
          f"{expr.shape[0]} genes")
    print(f"    {gct.name}, {cls.name}")


if __name__ == "__main__":
    for cancer, path in SOURCES.items():
        build(cancer, path)
    print(f"Done. Files are in {OUT_DIR}")
