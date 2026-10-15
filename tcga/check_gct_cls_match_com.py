# -*- coding: utf-8 -*-
"""
tcga/check_gct_cls_match_com.py
=======================================
GCT(발현 매트릭스)와 CLS(phenotype) 파일의 샘플 수가 일치하는지 확인하는
공통 점검 스크립트. GSEA는 두 파일의 샘플 순서가 어긋나면 조용히 잘못된
결과를 내므로, build_gsea_cls_specific.py로 새 .cls를 만들 때마다 항상 이걸로 확인하세요.

기존 파일과의 대응: check_gct_cls_match.py
  (기존 버전은 "LUAD_matrix.gct"/"LUAD_PTPN6.cls"가 파일명에 그대로 박혀 있어서
   유전자를 CD300LF로 바꾼 뒤에도 실수로 PTPN6 결과를 검증하게 될 위험이 있었습니다.
   아래는 파일 경로를 인자로 받는 공통 버전입니다.)
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths_config as pcfg  # noqa: E402

# ===== 여기만 바꾸면 됩니다 =====
CANCER_TYPE = "LUAD"
TARGET_GENE = "CD300LF"

GCT_FILE = pcfg.raw(f"TCGA_{CANCER_TYPE}") / f"{CANCER_TYPE}_matrix.gct"
CLS_FILE = pcfg.processed(f"TCGA_{CANCER_TYPE}") / f"{CANCER_TYPE}_{TARGET_GENE}.cls"
# ================================


def check_match(gct_file: Path, cls_file: Path):
    with open(gct_file, "r") as f:
        f.readline()  # #1.2
        dims = f.readline().strip().split("\t")
    n_genes, n_samples_gct = int(dims[0]), int(dims[1])
    print(f"GCT ({gct_file.name}): 유전자 {n_genes}개, 샘플 {n_samples_gct}개")

    with open(cls_file, "r") as f:
        f.readline()  # #numeric
        gene_line = f.readline().strip()  # #유전자명
        values_line = f.readline().strip()
    n_samples_cls = len(values_line.split())
    gene_name = gene_line.lstrip("#")
    print(f"CLS ({cls_file.name}, phenotype={gene_name}): 샘플 {n_samples_cls}개")

    if n_samples_gct == n_samples_cls:
        print(f"OK: 샘플 수 일치 ({n_samples_gct}개) — GSEA에 바로 사용 가능합니다.")
    else:
        print(f"불일치! GCT={n_samples_gct}개 vs CLS={n_samples_cls}개 — "
              f"같은 발현매트릭스에서 만든 .cls가 맞는지, 샘플 순서가 같은지 다시 확인하세요.")


if __name__ == "__main__":
    check_match(GCT_FILE, CLS_FILE)
