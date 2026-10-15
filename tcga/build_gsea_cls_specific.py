# -*- coding: utf-8 -*-
"""
tcga/build_gsea_cls_specific.py
=======================================
TCGA 발현매트릭스에서 TARGET_GENE의 발현값을 GSEA continuous phenotype(.cls) 파일로
만든다. GSEA(GenePattern/desktop)에서 이 .cls + 기존 GCT 파일을 continuous phenotype
correlation 모드로 돌리면, TARGET_GENE 발현과 상관된 pathway/gene set을 찾을 수 있다.

기존 파일과의 대응: rerun_PTPN6_tcga.py -> TARGET_GENE/암종을 파라미터화해서 일반화.
(원래 이름이 "PTPN6"에 고정돼 있던 게 이번 혼동의 원인 중 하나였습니다.
 이제 CANCER_TYPE + TARGET_GENE만 바꾸면 LUAD/LUSC 어느 쪽이든, 어떤 유전자든 재사용 가능.)

주의: GCT 파일은 이미 만들어져 있는 것을 그대로 재사용합니다 (여기서 새로 만들지 않음).
      만든 .cls 파일의 샘플 순서가 그 GCT의 샘플 순서와 반드시 일치해야 하므로,
      만든 뒤 check_gct_cls_match_com.py로 샘플 수가 맞는지 검증하세요.

실행 전 필요한 원본 파일 (raw_data/TCGA_LUAD/ 에 둘 것 — paths_config.py 참고.
TCGA_runner_com.py와 같은 폴더를 보므로 이 매트릭스 파일을 두 번 둘 필요 없음):
  LUAD_expr_matrix_clean.txt (또는 LUSC)
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths_config as pcfg  # noqa: E402

# ===== 여기만 바꾸면 됩니다 =====
CANCER_TYPE = "LUAD"          # "LUAD" 또는 "LUSC"
TARGET_GENE = "CD300LF"       # 예: "CD300LF", "PTPN6" 등 GSEA로 보고 싶은 아무 유전자

DATA_DIR = pcfg.raw(f"TCGA_{CANCER_TYPE}")
EXPR_MATRIX_FILE = DATA_DIR / f"{CANCER_TYPE}_expr_matrix_clean.txt"
OUTPUT_PREFIX = pcfg.processed(f"TCGA_{CANCER_TYPE}") / f"{CANCER_TYPE}_{TARGET_GENE}"
# ================================


def main():
    expr = pd.read_csv(EXPR_MATRIX_FILE, sep="\t", index_col=0)

    if TARGET_GENE not in expr.index:
        candidates = [g for g in expr.index if TARGET_GENE.upper() in str(g).upper()][:10]
        print(f"'{TARGET_GENE}' 매트릭스에 없음. 유사 후보: {candidates}")
        return

    target_expr = expr.loc[TARGET_GENE]
    cls_path = f"{OUTPUT_PREFIX}.cls"
    with open(cls_path, "w") as f:
        f.write("#numeric\n")
        f.write(f"#{TARGET_GENE}\n")
        f.write(" ".join(target_expr.astype(str).tolist()) + "\n")

    print(f"저장 완료: {cls_path} ({len(target_expr)} samples)")
    print(f"주의: 기존 GCT 파일({CANCER_TYPE}_matrix.gct)의 샘플 순서와 반드시 일치해야 합니다.")
    print("-> tcga/check_gct_cls_match_com.py 로 샘플 수 일치 여부를 확인하세요.")


if __name__ == "__main__":
    main()
