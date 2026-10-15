# -*- coding: utf-8 -*-
"""
tcga/TCGA_runner_com.py
==========================
TCGA(LUAD 또는 LUSC) 발현매트릭스 + 임상 파일 -> GEO 코호트와 동일한 형태의
merged_final.csv / merged_valid.csv 를 만든다.

핵심: 여기서 만든 csv는 analysis_specific/*.py가 기대하는 형태
(TARGET_GENE + SIGNATURE_GENES 컬럼, survival_months, survival_status)와
완전히 같습니다. 그래서 이 스크립트를 실행한 뒤에는 새 분석 스크립트를 짤
필요 없이, analysis_specific/cox_regression_specific.py 등의
COHORT_RUNS 리스트에 아래처럼 한 줄만 추가하면 TCGA도 GSE31210/GSE72094/
GSE36471과 나란히 forest plot에 들어갑니다:

    COHORT_RUNS = [
        ("GSE31210", "../data_processed/GSE31210/GSE31210_merged_valid.csv"),
        ("GSE72094", "../data_processed/GSE72094/GSE72094_merged_valid.csv"),
        ("GSE36471", "../data_processed/GSE36471/GSE36471_merged_valid.csv"),
        ("TCGA_LUAD", "../data_processed/TCGA_LUAD/TCGA_LUAD_merged_valid.csv"),  # <- 추가
    ]

주의 (실제 파일에 맞게 꼭 확인/수정하세요):
  - CLINICAL_* 컬럼명은 임상 파일 출처(GDC/cBioPortal/TCGA-CDR)마다 다릅니다.
    아래는 GDC clinical.tsv 스타일의 흔한 컬럼명을 기본값으로 넣어뒀을 뿐이니,
    실제 파일을 열어서 컬럼명을 맞춰주세요. (스크립트 실행 시 KeyError가 나면
    그게 컬럼명이 다르다는 뜻입니다 - 에러 메시지에 실제 컬럼 목록이 함께 출력됩니다)
  - EXPR_MATRIX_FILE은 build_gsea_cls_specific.py에서 쓰는 것과 같은,
    이미 유전자 심볼로 인덱싱된 매트릭스라고 가정합니다 (probe 매핑 불필요).
  - TUMOR_TYPE_CODE="01"(원발암)로 기본 필터링합니다. 매트릭스에 이미
    종양 샘플만 있다면 None으로 바꾸세요.
"""

import sys
from pathlib import Path

# tcga/ 폴더 자신 (tcga_common.py가 여기 있음 - 스크립트를 직접 실행하면 파이썬이
# 자동으로 넣어주지만, 다른 위치에서 import될 때도 되도록 명시적으로 추가)
sys.path.insert(0, str(Path(__file__).resolve().parent))
# pipeline_v2 루트 (gene_config 패키지를 찾기 위해 필요)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import tcga_common as tc  # noqa: E402
import paths_config as pcfg  # noqa: E402

# ===== 여기만 바꾸면 됩니다 =====
CANCER_TYPE = "LUAD"          # "LUAD" 또는 "LUSC"
GENE_CONFIG_MODULE = "gene_config.CD300LF_panel"

# raw_data/TCGA_LUAD/ 하나에 GSEA용 매트릭스/GCT와 임상 파일을 같이 둡니다
# (build_gsea_cls_specific.py도 같은 폴더를 봅니다 - 매트릭스 파일을 두 번 둘 필요 없음)
DATA_DIR = pcfg.raw(f"TCGA_{CANCER_TYPE}")
EXPR_MATRIX_FILE = DATA_DIR / f"{CANCER_TYPE}_expr_matrix_clean.txt"
CLINICAL_FILE = DATA_DIR / f"{CANCER_TYPE}_clinical.tsv"   # 실제 파일명/확장자에 맞게 수정

OUT_DIR = pcfg.processed(f"TCGA_{CANCER_TYPE}")

# --- 임상 파일 컬럼명 (실제 파일 열어서 맞는지 꼭 확인) ---
CLINICAL_SEP = "\t"
PATIENT_ID_COL = "case_submitter_id"      # 예: "TCGA-05-4244" 형태
VITAL_STATUS_COL = "vital_status"
DEAD_VALUES = ["Dead", "dead", "1", "Deceased", "DECEASED"]
TIME_COL = "days_to_death_or_last_followup"  # 없으면 days_to_death/days_to_last_follow_up 중
                                              # 하나로 직접 만들어서 넣어주세요
TIME_UNIT = "days"
EXTRA_FIELD_MAP = {
    "age": "age_at_index",
    "gender": "gender",
    "stage": "ajcc_pathologic_stage",
}

TUMOR_TYPE_CODE = "01"   # TCGA barcode의 sample-type 코드. None이면 필터 안 함.
# ================================


def main():
    import importlib

    gene_config = importlib.import_module(GENE_CONFIG_MODULE)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cohort_name = f"TCGA_{CANCER_TYPE}"

    print(f"=== [{cohort_name}] 1. 임상 데이터 파싱 ===")
    clinical = tc.load_tcga_clinical(
        str(CLINICAL_FILE),
        patient_id_col=PATIENT_ID_COL,
        vital_status_col=VITAL_STATUS_COL,
        dead_values=DEAD_VALUES,
        time_col=TIME_COL,
        time_unit=TIME_UNIT,
        extra_field_map=EXTRA_FIELD_MAP,
        sep=CLINICAL_SEP,
    )
    clinical.to_csv(OUT_DIR / f"{cohort_name}_clinical.csv", index=False)

    print(f"\n=== [{cohort_name}] 2. {gene_config.TARGET_GENE} 패널 발현 추출 ===")
    panel_expr = tc.extract_panel_expression(str(EXPR_MATRIX_FILE), gene_config.FULL_PANEL)
    panel_expr.to_csv(OUT_DIR / f"{cohort_name}_{gene_config.TARGET_GENE}_panel_expression.csv")

    print(f"\n=== [{cohort_name}] 3. barcode 매칭 + 병합 ===")
    merged, merged_valid = tc.match_and_merge(clinical, panel_expr, tumor_type_code=TUMOR_TYPE_CODE)
    merged.to_csv(OUT_DIR / f"{cohort_name}_merged_final.csv", index=False)
    merged_valid.to_csv(OUT_DIR / f"{cohort_name}_merged_valid.csv", index=False)
    print(f"저장 완료: {OUT_DIR}/{cohort_name}_merged_final.csv, {cohort_name}_merged_valid.csv")
    print("\n이제 analysis_specific/*.py의 COHORT_RUNS에 아래 줄을 추가하면 됩니다:")
    print(f'    ("{cohort_name}", "../data_processed/{cohort_name}/{cohort_name}_merged_valid.csv"),')


if __name__ == "__main__":
    main()
