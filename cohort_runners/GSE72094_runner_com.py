# -*- coding: utf-8 -*-
"""
cohort_runners/GSE72094_runner_com.py
========================================
GSE72094 (LUAD, GPL15048) 코호트 실행 스크립트.

기존 파일과의 대응:
  extract_GSE72094_expression.py, map_GSE72094_genes.py, parse_GSE72094_clinical.py,
  parse_GSE72094_metadata.py, merge_GSE72094.py, check_gpl15048_structure.py
  -> 이 스크립트 하나로 대체.

실행 전 필요한 원본 파일 (raw_data/GSE72094/ 에 둘 것 — paths_config.py 참고):
  GSE72094_series_matrix.txt, GPL15048.txt
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import pipeline_common as pc
import paths_config as pcfg

# ===== 여기만 바꾸면 됩니다 =====
COHORT_NAME = "GSE72094"
GENE_CONFIG_MODULE = "gene_config.CD300LF_panel"

DATA_DIR = pcfg.raw(COHORT_NAME)
OUT_DIR = pcfg.processed(COHORT_NAME)

SERIES_MATRIX_FILE = DATA_DIR / "GSE72094_series_matrix.txt"
PLATFORM_FILE = DATA_DIR / "GPL15048.txt"
SYMBOL_COL = "GeneSymbol"   # GPL15048은 컬럼명에 공백이 없음 (map_GSE72094_genes.py 확인)
PROBE_ID_COL = "ID"

# 원본 텍스트가 "age_at_diagnosis: 65"처럼 언더스코어로 들어있는 코호트
CLINICAL_FIELD_MAP = {
    "patient_id": "patient_id",
    "gender": "gender",
    "age": "age_at_diagnosis",
    "smoking_status": "smoking_status",
    "stage": "stage",
    "kras_status": "kras_status",
    "egfr_status": "egfr_status",
    "stk11_status": "stk11_status",
    "tp53_status": "tp53_status",
}
VITAL_STATUS_FIELD = "vital_status"
SURVIVAL_DAYS_KEYWORD = "survival_time_in_days"
EXCLUDE_FIELD_KEYWORD = None
# ================================


def main():
    import importlib

    gene_config = importlib.import_module(GENE_CONFIG_MODULE)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"=== [{COHORT_NAME}] 0. 원본 파일 확인 ===")
    missing = [p for p in (SERIES_MATRIX_FILE, PLATFORM_FILE) if not p.exists()]
    if missing:
        print(f"[{COHORT_NAME}] 다음 파일을 찾을 수 없어서 중단합니다:")
        for m in missing:
            print(f"  - {m}")
        if DATA_DIR.exists():
            found = [p.name for p in DATA_DIR.iterdir()]
            print(f"[{COHORT_NAME}] 현재 {DATA_DIR} 안에 있는 파일: {found if found else '(비어있음)'}")
        else:
            print(f"[{COHORT_NAME}] {DATA_DIR} 폴더 자체가 없습니다.")
        print(f"[{COHORT_NAME}] 위 파일을 raw_data/{COHORT_NAME}/ 안에 넣고 다시 실행하세요.")
        return
    print(f"사용할 파일: {SERIES_MATRIX_FILE.name}, {PLATFORM_FILE.name}")

    print(f"\n=== [{COHORT_NAME}] 1. 임상 데이터 파싱 ===")
    meta_df, _ = pc.parse_series_matrix_metadata(str(SERIES_MATRIX_FILE))
    raw_char = pc.extract_raw_characteristics(meta_df)
    clinical = pc.build_clinical(
        raw_char,
        field_keyword_map=CLINICAL_FIELD_MAP,
        vital_status_field=VITAL_STATUS_FIELD,
        survival_days_keyword=SURVIVAL_DAYS_KEYWORD,
        exclude_field_keyword=EXCLUDE_FIELD_KEYWORD,
    )
    clinical_path = OUT_DIR / f"{COHORT_NAME}_clinical.csv"
    clinical.to_csv(clinical_path, index=False)
    print(f"저장 완료: {clinical_path}")

    print(f"\n=== [{COHORT_NAME}] 2. probe -> 유전자 매핑 ({gene_config.TARGET_GENE} 패널) ===")
    probe_map = pc.map_probes_to_genes(
        str(PLATFORM_FILE), gene_config.ALIAS_MAP, symbol_col=SYMBOL_COL, id_col=PROBE_ID_COL
    )
    probe_map_path = OUT_DIR / f"{COHORT_NAME}_probe_mapping.csv"
    probe_map.to_csv(probe_map_path, index=False)
    print(f"저장 완료: {probe_map_path}")

    print(f"\n=== [{COHORT_NAME}] 3. 유전자 발현 매트릭스 추출 ===")
    gene_expr = pc.extract_gene_expression(str(SERIES_MATRIX_FILE), probe_map, probe_id_col=PROBE_ID_COL)
    expr_path = OUT_DIR / f"{COHORT_NAME}_{gene_config.TARGET_GENE}_panel_expression.csv"
    gene_expr.to_csv(expr_path)
    print(f"저장 완료: {expr_path}")

    print(f"\n=== [{COHORT_NAME}] 4. 임상 + 발현 병합 ===")
    merged, merged_valid = pc.merge_clinical_expression(clinical, gene_expr)
    merged.to_csv(OUT_DIR / f"{COHORT_NAME}_merged_final.csv", index=False)
    merged_valid.to_csv(OUT_DIR / f"{COHORT_NAME}_merged_valid.csv", index=False)
    print(f"저장 완료: {OUT_DIR}/{COHORT_NAME}_merged_final.csv, {COHORT_NAME}_merged_valid.csv")


if __name__ == "__main__":
    main()
