# -*- coding: utf-8 -*-
"""
cohort_runners/GSE68465_runner_com.py
========================================
GSE68465 (LUAD, GPL96 = HG-U133A) 코호트 실행 스크립트.

주의: 이 코호트는 예전에 4번째 검증 코호트 후보로 시도했다가 보류된 상태입니다.
GPL96은 2000년대 초반에 설계된 구형 칩이라 CD300LF, HLA-DRA 등 패널 유전자
일부가 아예 probe로 존재하지 않을 수 있습니다 (map_GSE68465_genes.py에서
이미 이 문제를 확인함). forest plot에는 포함되지 않았던 코호트이니,
다시 시도할 경우 아래 실행 후 콘솔에 뜨는 "MISSING" 유전자 목록부터 확인하세요.

기존 파일과의 대응:
  map_GSE68465_genes.py, parse_GSE68465_clinical.py, parse_GSE68465_metadata.py
  -> 이 스크립트로 대체 (단, extract/merge 단계까지 실제로 실행된 적은 없었음).

실행 전 필요한 원본 파일 (raw_data/GSE68465/ 에 둘 것 — paths_config.py 참고):
  GSE68465_series_matrix.txt, GPL96.txt
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import pipeline_common as pc
import paths_config as pcfg

# ===== 여기만 바꾸면 됩니다 =====
COHORT_NAME = "GSE68465"
GENE_CONFIG_MODULE = "gene_config.CD300LF_panel"

DATA_DIR = pcfg.raw(COHORT_NAME)
OUT_DIR = pcfg.processed(COHORT_NAME)

SERIES_MATRIX_FILE = DATA_DIR / "GSE68465_series_matrix.txt"
PLATFORM_FILE = DATA_DIR / "GPL96.txt"
SYMBOL_COL = "Gene Symbol"
PROBE_ID_COL = "ID"

CLINICAL_FIELD_MAP = {
    "disease_state": "disease_state",
    "gender": "sex",
    "age": "age",
    "stage": "disease_stage",
    "smoking_status": "smoking_history",
}
VITAL_STATUS_FIELD = "vital_status"
SURVIVAL_MONTHS_KEYWORD = "months_to_last_contact_or_death"   # 이미 개월 단위
EXCLUDE_FIELD_KEYWORD = None
# ================================


def main():
    import importlib

    gene_config = importlib.import_module(GENE_CONFIG_MODULE)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"=== [{COHORT_NAME}] 1. 임상 데이터 파싱 ===")
    meta_df, _ = pc.parse_series_matrix_metadata(str(SERIES_MATRIX_FILE))
    raw_char = pc.extract_raw_characteristics(meta_df)
    clinical = pc.build_clinical(
        raw_char,
        field_keyword_map=CLINICAL_FIELD_MAP,
        vital_status_field=VITAL_STATUS_FIELD,
        survival_months_keyword=SURVIVAL_MONTHS_KEYWORD,
        exclude_field_keyword=EXCLUDE_FIELD_KEYWORD,
    )
    clinical_path = OUT_DIR / f"{COHORT_NAME}_clinical.csv"
    clinical.to_csv(clinical_path, index=False)
    print(f"저장 완료: {clinical_path}")

    print(f"\n=== [{COHORT_NAME}] 2. probe -> 유전자 매핑 ({gene_config.TARGET_GENE} 패널) ===")
    print("!! GPL96은 구형 칩이라 일부 유전자가 MISSING으로 뜰 수 있습니다 !!")
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
