# -*- coding: utf-8 -*-
"""
cohort_runners/GSE31210_runner_com.py   ("_com" = 코호트 파이프라인, 공통 엔진 사용)
==========================================================================
GSE31210 (LUAD, GPL570) 코호트를 raw GEO 파일 -> 표준화된 clinical/expression/merged
csv까지 만드는 실행 스크립트. 로직 자체는 전부 pipeline_common.py에 있고,
여기는 "이 코호트만의 값"만 채운다.

기존 파일과의 대응:
  extract_GSE31210_expression.py, map_GSE31210_genes.py(v1~v4),
  parse_GSE31210_clinical.py(v1~v3), parse_GSE31210_metadata.py, merge_GSE31210.py
  -> 이 스크립트 하나로 대체.

주의: parse_GSE31210_clinical 은 v1 -> v2 -> v3로 가면서 로직이 바뀌었고,
실제 최종 결과(GSE31210_merged_valid.csv)는 v3 로직(= days -> months 변환,
"age (years)"/"pathological stage"/"days before death/censor" 키워드,
"exclude for prognosis" 필터)으로 만들어졌음을 확인함. 아래 설정은 v3을 따른다.

실행 전 필요한 원본 파일 (이 스크립트 폴더가 아니라 raw_data/GSE31210/ 에 둘 것 —
paths_config.py 참고. 어느 위치에서 이 스크립트를 실행하든 항상 그 폴더를 봅니다):
  GSE31210_series_matrix.txt, GPL570.txt
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import pipeline_common as pc
import paths_config as pcfg

# ===== 여기만 바꾸면 됩니다 =====
COHORT_NAME = "GSE31210"
GENE_CONFIG_MODULE = "gene_config.CD300LF_panel"   # 다른 유전자로 바꾸려면 이 문자열만 교체

DATA_DIR = pcfg.raw(COHORT_NAME)      # raw_data/GSE31210/ (실행 위치와 무관하게 항상 이 폴더)
OUT_DIR = pcfg.processed(COHORT_NAME)  # data_processed/GSE31210/

SERIES_MATRIX_FILE = DATA_DIR / "GSE31210_series_matrix.txt"
PLATFORM_FILE = DATA_DIR / "GPL570.txt"
SYMBOL_COL = None          # None = 자동 감지 (GPL570은 "Gene Symbol")
PROBE_ID_COL = "ID"

CLINICAL_FIELD_MAP = {
    "age": "age (years)",
    "gender": "gender",
    "smoking_status": "smoking status",
    "stage": "pathological stage",
    "relapse": "relapse",
}
VITAL_STATUS_FIELD = "death"
SURVIVAL_DAYS_KEYWORD = "days before death/censor"
EXCLUDE_FIELD_KEYWORD = "exclude for prognosis"
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
