# -*- coding: utf-8 -*-
"""
cohort_runners/GSE36471_runner_com.py
========================================
GSE36471 (LUAD, GPL9053 커스텀 어레이) 코호트 실행 스크립트.
가장 먼저 만들었던 "프로토타입" 코호트라 원본 파일들이
GSE36471이라는 이름 없이 범용 이름(parse_geo_metadata.py, extract_expression.py,
correlation_analysis.py, cox_regression.py, km_plot.py, km_plot_styled(_v2).py,
make_heatmap.py, make_clustermap.py, find_target_genes*.py, find_aliases.py,
map_final15.py, check_cd300lf.py)으로 흩어져 있던 것이 이번 정리의 계기였던
그 코호트입니다. 이제 다른 코호트와 동일한 이름 규칙으로 통일합니다.

주의(다른 코호트와의 차이점): GSE36471의 series_matrix characteristics 필드에는
"survival_status"와 "survival_months"가 이미 0/1과 개월 단위 숫자로 직접 들어있어서
(다른 코호트처럼 vital_status 텍스트 -> 0/1 변환이 필요 없음), 아래에서
survival_status_keyword를 쓰고 vital_status_field는 쓰지 않는다.

실행 전 필요한 원본 파일 (raw_data/GSE36471/ 에 둘 것 — paths_config.py 참고):
  GSE36471-GPL9053_series_matrix.txt, GPL9053.txt
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import pipeline_common as pc
import paths_config as pcfg

# ===== 여기만 바꾸면 됩니다 =====
COHORT_NAME = "GSE36471"
GENE_CONFIG_MODULE = "gene_config.CD300LF_panel"

DATA_DIR = pcfg.raw(COHORT_NAME)
OUT_DIR = pcfg.processed(COHORT_NAME)

SERIES_MATRIX_FILE = DATA_DIR / "GSE36471-GPL9053_series_matrix.txt"
PLATFORM_FILE = DATA_DIR / "GPL9053.txt"
SYMBOL_COL = "ORF"          # GPL9053 커스텀 어레이는 gene symbol 컬럼명이 "ORF"
PROBE_ID_COL = "ID"

CLINICAL_FIELD_MAP = {
    "subtype": "subtype",
    "age": "age",
    "gender": "sex",
    "smoking_status": "smoking_status",
    "stage": "stage",
}
VITAL_STATUS_FIELD = None
SURVIVAL_STATUS_KEYWORD = "survival_status"   # 다른 코호트와 달리 이미 0/1로 제공됨
SURVIVAL_MONTHS_KEYWORD = "survival_months"   # 이미 개월 단위로 제공됨
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
        survival_status_keyword=SURVIVAL_STATUS_KEYWORD,
        survival_months_keyword=SURVIVAL_MONTHS_KEYWORD,
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
