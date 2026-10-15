# -*- coding: utf-8 -*-
"""
cohort_runners/GSE13213_runner_com.py   ("_com" = 코호트 파이프라인, 공통 엔진 사용)
==========================================================================
GSE13213 (Okayama et al., LUAD, GPL6480 Agilent 2-channel) 코호트를 raw GEO
파일 -> 표준화된 clinical/expression/merged csv까지 만드는 실행 스크립트.

GSE13213_header_preview.txt (사용자가 PowerShell로 직접 뽑아준 메타데이터)를 보고
아래 설정을 채웠습니다. 확인한 내용:
  - 플랫폼: GPL6480 (Agilent Whole Human Genome 4x44K, 2-channel, VALUE=Cy5/Cy3 비율)
  - Series_overall_design: "Expression profiles in 117 patients with lung
    adenocarcinomas." -> 전원 ADC이므로 히스톨로지 subset 불필요.
  - characteristics_ch1 키(순서): Annotation, Cohort, Age, Sex, Histology,
    Smoking (BI), TNM (Pathological), Stage (Pathological ), Status(Dead/Alive),
    Survival (days), Evidence of relapse, Site of relapse, EGFR/K-ras/p53 status.
  - vital status는 "Status: Dead/Alive" 텍스트, survival은 "Survival (days)"
    (일 단위) -> GSE31210/GSE72094와 동일하게 30.4375로 나눠 개월 변환.
  - "exclude for prognosis" 같은 제외 플래그 필드는 없음 (GSE31210과 다른 점).

주의: pipeline_common._smart_open()이 .gz를 자동으로 처리하도록 이번에 고쳐서,
series_matrix.txt.gz / GPL6480.txt.gz를 압축 해제하지 않고 그대로 raw_data/에
넣어도 됩니다 (다른 코호트도 마찬가지로 이제 .gz를 바로 써도 됩니다).

실행 전 필요한 원본 파일 (raw_data/GSE13213/ 에 둘 것 — paths_config.py 참고):
  GSE13213_series_matrix.txt(.gz), GPL6480.txt(.gz)
  -> GPL6480.txt(.gz)는 GEO의 GPL6480 accession 페이지에서 "Download full table..."로 받으면 됩니다.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import pipeline_common as pc
import paths_config as pcfg

# ===== 여기만 바꾸면 됩니다 =====
COHORT_NAME = "GSE13213"
GENE_CONFIG_MODULE = "gene_config.CD300LF_panel"

DATA_DIR = pcfg.raw(COHORT_NAME)       # raw_data/GSE13213/
OUT_DIR = pcfg.processed(COHORT_NAME)  # data_processed/GSE13213/

# .gz를 압축 풀지 않고 그대로 둔 경우와 이미 풀어둔 경우 둘 다 자동으로 찾는다.
def _find_input(basename_no_ext: str) -> Path:
    for suffix in (".txt", ".txt.gz"):
        p = DATA_DIR / f"{basename_no_ext}{suffix}"
        if p.exists():
            return p
    # 못 찾아도 에러 대신 .txt 경로를 돌려줘서 아래 실행 시 안내 메시지가 뜨게 한다.
    return DATA_DIR / f"{basename_no_ext}.txt"


SERIES_MATRIX_FILE = _find_input("GSE13213_series_matrix")
PLATFORM_FILE = _find_input("GPL6480")
SYMBOL_COL = None          # None = 자동 감지 (GPL6480 Agilent annotation은 보통 "GeneSymbol")
PROBE_ID_COL = "ID"

CLINICAL_FIELD_MAP = {
    "age": "age",
    "gender": "sex",
    "histology": "histology",              # 전원 AD(adenocarcinoma)일 것으로 예상되지만 확인용으로 남겨둠
    "smoking_bi": "smoking (bi",           # Brinkman Index (숫자), "smoking (bi)"의 안전한 prefix
    "stage": "stage (patholog",            # "Stage (Pathological ):" - 괄호/공백 이슈 피하려 앞부분만 사용
    "tnm": "tnm (patholog",
    "cohort_subset": "cohort",             # dataset I/II/III (training/test/additional) 배치 정보
    "relapse": "evidence of relapse",
    "egfr_status": "egfr status",
    "kras_status": "k-ras status",
    "p53_status": "p53 status",
}
VITAL_STATUS_FIELD = "status"              # "Status: Dead"/"Status: Alive"
SURVIVAL_DAYS_KEYWORD = "survival (days"   # "Survival (days): 1326" (일 단위 -> 개월 변환)
EXCLUDE_FIELD_KEYWORD = None                # GSE31210과 달리 별도 제외 플래그 없음
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
        print(f"[{COHORT_NAME}] 위 파일을 raw_data/{COHORT_NAME}/ 안에 넣고 다시 실행하세요. "
              f"(.txt / .txt.gz 둘 다 인식합니다)")
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
