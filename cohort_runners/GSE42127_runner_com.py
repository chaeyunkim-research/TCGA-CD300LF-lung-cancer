# -*- coding: utf-8 -*-
"""
cohort_runners/GSE42127_runner_com.py   ("_com" = 코호트 파이프라인, 공통 엔진 사용)
==========================================================================
GSE42127 (MD Anderson Cancer Center, NSCLC, GPL6884 Illumina HumanWG-6 v3)
코호트를 raw GEO 파일 -> 표준화된 clinical/expression/merged csv까지 만드는 실행 스크립트.

실제 series_matrix.txt.gz 176명 전체 메타데이터로 clinical 파싱을 미리 검증했습니다.
확인한 내용 (기존 코호트와 다른 점 위주):
  - 이 코호트는 ADC(133) + Squamous(43) 혼합입니다. LUAD 검증 트랙이므로 이 러너에서
    histology == "Adenocarcinoma"인 행만 남기고 SCC(43예)는 제외합니다 (아래
    "ADC subset" 단계). SCC 43예는 나중에 LUSC 검증할 때 재사용할 수 있습니다.
  - vital status가 "Alive"/"Dead" 전체 단어가 아니라 "A"/"D" 약어로 들어있어서,
    build_clinical()에 넘기기 전에 이 러너에서 미리 "Dead"/"Alive"로 변환합니다
    (pipeline_common.build_clinical()은 "dead"/"alive" 전체 단어만 인식하므로).
  - survival은 "overall survival months"로 이미 개월 단위 제공 (일 단위 변환 불필요).
  - "had_adjuvant_chemo" 필드가 있음 (TRUE/FALSE) — 무작위 배정 임상시험이 아니라
    관찰 코호트라서 GSE14814(JBR.10)처럼 강한 치료군 교란은 아니지만, 필요하면
    Cox 모델에 공변량으로 추가 고려 가능 (지금은 clinical에 컬럼만 남겨둠).
  - ADC subset 후 133명, survival_status 분포 Alive 90 / Dead 43, survival_months
    전원 유효함을 실제 데이터로 확인함.

주의: 아직 probe 매핑~발현 추출 단계는 GPL6884.txt(Illumina annotation) 없이는
검증 못 함 — 처음 실행 시 콘솔의 symbol 컬럼 자동감지 결과와 MISSING 유전자를 확인하세요
(Illumina 플랫폼은 흔히 "Symbol"이나 "ILMN_Gene" 컬럼을 씁니다 — 자동감지 실패하면
알려주시면 SYMBOL_COL을 직접 지정해드립니다).

실행 전 필요한 원본 파일 (raw_data/GSE42127/ 에 둘 것 — paths_config.py 참고):
  GSE42127_series_matrix.txt(.gz), GPL6884.txt(.gz)
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import pipeline_common as pc
import paths_config as pcfg

# ===== 여기만 바꾸면 됩니다 =====
COHORT_NAME = "GSE42127"
GENE_CONFIG_MODULE = "gene_config.CD300LF_panel"

DATA_DIR = pcfg.raw(COHORT_NAME)       # raw_data/GSE42127/
OUT_DIR = pcfg.processed(COHORT_NAME)  # data_processed/GSE42127/


def _find_input(basename_no_ext: str) -> Path:
    for suffix in (".txt", ".txt.gz"):
        p = DATA_DIR / f"{basename_no_ext}{suffix}"
        if p.exists():
            return p
    return DATA_DIR / f"{basename_no_ext}.txt"


SERIES_MATRIX_FILE = _find_input("GSE42127_series_matrix")
PLATFORM_FILE = _find_input("GPL6884")
SYMBOL_COL = None          # None = 자동 감지 (Illumina는 보통 "Symbol" 계열 컬럼)
PROBE_ID_COL = "ID"

CLINICAL_FIELD_MAP = {
    "age": "age at surgery",
    "gender": "gender",
    "stage": "final.pat.stage",
    "histology": "histology",           # ADC subset 확인용으로 남겨둠
    "adjuvant_chemo": "had_adjuvant_chemo",
}
VITAL_STATUS_FIELD = "survival status"       # 전처리 후 "Dead"/"Alive"
SURVIVAL_MONTHS_KEYWORD = "overall survival months"   # 이미 개월 단위
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

    # --- 이 코호트만의 전처리 두 가지 ---
    # (a) LUAD 검증 트랙이므로 ADC만 남긴다 (SCC 43예는 LUSC용으로 나중에 재사용 가능)
    n_before = len(raw_char)
    raw_char = raw_char[raw_char["histology"].str.strip().str.lower() == "adenocarcinoma"].copy()
    print(f"[{COHORT_NAME}] histology ADC subset: {n_before} -> {len(raw_char)} "
          f"(Squamous {n_before - len(raw_char)}예 제외)")
    # (b) vital status 약어("D"/"A") -> build_clinical이 인식하는 "Dead"/"Alive"로 변환
    raw_char["survival status"] = raw_char["survival status"].replace({"D": "Dead", "A": "Alive"})

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
