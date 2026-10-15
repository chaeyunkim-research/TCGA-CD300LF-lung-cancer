# -*- coding: utf-8 -*-
"""
tcga/tcga_common.py  (_com)
==============================
TCGA 데이터(발현매트릭스 + 임상 파일)를 GEO 코호트와 "같은 모양"의
merged_final.csv / merged_valid.csv로 만들어주는 공통 함수.

왜 이게 필요한가: analysis_specific/ 밑의 cox_regression_specific.py,
km_plot_specific.py, coexpression_clustermap_specific.py,
correlation_table_specific.py는 전부 "TARGET_GENE, survival_months,
survival_status 컬럼이 있는 merged_valid.csv"만 있으면 코호트가 무엇이든
그대로 동작하도록 만들어져 있습니다. 즉 TCGA도 이 모양의 csv만 만들면
새 분석 스크립트를 짤 필요 없이 기존 GEO용 _specific 스크립트를 그대로
재사용할 수 있습니다 (TCGA_runner_com.py의 COHORT_RUNS 예시 참고).

TCGA는 GEO series_matrix와 완전히 다른 두 가지 특징이 있어서 GEO용
pipeline_common.py를 그대로 쓸 수 없습니다:
  1. 발현매트릭스가 이미 probe가 아니라 유전자 심볼로 인덱싱되어 있음
     (rerun_PTPN6_tcga.py에서 `expr.loc[GENE_NAME]`으로 바로 접근했던 것과 동일).
     -> probe->gene 매핑 단계가 필요 없음.
  2. 발현매트릭스의 샘플 컬럼명은 TCGA barcode
     (예: "TCGA-05-4244-01A-01R-1107-07")이고, 임상 파일의 환자 ID는
     보통 앞 12자리("TCGA-05-4244")만 있음. 이 둘을 앞 12자리 기준으로
     매칭해야 하고, barcode 14~15번째 자리가 "01"인 원발암(primary tumor)
     샘플만 남겨야 정상조직/재발암 샘플이 섞이지 않습니다.

주의: 실제 TCGA 임상 파일(GDC, cBioPortal, TCGA-CDR 등)의 컬럼명은
출처마다 다릅니다. 아래 함수들은 컬럼명을 전부 인자로 받으므로,
TCGA_runner_com.py 상단 CONFIG에서 실제 파일의 컬럼명에 맞게 지정하세요.
"""

from __future__ import annotations

import pandas as pd


def load_tcga_clinical(
    filepath: str,
    patient_id_col: str,
    vital_status_col: str,
    dead_values: list,
    time_col: str,
    time_unit: str = "days",
    extra_field_map: dict | None = None,
    sep: str = "\t",
) -> pd.DataFrame:
    """TCGA 임상 파일(GDC clinical.tsv, cBioPortal clinical export,
    TCGA-CDR 등 flat table)을 표준 컬럼(patient_id, survival_months,
    survival_status)으로 정리한다.

    Parameters
    ----------
    dead_values : list
        vital_status_col 값 중 "사망"을 의미하는 값들(대소문자 무시).
        예: ["Dead", "1", "DECEASED"]
    time_unit : "days" | "months"
        time_col이 일 단위면 "days"(30.4375로 나눠 변환), 이미 개월이면 "months".
    extra_field_map : dict, optional
        {출력컬럼명: 원본컬럼명} 형태로 age/gender/stage 등 추가 필드도 같이 가져오고 싶을 때.
    """
    raw = pd.read_csv(filepath, sep=sep, dtype=str, low_memory=False)
    print(f"[load_tcga_clinical] {filepath}: {raw.shape[0]}명, 컬럼 {raw.shape[1]}개")

    clinical = pd.DataFrame()
    clinical["patient_id"] = raw[patient_id_col].astype(str).str.strip().str.upper()

    dead_values_upper = {v.upper() for v in dead_values}
    clinical["vital_status_raw"] = raw[vital_status_col]
    clinical["survival_status"] = raw[vital_status_col].apply(
        lambda x: 1 if str(x).strip().upper() in dead_values_upper
        else (0 if pd.notna(x) and str(x).strip() != "" else None)
    )

    time_numeric = pd.to_numeric(raw[time_col], errors="coerce")
    clinical["survival_months"] = time_numeric / 30.4375 if time_unit == "days" else time_numeric

    if extra_field_map:
        for out_col, src_col in extra_field_map.items():
            clinical[out_col] = raw[src_col]

    print(f"[load_tcga_clinical] survival_status 분포: "
          f"{clinical['survival_status'].value_counts(dropna=False).to_dict()}, "
          f"survival_months 결측 아님: {clinical['survival_months'].notna().sum()}")
    return clinical


def extract_panel_expression(expr_matrix_filepath: str, gene_list: list, sep: str = "\t") -> pd.DataFrame:
    """이미 유전자 심볼로 인덱싱된 TCGA 발현매트릭스(예: LUAD_expr_matrix_clean.txt)에서
    gene_list에 해당하는 행만 꺼낸다. (probe 매핑이 필요 없는 TCGA 전용 경로)
    """
    expr = pd.read_csv(expr_matrix_filepath, sep=sep, index_col=0)
    print(f"[extract_panel_expression] {expr_matrix_filepath}: {expr.shape[0]}개 유전자, {expr.shape[1]}개 샘플")

    found = [g for g in gene_list if g in expr.index]
    missing = [g for g in gene_list if g not in expr.index]
    if missing:
        print(f"[extract_panel_expression] 경고 - 매트릭스에 없는 유전자: {missing}")
        for g in missing:
            candidates = [idx for idx in expr.index if g.upper() in str(idx).upper()][:5]
            print(f"  '{g}' 유사 후보: {candidates}")

    panel_expr = expr.loc[found]
    print(f"[extract_panel_expression] 최종 추출: {panel_expr.shape}")
    return panel_expr


def match_and_merge(
    clinical: pd.DataFrame,
    panel_expr: pd.DataFrame,
    barcode_length: int = 12,
    tumor_type_code: str | None = "01",
):
    """clinical(patient_id 컬럼 보유) + panel_expr(유전자 x TCGA barcode)를
    barcode 앞 12자리 기준으로 병합한다.

    tumor_type_code: barcode의 sample-type 자리(14~15번째 문자, 0-index 13:15)가
    이 값인 샘플만 남긴다. "01"=Primary Solid Tumor. None이면 필터링 안 함
    (매트릭스에 이미 종양 샘플만 있는 경우).
    """
    expr_t = panel_expr.T
    expr_t.index.name = "sample_barcode"
    expr_t = expr_t.reset_index()

    if tumor_type_code is not None:
        sample_type = expr_t["sample_barcode"].str.slice(13, 15)
        before = len(expr_t)
        expr_t = expr_t[sample_type == tumor_type_code].copy()
        print(f"[match_and_merge] tumor_type_code='{tumor_type_code}' 필터: {before}개 -> {len(expr_t)}개 샘플")

    expr_t["patient_id"] = expr_t["sample_barcode"].str.slice(0, barcode_length).str.upper()

    dup = expr_t["patient_id"].duplicated().sum()
    if dup:
        print(f"[match_and_merge] 경고: 같은 환자에서 중복 샘플 {dup}건 발견 -> 첫 번째만 사용")
        expr_t = expr_t.drop_duplicates(subset="patient_id", keep="first")

    merged = pd.merge(clinical, expr_t, on="patient_id", how="inner")
    merged_valid = merged.dropna(subset=["survival_status", "survival_months"])

    print(f"[match_and_merge] 임상 {len(clinical)}명 x 발현 {len(expr_t)}건 -> 병합 {len(merged)}건")
    print(f"[match_and_merge] 생존정보 유효한 샘플 수: {len(merged_valid)}")
    return merged, merged_valid
