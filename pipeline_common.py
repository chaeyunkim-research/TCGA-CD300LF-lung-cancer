# -*- coding: utf-8 -*-
"""
pipeline_common.py  (공통 엔진, "_com")
=========================================
GEO 시리즈(series matrix) 기반 발현/임상 데이터 파이프라인의 코호트-무관 공통 로직.

이 파일은 "어떤 유전자를 보고 있는가"와 "어떤 GEO 코호트를 보고 있는가"에
관계없이 항상 똑같이 동작해야 하는 부분만 모아둔 것입니다.
- 코호트마다 달라지는 값(파일명, 플랫폼 symbol 컬럼명, 임상 필드 키워드 등)은
  cohort_runners/<GSE번호>_runner_com.py 안의 COHORT config 딕셔너리로 넘깁니다.
- 유전자 패널(타겟 유전자 + 같이 볼 시그니처 유전자들 + 동의어)은
  gene_config/<유전자명>_panel.py 로 분리되어 있습니다.

기존에 흩어져 있던 스크립트들과의 대응 관계 (레거시 -> 이 파일의 함수):
  parse_geo_metadata.py, parse_GSE31210_metadata.py, parse_GSE72094_metadata.py,
  parse_GSE68465_metadata.py                          -> parse_series_matrix_metadata()
  parse_GSE31210_clinical.py (v1~v3), parse_GSE72094_clinical.py,
  parse_GSE68465_clinical.py, correlation_analysis.py의 1번 섹션
                                                        -> extract_raw_characteristics(), build_clinical()
  find_target_genes*.py, map_final15.py, map_GSE31210_genes.py(v1~v4),
  map_GSE72094_genes.py, map_GSE68465_genes.py, check_gpl15048_structure.py
                                                        -> find_platform_header_line(), map_probes_to_genes()
  extract_expression.py, extract_GSE31210_expression.py, extract_GSE72094_expression.py
                                                        -> extract_gene_expression()
  merge_GSE31210.py, merge_GSE72094.py                 -> merge_clinical_expression()

각 GSE31210 map_*_genes 버전(v1~v4)과 parse_GSE31210_clinical 버전(v1~v3)을
비교해서, 실제로 최종 결과(HR=1.100, N=204 등)를 만들어낸 가장 견고한 로직만
아래에 채택했습니다:
  - 플랫폼 주석 파일 헤더 찾기: "!platform_table_begin" 마커가 있으면 그걸 쓰고,
    없으면(GPL570처럼) '#'으로 시작하는 줄을 건너뛰고 첫 데이터 줄을 헤더로 씀 (v3/v4 방식)
  - 유전자 심볼 매칭: "///"로 여러 심볼이 묶여 있는 경우까지 처리 (v4 방식)
  - 임상 필드 파싱: 샘플마다 characteristics 필드 "순서"가 다를 수 있다고 가정하고,
    각 샘플별로 "key: value"를 직접 파싱해서 keyword로 찾는 방식 (v2/v3 방식).
    generation time 이 아니라 실제 GSE31210_merged_valid.csv의 survival_months 값과
    대조해서 이 로직(=v3)이 맞다는 것을 확인함 (예: 743일 / 30.4375 = 24.41066... 일치).
"""

from __future__ import annotations

import gzip

import pandas as pd


def _smart_open(filepath: str, mode: str = "r", encoding: str = "utf-8", errors: str = "replace"):
    """파일명이 .gz로 끝나면 gzip으로, 아니면 일반 텍스트로 연다.

    GSE13213처럼 원본을 압축 해제하지 않고 series_matrix.txt.gz / GPLxxxx.txt.gz를
    그대로 raw_data/에 넣어도 동작하도록 하기 위함 (pd.read_csv는 이미 .gz를 자동
    인식하지만, 메타데이터/헤더를 한 줄씩 읽는 순수 open()은 그렇지 않아서 추가함).
    """
    if str(filepath).endswith(".gz"):
        return gzip.open(filepath, mode + "t", encoding=encoding, errors=errors)
    return open(filepath, mode, encoding=encoding, errors=errors)


# ---------------------------------------------------------------------------
# 1. Series matrix 메타데이터 파싱
# ---------------------------------------------------------------------------

def parse_series_matrix_metadata(filepath: str):
    """GEO series_matrix.txt 파일의 '!Sample_...' 헤더 라인들을 파싱한다.

    Returns
    -------
    meta_df : pd.DataFrame
        컬럼 = "!Sample_geo_accession", "!Sample_characteristics_ch1",
        "!Sample_characteristics_ch1__1", ... 등. (같은 key가 여러 번 나오면 __N 접미사)
        행 = 샘플.
    table_start_line : int
        실제 발현 매트릭스(!series_matrix_table_begin 다음 줄)가 시작하는 줄 번호.
        pd.read_csv(..., skiprows=table_start_line) 에 바로 쓸 수 있음.
    """
    sample_meta = {}
    table_start_line = None

    with _smart_open(filepath) as f:
        for i, line in enumerate(f):
            if line.startswith("!Sample_"):
                parts = line.strip().split("\t")
                key = parts[0].replace('"', "")
                values = [v.replace('"', "") for v in parts[1:]]
                existing_count = len(
                    [k for k in sample_meta if k == key or k.startswith(key + "__")]
                )
                unique_key = key if existing_count == 0 else f"{key}__{existing_count}"
                sample_meta[unique_key] = values
            if line.startswith("!series_matrix_table_begin"):
                table_start_line = i + 1
                break

    meta_df = pd.DataFrame(sample_meta)
    return meta_df, table_start_line


def extract_raw_characteristics(meta_df: pd.DataFrame) -> pd.DataFrame:
    """!Sample_characteristics_ch1* 컬럼들을 샘플별로 "key: value" 단위로 풀어서
    geo_accession + 각 key를 컬럼으로 갖는 넓은(wide) 테이블을 만든다.

    샘플마다 characteristics 필드가 나열되는 "순서"가 다를 수 있다는 전제 하에
    (실제로 GSE31210이 그랬음), 컬럼 위치가 아니라 각 샘플의 텍스트를 직접
    "key: value"로 쪼개서 매핑한다. -> parse_GSE31210_clinical_v2/v3 방식.
    """
    char_cols = [c for c in meta_df.columns if c.startswith("!Sample_characteristics_ch1")]
    n_samples = len(meta_df)
    records = []

    for i in range(n_samples):
        record = {"geo_accession": meta_df["!Sample_geo_accession"].iloc[i]}
        for col in char_cols:
            val = str(meta_df[col].iloc[i])
            if ":" in val:
                key_part, value_part = val.split(":", 1)
                record[key_part.strip().lower()] = value_part.strip()
        records.append(record)

    return pd.DataFrame(records)


def get_field_by_keyword(df: pd.DataFrame, keyword: str) -> pd.Series:
    """컬럼명이 keyword로 시작하는 첫 번째 컬럼을 반환. 없으면 전부 None인 Series."""
    matching_cols = [c for c in df.columns if c.startswith(keyword)]
    if not matching_cols:
        return pd.Series([None] * len(df))
    return df[matching_cols[0]]


# ---------------------------------------------------------------------------
# 2. 임상 데이터(생존정보 포함) 표준화
# ---------------------------------------------------------------------------

def build_clinical(
    raw_df: pd.DataFrame,
    field_keyword_map: dict,
    vital_status_field: str | None = None,
    survival_status_keyword: str | None = None,
    survival_days_keyword: str | None = None,
    survival_months_keyword: str | None = None,
    exclude_field_keyword: str | None = None,
    exclude_value: str = "exclude",
) -> pd.DataFrame:
    """extract_raw_characteristics() 결과를 표준 임상 테이블로 정리한다.

    Parameters
    ----------
    field_keyword_map : dict
        {출력컬럼명: raw_df에서 찾을 키워드}. 예: {"age": "age (years)", "gender": "gender"}
        geo_accession, vital_status, survival 관련은 별도 파라미터로 처리하므로 여기 넣지 않는다.
    vital_status_field : str, optional
        raw_df에서 생존 여부(Dead/Alive 텍스트)를 담은 필드를 찾을 키워드.
        GSE31210/GSE72094/GSE68465처럼 텍스트로 제공되는 경우 이걸 쓴다.
    survival_status_keyword : str, optional
        vital_status_field 대신, survival_status(0/1)가 이미 숫자로 제공되는 코호트용
        (예: GSE36471은 characteristics에 "survival_status: 0/1"이 직접 들어있음).
        둘 중 정확히 하나만 지정해야 한다.
    survival_days_keyword / survival_months_keyword : str
        둘 중 하나만 지정. days로 주어지면 30.4375로 나눠 개월로 통일한다.
        (GSE72094, GSE31210(v3 기준)은 days만 존재 -> 변환 필요.
         GSE68465/GSE36471은 이미 개월 단위 필드를 직접 제공.)
    exclude_field_keyword : str, optional
        "exclude for prognosis" 같은 제외 플래그 필드가 있으면 지정.
        exclude_value(기본 "exclude")와 일치하는 행은 걸러낸다. (GSE31210 전용 케이스)
    """
    if bool(vital_status_field) == bool(survival_status_keyword):
        raise ValueError(
            "vital_status_field 또는 survival_status_keyword 중 정확히 하나만 지정하세요."
        )

    clinical = pd.DataFrame()
    clinical["geo_accession"] = raw_df["geo_accession"]

    for out_col, keyword in field_keyword_map.items():
        clinical[out_col] = get_field_by_keyword(raw_df, keyword)

    if vital_status_field:
        clinical["vital_status"] = get_field_by_keyword(raw_df, vital_status_field)
        clinical["survival_status"] = clinical["vital_status"].apply(
            lambda x: 1 if str(x).strip().lower() == "dead"
            else (0 if str(x).strip().lower() == "alive" else None)
        )
    else:
        clinical["survival_status"] = pd.to_numeric(
            get_field_by_keyword(raw_df, survival_status_keyword), errors="coerce"
        )

    if survival_days_keyword:
        clinical["survival_days"] = pd.to_numeric(
            get_field_by_keyword(raw_df, survival_days_keyword), errors="coerce"
        )
        clinical["survival_months"] = clinical["survival_days"] / 30.4375
    elif survival_months_keyword:
        clinical["survival_months"] = pd.to_numeric(
            get_field_by_keyword(raw_df, survival_months_keyword), errors="coerce"
        )
    else:
        raise ValueError("survival_days_keyword 또는 survival_months_keyword 중 하나는 필수입니다.")

    if exclude_field_keyword:
        clinical["exclude_flag"] = get_field_by_keyword(raw_df, exclude_field_keyword)
        before = len(clinical)
        clinical = clinical[
            clinical["exclude_flag"].astype(str).str.strip().str.lower() != exclude_value
        ].copy()
        print(f"[build_clinical] exclude_flag 적용: {before}개 -> {len(clinical)}개")

    print(f"[build_clinical] 최종 샘플 수: {len(clinical)}, "
          f"survival_status 분포: {clinical['survival_status'].value_counts(dropna=False).to_dict()}, "
          f"survival_months 결측 아님: {clinical['survival_months'].notna().sum()}")

    return clinical


# ---------------------------------------------------------------------------
# 3. 플랫폼 annotation -> probe/gene 매핑
# ---------------------------------------------------------------------------

def find_platform_header_line(filepath: str) -> int:
    """플랫폼 annotation(GPLxxxxx.txt) 파일에서 실제 데이터 헤더가 시작하는 줄 번호를 찾는다.

    1순위: "!platform_table_begin" 마커 (GPL9053, GPL15048, GPL96 등 대부분의 GEO 플랫폼 파일)
    2순위: 마커가 없는 파일(GPL570처럼 순수 주석파일 형태로 받은 경우)은
           '#'으로 시작하는 설명줄을 건너뛰고 첫 비-# 줄을 헤더로 사용.
    """
    with _smart_open(filepath) as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.startswith("!platform_table_begin"):
            j = i + 1
            while j < len(lines) and lines[j].startswith("#"):
                j += 1
            return j

    for i, line in enumerate(lines):
        if not line.startswith("#"):
            return i

    return 0


def _detect_symbol_column(columns, symbol_col: str | None):
    if symbol_col:
        if symbol_col not in columns:
            raise KeyError(
                f"지정한 symbol_col='{symbol_col}' 이 컬럼 목록에 없습니다: {list(columns)}"
            )
        return symbol_col

    # 우선순위: 정확히 알려진 이름들 -> "symbol"을 포함하는 컬럼 -> ORF(오래된 GEO 플랫폼 관례)
    known_exact = ["Gene Symbol", "GeneSymbol", "gene_symbol"]
    for name in known_exact:
        if name in columns:
            return name
    candidates = [c for c in columns if "symbol" in c.lower()]
    if candidates:
        return candidates[0]
    if "ORF" in columns:
        return "ORF"
    raise KeyError(
        f"gene symbol 컬럼을 자동으로 찾지 못했습니다. symbol_col을 직접 지정하세요. 컬럼 목록: {list(columns)}"
    )


def _symbol_matches(symbol_field, target_aliases_upper: set) -> bool:
    """'///'로 여러 심볼이 묶인 필드까지 처리 (map_GSE31210_genes_v4.py 방식을 전체 공통화)."""
    if pd.isna(symbol_field):
        return False
    symbols = [s.strip().upper() for s in str(symbol_field).split("///")]
    return any(s in target_aliases_upper for s in symbols)


def map_probes_to_genes(
    platform_filepath: str,
    alias_map: dict,
    symbol_col: str | None = None,
    id_col: str = "ID",
) -> pd.DataFrame:
    """플랫폼 annotation 파일에서 alias_map에 정의된 유전자들의 probe를 찾아
    표준 유전자명(GENE_STANDARD) 컬럼을 붙인 매핑 테이블을 만든다.

    alias_map: {표준유전자명: [동의어1, 동의어2, ...]} (gene_config/*_panel.py 에서 가져옴)
    """
    header_line = find_platform_header_line(platform_filepath)
    annot = pd.read_csv(
        platform_filepath, sep="\t", skiprows=header_line,
        dtype=str, on_bad_lines="skip", engine="python",
    )
    print(f"[map_probes_to_genes] {platform_filepath}: probe {len(annot)}개, "
          f"컬럼: {annot.columns.tolist()}")

    resolved_symbol_col = _detect_symbol_column(annot.columns, symbol_col)
    print(f"[map_probes_to_genes] symbol 컬럼으로 '{resolved_symbol_col}' 사용")

    all_matches = []
    report = {}
    for gene, aliases in alias_map.items():
        alias_set = {a.upper() for a in aliases}
        mask = annot[resolved_symbol_col].apply(lambda x: _symbol_matches(x, alias_set))
        m = annot[mask].copy()
        if len(m) > 0:
            m["GENE_STANDARD"] = gene
            all_matches.append(m)
        report[gene] = len(m)

    print("[map_probes_to_genes] 유전자별 매칭 probe 수:")
    for gene, count in report.items():
        print(f"  {'OK' if count else 'MISSING':7s} {gene}: {count}")
    missing = [g for g, c in report.items() if c == 0]
    if missing:
        print(f"[map_probes_to_genes] 경고 - 못 찾은 유전자: {missing}")

    combined = pd.concat(all_matches, ignore_index=True) if all_matches else pd.DataFrame()
    if id_col not in combined.columns:
        raise KeyError(
            f"id_col='{id_col}' 이 annotation 컬럼에 없습니다. 실제 컬럼: {annot.columns.tolist()}"
        )
    return combined


# ---------------------------------------------------------------------------
# 4. 발현 매트릭스 추출 (probe -> gene 평균)
# ---------------------------------------------------------------------------

def extract_gene_expression(
    series_matrix_filepath: str,
    probe_map: pd.DataFrame,
    probe_id_col: str = "ID",
) -> pd.DataFrame:
    """series_matrix 발현 매트릭스에서 probe_map에 있는 probe만 남기고,
    같은 유전자로 매핑된 probe는 평균내서 유전자 x 샘플 매트릭스를 만든다.
    """
    _, table_start_line = parse_series_matrix_metadata(series_matrix_filepath)
    expr = pd.read_csv(
        series_matrix_filepath, sep="\t", skiprows=table_start_line, low_memory=False
    )
    id_col = expr.columns[0]
    print(f"[extract_gene_expression] 발현 매트릭스 shape: {expr.shape}")

    expr[id_col] = expr[id_col].astype(str)
    probe_map = probe_map.copy()
    probe_map[probe_id_col] = probe_map[probe_id_col].astype(str)

    id_to_gene = dict(zip(probe_map[probe_id_col], probe_map["GENE_STANDARD"]))
    expr_filtered = expr[expr[id_col].isin(id_to_gene.keys())].copy()
    expr_filtered["GENE"] = expr_filtered[id_col].map(id_to_gene)

    print(f"[extract_gene_expression] 필터링된 probe 수: {len(expr_filtered)} "
          f"(매핑 테이블: {len(probe_map)}개 probe, {probe_map['GENE_STANDARD'].nunique()}개 유전자)")

    sample_cols = [c for c in expr.columns if c != id_col]
    gene_expr = expr_filtered.groupby("GENE")[sample_cols].mean()
    print(f"[extract_gene_expression] 최종 shape: {gene_expr.shape}, 유전자: {gene_expr.index.tolist()}")
    return gene_expr


# ---------------------------------------------------------------------------
# 5. 임상 + 발현 병합
# ---------------------------------------------------------------------------

def merge_clinical_expression(clinical: pd.DataFrame, gene_expr: pd.DataFrame):
    """clinical(샘플 단위 행) + gene_expr(유전자 x 샘플)을 geo_accession 기준으로 병합.

    Returns (merged_final, merged_valid) — merged_valid는 survival_status/survival_months가
    모두 있는 행만 남긴 것 (KM/Cox 분석용).
    """
    expr_t = gene_expr.T
    expr_t.index.name = "geo_accession"
    expr_t = expr_t.reset_index()

    merged = pd.merge(clinical, expr_t, on="geo_accession", how="inner")
    merged_valid = merged.dropna(subset=["survival_status", "survival_months"])

    print(f"[merge_clinical_expression] 전체 병합 샘플 수: {len(merged)}")
    print(f"[merge_clinical_expression] 생존정보 유효한 샘플 수: {len(merged_valid)}")
    return merged, merged_valid
