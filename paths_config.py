# -*- coding: utf-8 -*-
"""
paths_config.py
==================
이 파이프라인의 모든 스크립트(cohort_runners/, tcga/, analysis_specific/)가
공유하는 "원본 데이터가 어디 있는지 / 결과물을 어디에 쓸지" 설정.

왜 필요한가
------------
cohort_runners/GSE31210_runner_com.py, analysis_specific/cox_regression_specific.py,
tcga/TCGA_runner_com.py는 전부 다른 폴더에 있습니다. 예전처럼 "스크립트와 같은 폴더에
원본 파일을 둔다"는 규칙이면, 코호트 하나를 처리할 때마다 원본 GEO 파일을 여러 폴더에
복사해 넣어야 하는 번거로움이 생깁니다.

그래서 이 파일은 **딱 한 곳(raw_data/)**을 원본 데이터의 정본 위치로 정하고, 모든
스크립트가 "지금 어느 폴더에서 실행되든" 항상 같은 절대경로를 보도록 만듭니다.
경로 계산 기준이 "현재 작업 디렉토리(cwd)"가 아니라 **이 파일 자신의 위치**이기
때문에, 아래 세 가지 실행 방법이 전부 동일하게 동작합니다:

    cd pipeline_v2 && python cohort_runners/GSE31210_runner_com.py
    cd pipeline_v2/cohort_runners && python GSE31210_runner_com.py
    python /아무데서나/pipeline_v2/cohort_runners/GSE31210_runner_com.py

사용법
------
1. pipeline_v2/raw_data/ 밑에 코호트별 폴더를 만들고 원본 파일을 "한 번만" 넣습니다.
     raw_data/GSE31210/GSE31210_series_matrix.txt, GPL570.txt
     raw_data/GSE72094/GSE72094_series_matrix.txt, GPL15048.txt
     raw_data/GSE36471/GSE36471-GPL9053_series_matrix.txt, GPL9053.txt
     raw_data/TCGA_LUAD/LUAD_expr_matrix_clean.txt, LUAD_matrix.gct, LUAD_clinical.tsv
2. 각 스크립트는 이 파일의 raw()/processed()/results() 함수로 경로를 가져다 쓰므로,
   스크립트를 어디서 실행하든, raw_data 폴더 밑 원본 파일을 두 번 복사할 필요가 없습니다.
3. 원본 데이터를 다른 위치(예: 외장 드라이브, NAS)에 두고 싶다면 아래 RAW_DATA_DIR
   한 줄만 원하는 절대경로로 바꾸면 전체 파이프라인이 그 위치를 보게 됩니다.
"""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent

RAW_DATA_DIR = PROJECT_ROOT / "raw_data"              # 원본 GEO/TCGA 파일 (코호트별 하위 폴더). 필요하면 이 줄만 절대경로로 교체.
DATA_PROCESSED_DIR = PROJECT_ROOT / "data_processed"  # 코호트 러너의 출력 (clinical/expression/merged csv)
RESULTS_DIR = PROJECT_ROOT / "results"                # 분석 스크립트의 출력 (cox/km/clustermap/forest plot)


def raw(cohort_name: str) -> Path:
    """raw_data/<cohort_name>/ 경로를 돌려준다. 폴더가 없으면 안내만 하고 만들지는 않는다
    (원본 데이터를 실수로 빈 폴더로 덮어쓰지 않기 위해)."""
    p = RAW_DATA_DIR / cohort_name
    if not p.exists():
        print(f"[paths_config] 안내: {p} 폴더가 아직 없습니다. "
              f"이 코호트의 원본 파일을 여기에 한 번만 넣어두세요.")
    return p


def processed(cohort_name: str) -> Path:
    """data_processed/<cohort_name>/ 경로. 없으면 자동 생성."""
    p = DATA_PROCESSED_DIR / cohort_name
    p.mkdir(parents=True, exist_ok=True)
    return p


def results(gene_name: str, subfolder: str) -> Path:
    """results/<gene_name>/<subfolder>/ 경로. 없으면 자동 생성."""
    p = RESULTS_DIR / gene_name / subfolder
    p.mkdir(parents=True, exist_ok=True)
    return p
