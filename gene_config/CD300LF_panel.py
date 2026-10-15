# -*- coding: utf-8 -*-
"""
gene_config/CD300LF_panel.py
=============================
CD300LF 연구용 "유전자별(_specific)" 설정.

이 파일 하나만 바꾸면 다른 타겟 유전자 연구(예: PTPN6)로 그대로 전환할 수 있도록,
타겟 유전자 + 같이 볼 시그니처 패널 + probe 매칭용 동의어(alias)를 한 곳에 모아둔다.

SIGNATURE_GENES는 실제로 GSE31210/GSE72094/GSE36471 3개 코호트에서 공통으로
분석했던 "final15" 패널(면역/대식세포 마커 7개 + 증식/기타 마커 7개)이다.
"""

TARGET_GENE = "CD300LF"

# CD300LF를 제외한 co-expression / correlation 비교 대상 14개 유전자
SIGNATURE_GENES = [
    "CD68", "TYROBP", "HLA-DRA", "CXCL10", "TNF", "STAT1", "IL6",
    "CDK1", "CCNB1", "CDKN2A", "MYC", "VIM", "CD274", "PPARG",
]

# 전체 패널 (co-expression clustermap 등에서 TARGET_GENE을 포함해서 볼 때 사용)
FULL_PANEL = [TARGET_GENE] + SIGNATURE_GENES

# 플랫폼 annotation에서 probe를 찾을 때 쓰는 동의어 목록.
# map_GSE31210_genes*.py / map_GSE72094_genes.py / map_GSE68465_genes.py / map_final15.py
# 에 흩어져 있던 alias_map들을 통합 (합집합, 가장 포괄적인 버전 채택).
ALIAS_MAP = {
    "CD300LF": ["CD300LF"],
    "CD68": ["CD68"],
    "TYROBP": ["TYROBP", "DAP12", "KARAP"],
    "HLA-DRA": ["HLA-DRA", "HLADRA", "HLA_DRA"],
    "CXCL10": ["CXCL10", "IP10", "IP-10", "SCYB10"],
    "TNF": ["TNF", "TNFA", "TNFSF2"],
    "STAT1": ["STAT1"],
    "IL6": ["IL6", "IL-6"],
    "CDK1": ["CDK1", "CDC2", "CDC28A", "P34CDC2"],
    "CCNB1": ["CCNB1"],
    "CDKN2A": ["CDKN2A", "P16", "P16INK4A", "MTS1"],
    "MYC": ["MYC", "C-MYC"],
    "VIM": ["VIM"],
    "CD274": ["CD274", "PD-L1", "PDL1", "B7H1"],
    "PPARG": ["PPARG"],
}

# 2-그룹으로 나눠서 볼 때 쓰는 테마 (correlation_table_specific.py의 출력 순서용).
# 원래 correlation_analysis.py에서 "면역 관련" / "증식·기타"로 나눴던 것을 그대로 유지.
THEME_GROUPS = {
    "immune": ["CD68", "TYROBP", "HLA-DRA", "CXCL10", "TNF", "STAT1", "IL6"],
    "proliferation_other": ["CDK1", "CCNB1", "CDKN2A", "MYC", "VIM", "CD274", "PPARG"],
}
