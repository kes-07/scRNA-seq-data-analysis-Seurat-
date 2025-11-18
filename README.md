# scRNA-seq-data-analysis-Seurat-
Single-cell RNA-seq Workflow with Seurat – Fully Automated by Claude_Sonnet_4.5
scRNA_AutoPipeline/
├── README.md                    # 상세 사용 설명서
├── PIPELINE_OVERVIEW.md         # 파이프라인 전체 개요
├── main_pipeline.R              # 메인 실행 스크립트
├── example_usage.R              # 실행 예제
│
├── config/
│   └── config.yaml             # 모든 파라미터 설정
│
├── scripts/                    # 10단계 분석 스크립트
│   ├── 00_setup.R             # 환경 설정
│   ├── 01_data_loading.R      # 데이터 로딩 (10X/GEO/H5/등)
│   ├── 02_quality_control.R   # QC + Doublet 검출
│   ├── 03_normalization.R     # SCT/LogNormalize
│   ├── 04_dimensionality_reduction.R  # PCA/UMAP/tSNE
│   ├── 05_clustering.R        # Leiden/Louvain
│   ├── 06_cell_annotation.R   # Manual + SingleR + Azimuth
│   ├── 07_deg_analysis.R      # Wilcox + MAST + Pseudobulk
│   ├── 08_functional_analysis.R  # GO/KEGG/GSEA/Reactome
│   ├── 09_cellchat.R          # Cell-cell communication
│   └── 10_report_generation.R # HTML/PDF 리포트
│
└── utils/                      # 유틸리티 함수
    ├── plotting_functions.R    # Publication-quality 플롯
    ├── statistical_functions.R # 통계 함수
    └── helper_functions.R      # 헬퍼 함수

    ⭐ 주요 기능
1. 다양한 입력 형식 지원

✅ 10X Genomics
✅ GEO accession
✅ Seurat objects
✅ H5 files
✅ Matrix files (MTX)

2. High-Impact Journal 기준 분석

✅ SCTransform (Hafemeister & Satija, Genome Biology 2019)
✅ Pseudobulk DESeq2 (Squair et al., Nature Communications 2021) ← 가장 권장
✅ CellChat (Jin et al., Nature Communications 2021)
✅ GSEA (Subramanian et al., PNAS 2005)

3. DEG 분석 - 3가지 방법

✅ Wilcoxon - 빠르고 robust
✅ MAST - Dropout 고려
✅ DESeq2 Pseudobulk - High-impact journal 추천 방법!

4. Functional Analysis

✅ Gene Ontology (BP/MF/CC)
✅ KEGG Pathway
✅ GSEA (Hallmark, C2, C5 gene sets)
✅ Reactome Pathway

5. Cell Type Annotation - 이중 검증

✅ Manual (사용자 정의 마커)
✅ SingleR (자동)
✅ Azimuth (자동)
✅ scType (자동)
→ 비교하여 최적의 annotation 선택 가능!

6. 자동화 기능

✅ Checkpoint/Resume 기능
✅ Optimal PC 자동 선택 (Elbow/JackStraw)
✅ Optimal Resolution 자동 선택
✅ 병렬 처리 지원
✅ 상세한 로깅

1단계: Config 파일 수정
# config/config.yaml 편집
project:
  name: "My_Retina_Analysis"
  output_dir: "./results"

input:
  data_type: "10X"
  data_path: "/path/to/your/data"

species:
  organism: "mouse"

2단계: 실행
# R에서
source("main_pipeline.R")

# 또는 터미널에서
Rscript main_pipeline.R config/config.yaml
```

#### **3단계: 결과 확인**
```
results/
├── final_seurat_object.rds      # 최종 분석 객체
├── 07_DEGs_wilcox.csv           # Wilcoxon DEGs
├── 07_DEGs_MAST.csv             # MAST DEGs
├── 07_DEGs_DESeq2_pseudobulk.csv  # Pseudobulk DEGs ⭐
├── 08_GO_enrichment.csv         # GO 분석
├── 08_KEGG_enrichment.csv       # KEGG 분석
├── 08_GSEA_results.csv          # GSEA 결과
├── 09_cellchat_object.rds       # CellChat 객체
└── scRNA_Analysis_Report.html   # 종합 리포트

## High-Impact Journal 권장사항
Nature/Cell/Science 수준 분석을 위해:

✅ Pseudobulk DESeq2 사용 (필수!)
✅ SCTransform normalization
✅ Resolution optimization
✅ Multiple annotation methods 비교
✅ Comprehensive functional analysis
✅ Publication-quality figures (300 DPI)

## 다음 단계

파이프라인 다운로드:

압축 파일: scRNA_AutoPipeline.tar.gz


Config 파일 수정:

형님의 데이터 경로 입력
Retinal cell markers 확인/수정
Comparison groups 설정


실행:
Rscript main_pipeline.R config_retina.yaml
