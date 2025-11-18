# Automated scRNA-seq Analysis Pipeline

A comprehensive, publication-ready single-cell RNA-seq analysis pipeline based on best practices from high-impact journals (Nature, Cell, Science).

## 📋 Overview

This automated pipeline performs end-to-end scRNA-seq analysis including:

- **Quality Control** (scDblFinder/DoubletFinder)
- **Normalization** (SCTransform/LogNormalize)
- **Dimensionality Reduction** (PCA, UMAP, tSNE)
- **Clustering** (Leiden/Louvain with resolution optimization)
- **Cell Type Annotation** (Manual + Automated: SingleR, Azimuth, scType)
- **Differential Expression** (Wilcoxon, MAST, DESeq2 pseudobulk)
- **Functional Enrichment** (GO, KEGG, GSEA, Reactome)
- **Cell-Cell Communication** (CellChat)
- **Automated Reporting** (HTML/PDF)

## 🎯 Key Features

### Publication-Ready Analysis
- **Methods based on high-impact publications**:
  - Hao et al. Cell 2021 (Seurat v4)
  - Squair et al. Nature Communications 2021 (Pseudobulk DEG)
  - Hafemeister & Satija, Genome Biology 2019 (SCTransform)
  - Jin et al. Nature Communications 2021 (CellChat)

### Multiple Input Formats
- 10X Genomics (CellRanger output)
- GEO accession numbers
- Seurat objects (.rds)
- H5 files
- Matrix files (MTX format)

### Automated & Reproducible
- YAML-based configuration
- Checkpoint/resume capability
- Comprehensive logging
- All parameters documented

### Parallel Processing
- Multi-core support
- Memory-optimized operations
- Progress tracking

## 📦 Installation

### Prerequisites

```r
# Required R version
R >= 4.1.0

# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

### Install Pipeline

```bash
# Clone repository
git clone https://github.com/yourusername/scRNA_AutoPipeline.git
cd scRNA_AutoPipeline

# Run setup (installs all dependencies)
Rscript scripts/00_setup.R
```

### Manual Package Installation

If automatic installation fails:

```r
# Core packages
install.packages(c("Seurat", "tidyverse", "yaml", "here", "logger"))

# Bioconductor packages
BiocManager::install(c("scDblFinder", "SingleR", "celldex", "DESeq2", 
                       "clusterProfiler", "org.Mm.eg.db", "org.Hs.eg.db"))

# GitHub packages
devtools::install_github("sqjin/CellChat")
```

## 🚀 Quick Start

### 1. Configure Your Analysis

Edit `config/config.yaml`:

```yaml
project:
  name: "My_scRNA_Analysis"
  output_dir: "./results"

input:
  data_type: "10X"  # Options: "10X", "GEO", "Seurat", "H5", "MTX"
  data_path: "/path/to/your/data"

species:
  organism: "mouse"  # or "human"

# ... (see config.yaml for all options)
```

### 2. Run the Pipeline

```r
# From R console
source("main_pipeline.R")

# Or from command line
Rscript main_pipeline.R config/config.yaml
```

### 3. View Results

Results will be in your specified `output_dir`:
```
results/
├── 01_raw_data.rds
├── 02_qc_filtered.rds
├── 02_QC_plots.pdf
├── 03_normalized.rds
├── 04_dimred_plots.pdf
├── 05_clustered.rds
├── 06_annotated.rds
├── 06_marker_genes.csv
├── 07_DEGs_wilcox.csv
├── 07_DEGs_MAST.csv
├── 07_DEGs_DESeq2_pseudobulk.csv
├── 08_GO_enrichment.csv
├── 08_KEGG_enrichment.csv
├── 09_cellchat_object.rds
├── final_seurat_object.rds
└── scRNA_Analysis_Report.html
```

## 📚 Documentation

### Configuration Guide

#### Input Options

```yaml
input:
  # For 10X data
  data_type: "10X"
  data_path: "/path/to/cellranger/outs/filtered_feature_bc_matrix"
  
  # For GEO data
  data_type: "GEO"
  geo_accession: "GSE123456"
  
  # For multiple samples
  samples:
    - name: "Control"
      path: "/path/to/control"
      condition: "Control"
    - name: "Treatment"
      path: "/path/to/treatment"
      condition: "Treatment"
```

#### Quality Control

```yaml
quality_control:
  min_features: 200        # Minimum genes per cell
  max_features: 7500       # Maximum genes per cell
  mt_percent_max: 10       # Maximum mitochondrial %
  doublet_detection: TRUE
  doublet_method: "scDblFinder"  # or "DoubletFinder"
```

#### Normalization

```yaml
normalization:
  method: "SCT"  # or "LogNormalize"
  
  # SCTransform settings
  sct_params:
    vars_to_regress: ["percent.mt"]
    variable_features: 3000
```

#### Clustering

```yaml
clustering:
  algorithm: "leiden"  # or "louvain"
  test_resolutions: TRUE
  resolutions: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
```

#### Cell Type Annotation

```yaml
annotation:
  methods:
    manual: TRUE
    SingleR: TRUE
    Azimuth: TRUE
    scType: TRUE
  
  # Manual markers (for retinal cells example)
  manual_markers:
    Rods: ["Rho", "Nr2e3", "Nrl"]
    Cones: ["Opn1sw", "Opn1mw", "Arr3"]
    Muller: ["Rlbp1", "Glul", "Sox9"]
    # Add your cell type markers here
```

#### Differential Expression

```yaml
differential_expression:
  methods: ["wilcox", "MAST", "DESeq2_pseudobulk"]
  
  # Pseudobulk settings (RECOMMENDED for publications)
  pseudobulk:
    use: TRUE
    aggregate_by: "sample"
    min_cells: 10
  
  # Thresholds
  logfc_threshold: 0.25
  pval_cutoff: 0.05
```

## 🔬 Advanced Usage

### Resume from Checkpoint

```yaml
advanced:
  continue_from_step: 5  # Resume from Step 5 (Clustering)
```

### Custom Gene Sets

```yaml
functional_analysis:
  custom_genesets:
    use: TRUE
    gmt_file: "/path/to/custom_genesets.gmt"
```

### Memory Management

```yaml
computation:
  n_cores: 8
  memory_limit: 32  # GB
  future_globals: 4000  # MB
```

## 📊 Output Files

### Seurat Objects
- `01_raw_data.rds` - Raw data after loading
- `02_qc_filtered.rds` - After quality control
- `03_normalized.rds` - After normalization
- `04_dimred.rds` - After dimensionality reduction
- `05_clustered.rds` - After clustering
- `06_annotated.rds` - After cell type annotation
- `final_seurat_object.rds` - Complete analysis

### Data Tables
- `06_marker_genes.csv` - Marker genes for each cluster
- `07_DEGs_[method].csv` - Differential expression results
- `08_GO_enrichment.csv` - GO enrichment results
- `08_KEGG_enrichment.csv` - KEGG pathway results
- `08_GSEA_results.csv` - GSEA results

### Plots
- `02_QC_plots.pdf` - Quality control visualizations
- `03_normalization_plots.pdf` - Normalization diagnostics
- `04_dimred_plots.pdf` - PCA, UMAP, tSNE plots
- `05_clustering_plots.pdf` - Clustering visualizations
- `06_annotation_plots.pdf` - Cell type annotation
- `07_DEG_plots.pdf` - Volcano plots, heatmaps
- `08_enrichment_plots.pdf` - Pathway enrichment
- `09_cellchat_plots.pdf` - Cell-cell communication

### Reports
- `scRNA_Analysis_Report.html` - Comprehensive HTML report
- `pipeline_log.txt` - Detailed execution log

## 🛠️ Troubleshooting

### Common Issues

**1. Memory errors**
```yaml
# Reduce memory usage
computation:
  n_cores: 4  # Reduce cores
  
quality_control:
  max_features: 5000  # Stricter filtering
```

**2. Missing packages**
```r
# Re-run setup
source("scripts/00_setup.R")
setup_environment(config)
```

**3. GEO download fails**
```yaml
# Use manual download
input:
  data_type: "10X"  # Change to local data
  data_path: "/path/to/downloaded/data"

```

And the key methods:
- Seurat: Hao et al. Cell 2021
- SCTransform: Hafemeister & Satija, Genome Biology 2019
- Pseudobulk DEG: Squair et al. Nature Communications 2021
- CellChat: Jin et al. Nature Communications 2021

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## 🙏 Acknowledgments

This pipeline integrates best practices from:
- Satija Lab (Seurat framework)
- Broad Institute (single-cell methods)
- Multiple high-impact publications

---

**Version:** 1.0.0
**Last Updated:** November 2025
