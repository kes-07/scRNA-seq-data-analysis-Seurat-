# scRNA-seq Automated Analysis Pipeline - Complete Overview

## 📁 File Structure

```
scRNA_AutoPipeline/
│
├── README.md                          # Main documentation
├── main_pipeline.R                    # Main orchestration script
├── example_usage.R                    # Example configurations
│
├── config/
│   └── config.yaml                    # Main configuration file
│
├── scripts/                           # Analysis step scripts
│   ├── 00_setup.R                    # Environment setup
│   ├── 01_data_loading.R             # Data input (10X/GEO/H5/etc)
│   ├── 02_quality_control.R          # QC + doublet detection
│   ├── 03_normalization.R            # SCT/LogNormalize + integration
│   ├── 04_dimensionality_reduction.R # PCA/UMAP/tSNE
│   ├── 05_clustering.R               # Leiden/Louvain clustering
│   ├── 06_cell_annotation.R          # Manual + automated annotation
│   ├── 07_deg_analysis.R             # Wilcox/MAST/pseudobulk DESeq2
│   ├── 08_functional_analysis.R      # GO/KEGG/GSEA/Reactome
│   ├── 09_cellchat.R                 # Cell-cell communication
│   └── 10_report_generation.R        # HTML/PDF report
│
└── utils/                             # Utility functions
    ├── plotting_functions.R          # Publication-quality plots
    ├── statistical_functions.R       # Statistical methods
    └── helper_functions.R            # Helper utilities
```

## 🔄 Pipeline Workflow

```
┌─────────────────────┐
│  Step 0: Setup      │ ← Install packages, configure environment
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│  Step 1: Load Data  │ ← 10X/GEO/H5/Seurat/MTX formats
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│  Step 2: QC         │ ← Filter cells, detect doublets
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│  Step 3: Normalize  │ ← SCTransform or LogNormalize
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│  Step 4: DimRed     │ ← PCA, UMAP, tSNE
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│  Step 5: Cluster    │ ← Leiden/Louvain with optimization
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│  Step 6: Annotate   │ ← Manual + SingleR + Azimuth
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│  Step 7: DEG        │ ← Wilcoxon + MAST + Pseudobulk
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│  Step 8: Enrichment │ ← GO + KEGG + GSEA + Reactome
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│  Step 9: CellChat   │ ← Cell-cell communication
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│  Step 10: Report    │ ← HTML/PDF comprehensive report
└─────────────────────┘
```

## 🎓 Methods & Citations

### Core Framework
**Seurat v5**
- Hao et al. *Cell* 2021
- Stuart et al. *Cell* 2019

### Quality Control
**Doublet Detection**
- scDblFinder: Germain et al. *F1000Research* 2020
- DoubletFinder: McGinnis et al. *Cell Systems* 2019

### Normalization
**SCTransform**
- Hafemeister & Satija *Genome Biology* 2019

**Integration**
- RPCA/CCA: Hao et al. *Cell* 2021
- Harmony: Korsunsky et al. *Nature Methods* 2019

### Clustering
**Leiden Algorithm**
- Traag et al. *Scientific Reports* 2019

### Cell Type Annotation
**SingleR**
- Aran et al. *Nature Immunology* 2019

**Azimuth**
- Hao et al. *Cell* 2021

### Differential Expression
**Pseudobulk DESeq2 (RECOMMENDED)**
- Squair et al. *Nature Communications* 2021
- Soneson & Robinson *Nature Methods* 2018

**MAST**
- Finak et al. *Genome Biology* 2015

### Functional Enrichment
**clusterProfiler**
- Yu et al. *OMICS* 2012
- Wu et al. *Innovation* 2021

**GSEA**
- Subramanian et al. *PNAS* 2005
- Korotkevich et al. *bioRxiv* 2021 (fgsea)

### Cell-Cell Communication
**CellChat**
- Jin et al. *Nature Communications* 2021

## 📊 Output Files Guide

### Intermediate Checkpoints
All `.rds` files are Seurat objects saved at each step for resuming analysis.

| File | Description | Use Case |
|------|-------------|----------|
| `01_raw_data.rds` | Initial loaded data | Resume from data loading |
| `02_qc_filtered.rds` | Post-QC data | Adjust normalization params |
| `03_normalized.rds` | Normalized data | Adjust clustering params |
| `04_dimred.rds` | With PCA/UMAP | Adjust clustering resolution |
| `05_clustered.rds` | Clustered data | Re-annotate cell types |
| `06_annotated.rds` | Annotated data | Re-run DEG with different parameters |
| `final_seurat_object.rds` | Complete analysis | Publication, further analysis |

### Data Tables
| File | Description | Columns |
|------|-------------|---------|
| `06_marker_genes.csv` | Cluster markers | cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2 |
| `07_DEGs_wilcox.csv` | Wilcoxon DEGs | gene, celltype, avg_log2FC, p_val, p_val_adj, pct.1, pct.2 |
| `07_DEGs_MAST.csv` | MAST DEGs | Same as above |
| `07_DEGs_DESeq2_pseudobulk.csv` | Pseudobulk DEGs (BEST) | Same as above |
| `08_GO_enrichment.csv` | GO terms | celltype, ontology, Description, GeneRatio, p.adjust |
| `08_KEGG_enrichment.csv` | KEGG pathways | celltype, Description, GeneRatio, p.adjust |
| `08_GSEA_results.csv` | GSEA gene sets | pathway, NES, pval, padj, leadingEdge |

### Plots (PDF format)
All publication-quality at 300 DPI.

| File | Contents |
|------|----------|
| `02_QC_plots.pdf` | Pre/post QC histograms, scatter plots, doublet detection |
| `03_normalization_plots.pdf` | Variable features, variance plots |
| `04_dimred_plots.pdf` | Elbow plot, PCA, UMAP, tSNE, QC metrics |
| `05_clustering_plots.pdf` | UMAP clusters, dendrogram, QC by cluster |
| `06_annotation_plots.pdf` | Cell type UMAP, marker dotplots, heatmaps |
| `07_DEG_plots.pdf` | Volcano plots, heatmaps, feature plots |
| `08_enrichment_plots.pdf` | GO/KEGG bar plots, dot plots, heatmaps |
| `09_cellchat_plots.pdf` | Circle plots, network plots, pathway activity |

## ⚙️ Configuration Examples

### Minimal Configuration (10X data)
```yaml
project:
  name: "My_Analysis"
  output_dir: "./results"

input:
  data_type: "10X"
  data_path: "/path/to/cellranger/output"

species:
  organism: "mouse"
```

### Full-Featured Multi-Sample Analysis
```yaml
project:
  name: "Control_vs_Treatment"
  output_dir: "./results"

input:
  data_type: "10X"
  samples:
    - {name: "Ctrl1", path: "/path/ctrl1", condition: "Control"}
    - {name: "Ctrl2", path: "/path/ctrl2", condition: "Control"}
    - {name: "Trt1", path: "/path/trt1", condition: "Treatment"}
    - {name: "Trt2", path: "/path/trt2", condition: "Treatment"}

species:
  organism: "mouse"

normalization:
  method: "SCT"

integration:
  perform: TRUE
  method: "RPCA"

differential_expression:
  methods: ["wilcox", "MAST", "DESeq2_pseudobulk"]
  group_by: "condition"
  ident_1: "Treatment"
  ident_2: "Control"
  pseudobulk:
    use: TRUE
    aggregate_by: "orig.ident"

cellchat:
  compare_conditions: TRUE
  conditions: ["Control", "Treatment"]
```

## 🚀 Quick Start Commands

### 1. Basic Run
```bash
Rscript main_pipeline.R config/config.yaml
```

### 2. Resume from Step
Edit config.yaml:
```yaml
advanced:
  continue_from_step: 5  # Resume from clustering
```

### 3. Custom Config
```bash
Rscript main_pipeline.R my_custom_config.yaml
```

## 💡 Best Practices

### For Publication-Quality Results

1. **Use Pseudobulk DEG**
   ```yaml
   differential_expression:
     methods: ["DESeq2_pseudobulk"]
     pseudobulk:
       use: TRUE
   ```

2. **Enable SCTransform**
   ```yaml
   normalization:
     method: "SCT"
   ```

3. **Test Multiple Resolutions**
   ```yaml
   clustering:
     test_resolutions: TRUE
     resolutions: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
   ```

4. **Use Both Manual and Automated Annotation**
   ```yaml
   annotation:
     methods:
       manual: TRUE
       SingleR: TRUE
   ```

5. **Comprehensive Functional Analysis**
   ```yaml
   functional_analysis:
     GO: {run: TRUE}
     KEGG: {run: TRUE}
     GSEA: {run: TRUE}
     Reactome: {run: TRUE}
   ```

### For Large Datasets (>100K cells)

1. **Increase computational resources**
   ```yaml
   computation:
     n_cores: 16
     memory_limit: 64
   ```

2. **Use downsampling for visualization**
   ```r
   # Add to config or use helper functions
   seurat_downsampled <- downsample_seurat(seurat_obj, max_cells_per_cluster = 1000)
   ```

3. **Use Harmony instead of RPCA for faster integration**
   ```yaml
   integration:
     method: "harmony"
   ```

## 🔧 Customization

### Adding Custom Cell Type Markers
Edit `config.yaml`:
```yaml
annotation:
  manual_markers:
    MyNewCellType:
      - "Marker1"
      - "Marker2"
      - "Marker3"
```

### Adding Custom Gene Sets for GSEA
```yaml
functional_analysis:
  custom_genesets:
    use: TRUE
    gmt_file: "/path/to/my_genesets.gmt"
```

### Changing Color Palettes
Edit `utils/plotting_functions.R`:
```r
get_color_palette <- function(n, palette = "your_favorite_palette") {
  # Customize here
}
```

## 📝 Logging and Debugging

### Pipeline Log
Check `pipeline_log.txt` for detailed execution information.

### Debug Mode
```r
# In main_pipeline.R
options(error = recover)  # Add this line for interactive debugging
```

### Check Intermediate Results
```r
# Load any checkpoint
seurat_obj <- readRDS("05_clustered.rds")
DimPlot(seurat_obj)
```

## 🎯 Typical Analysis Timeline

| Dataset Size | Expected Runtime | Memory Required |
|--------------|------------------|-----------------|
| 5K cells | 30-60 minutes | 8 GB |
| 20K cells | 1-2 hours | 16 GB |
| 50K cells | 2-4 hours | 32 GB |
| 100K+ cells | 4-8 hours | 64+ GB |

*Times are approximate and depend on:*
- Number of samples
- Integration method
- Number of DEG comparisons
- Enabled analysis modules


- **Questions**: Read the README.md and this overview
- **Custom modifications**: Edit the modular scripts in `scripts/`

---
