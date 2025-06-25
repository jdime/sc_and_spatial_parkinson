
# Single cell RNA-seq analysis of Parkinson's data

- This repo contains code to integrate scRNA-seq samples from multiple studies using `Scanpy` and other Python modules. Code is stored at `src/python`.

- The following steps are conducted:
    1. Data management
        - Load data from original sources (GEO and CELLxGENE)
        - Standardize sample metadata
        - Downsample studies to balance number of cells per study
        - Save checkpoints
    2. Integration
        - Load standardized files
        - QC filtering
        - Normalization
        - Scaling
        - Dimension reduction
        - Unsupervised clustering
        - Batch correction
        - Save checkpoints
    3. Cell type labeling
        - Identify marker genes
        - Contrast markers vs. clusters
        - Manual refinement of clusters

- Requirements specified at: `requirements.txt`
- A copy if inputs can be downloaded from: https://zenodo.org/uploads/15734764 (v2 upload)

## Installation
```bash
  git clone https://github.com/jdime/sc_and_spatial_parkinson.git
```
    
## Authors

- [@Javier Diaz-Mejia](https://www.github.com/jdime)

