# MDS_subtypes_RNAseq

Deconvolution of bulk RNA-seq data (from Hershberger et al.)
using **CIBERSORTx** and comparison of
cell-type fractions across molecular subtypes of MDS.

The workflow prepares a single-cell bone marrow reference (from whatever et al.),
normalizes bulk RNA-seq
data, performs deconvolution with CIBERSORTx, and analyzes inferred
cell fractions and cell-type-specific expression profiles.

> ⚠️ **Note:** CIBERSORTx must be run manually through the web interface.

---

# Workflow

The analysis consists of the following steps:

1. Subsample the reference single-cell dataset  
2. Normalize the single-cell reference and bulk RNA-seq data  
3. Run CIBERSORTx (manual step)  
4. Analyze inferred cell fractions  
5. Perform differential expression and pathway enrichment

---

# 1. Prepare input data

Run the preprocessing scripts:

```bash
Rscript R/01_subsample_reference_sc.R
Rscript R/02_normalize_and_subset.R
Deconvolution of a bulk dataset (citation) using cibersortx, and comparison 
of the cell fractions across defined subtypes
```

These scripts generate the files required for CIBERSORTx:

```
data/intermediate/
  sc_reference_subsetted.txt
  bulk_mixture.txt
  gene_subset.txt
```

# 2. Upload files to CIBERSORTx

Go to 

https://cibersortx.stanford.edu/

Navigate to 

Menu -> Upload Files 

 `data/intermediate/sc_reference_subsetted.txt` 
 `data/intermediate/bulk_mixture.txt` 
 `data/intermediate/gene_subset.txt`
 
# 3. Create signature matrix

Menu → Run CIBERSORTx
 
Parameters: 
```
Analysis Module: Create Signature Matrix
Analysis Mode: Custom
Input Data Type: scRNA-Seq

Single cell reference matrix file:
    sc_reference_subsetted.txt

Custom signature file name:
    signature_matrix.txt

Other parameters:
    default
```

# 4. Impute cell-type expression

Run CIBERSORTx again with the following settings:

```Analysis Module: Impute Cell Expression
Analysis Mode: Custom
Expression Analysis Type: High-Resolution

Signature matrix file:
    signature_matrix.txt

Mixture file:
    bulk_mixture.txt

Gene subset file:
    gene_subset.txt

Batch correction:
    enabled

Batch correction mode:
    S-mode

Single cell reference matrix file:
    sc_reference_subsetted.txt
```

# 5. Download results





