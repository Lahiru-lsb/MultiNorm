# Example multislice SpatialExperiment dataset

A toy dataset derived from the DLPFC dataset provided by the spatialLIBD
package. The data have been subset and simplified for demonstration
purposes.

## Usage

``` r
spe_list
```

## Format

A list of SpatialExperiment objects

## Source

Derived from the DLPFC dataset in the spatialLIBD package

## Details

The dataset consists of a list of SpatialExperiment objects, each
representing a tissue slice. Spot-level niche annotations are stored in
the `niche` column of `colData`.

## References

Pardo B, Spangler A, Weber LM, Hicks SC, Jaffe AE, Martinowich K,
Maynard KR, Collado-Torres L (2022). spatialLIBD: an R/Bioconductor
package to visualize spatially-resolved transcriptomics data. *BMC
Genomics*.
[doi:10.1186/s12864-022-08601-w](https://doi.org/10.1186/s12864-022-08601-w)
