

## Loading the Package




## Loading the DataSets

The datasets used here are only for demonstration, they do not constitute any meaningful result.


```r
sample_names <- c("19-092_RLL","22-221_LLL","22-295_LUL","22-313_LUL","22-401_LUL","22-221_RLL")
ind_sample_name <- 1
sample_name <- sample_names[ind_sample_name]
print(ind_sample_name)
```

```
## [1] 1
```

```r
data.dir = paste0("/work/tm389/Visium_Lung/",sample_name) ## Change to Your Location
Lung <- Load10X_Spatial(data.dir,
                        filename = "raw_feature_bc_matrix.h5",
                        assay = "Spatial",
                        slice = "slice",
                        filter.matrix = TRUE,
                        to.upper = FALSE,
                        image = NULL)
```

```
## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
```

```r
## scRNA Data
dat_scRNA <- pbmc_small

## Visium Data
dat_visium <- Lung
```


## Transform the Datasets


```r
dat_scRNA <- InterSpatial::transform_scRNA_data(dat_scRNA)
dat_visium <- InterSpatial::transform_spatial_data(dat_visium)
```
