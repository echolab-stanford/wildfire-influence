# wildfire-influence
Repo supporting [Burke et al 2023 "The contribution of wildfire to PM<sub>2.5</sub> trends in the USA"](https://doi.org/10.1038/s41586-023-06522-6).

Results from the paper are in the `figures/clean` and `tables/clean` folders. Code to replicate results are in the `scripts` folder. Data are in [Dropbox](https://www.dropbox.com/sh/3zz7ri3uzc5uf6t/AAAcwLegWlEkA31EkDXuEPZna?dl=0).

## How to replicate results
1. Download this repository.
2. Download the [Dropbox](https://www.dropbox.com/sh/3zz7ri3uzc5uf6t/AAAcwLegWlEkA31EkDXuEPZna?dl=0) folder. Place files downloaded from Dropbox in the same folder as the downloaded GitHub repository.
3. Change settings in `scripts/setup/00_03_load_settings.R`:
    1. Set `path_dropbox` to the location of the data downloaded from Dropbox.
    2. Set `path_github` to the location of this downloaded repository's root.
4. Install packages by running `scripts/setup/00_00_install_packages.R`.
5. Set working directory to this downloaded repository's root.
6. Run scripts in `scripts/main`.

## Computational environment
```
R version 4.3.0 (2023-04-21)
Platform: x86_64-apple-darwin20 (64-bit)
Running under: macOS Ventura 13.3.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

GEOS 3.11.2
GDAL 3.6.4
PROJ 9.2.0
```

R packages versions are specified in `scripts/setup/00_00_install_packages.R` and below.
```
attached base packages:
[1] parallel  stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
 [1] callr_3.7.0         tictoc_1.2          doFuture_1.0.0      doRNG_1.8.6         rngtools_1.5.2     
 [6] future.callr_0.8.1  future.apply_1.11.0 future_1.32.0       doParallel_1.0.17   iterators_1.0.14   
[11] foreach_1.5.2       scales_1.2.1        usmap_0.6.1         MetBrewer_0.2.0     ggridges_0.5.3     
[16] ggforce_0.4.1       geofacet_0.2.0      cowplot_1.1.1       gtable_0.3.1        strucchange_1.5-3  
[21] sandwich_3.0-1      zoo_1.8-9           robslopes_1.1.3     fixest_0.11.1       feather_0.3.5      
[26] forcats_1.0.0       stringr_1.5.0       dplyr_1.1.2         purrr_1.0.1         readr_2.1.4        
[31] tidyr_1.3.0         tibble_3.2.1        ggplot2_3.4.2       tidyverse_2.0.0     magrittr_2.0.3     
[36] tigris_2.0.3        sf_0.9-8            lubridate_1.9.2    

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.0    farver_2.1.1        tweenr_2.0.2        digest_0.6.27       timechange_0.2.0   
 [6] lifecycle_1.0.3     dreamerr_1.2.3      processx_3.5.2      compiler_4.3.0      rlang_1.1.1        
[11] tools_4.3.0         utf8_1.2.2          sp_1.6-0            classInt_0.4-8      plyr_1.8.6         
[16] imguR_1.0.3         geogrid_0.1.1       KernSmooth_2.23-20  withr_2.5.0         numDeriv_2016.8-1.1
[21] grid_4.3.0          polyclip_1.10-4     fansi_1.0.4         e1071_1.7-12        colorspace_2.1-0   
[26] globals_0.16.2      MASS_7.3-58.2       cli_3.6.0           generics_0.1.0      rstudioapi_0.14    
[31] httr_1.4.6          tzdb_0.1.2          DBI_1.1.3           proxy_0.4-27        splines_4.3.0      
[36] rnaturalearth_0.3.2 vctrs_0.6.2         jsonlite_1.8.4      hms_1.1.3           ggrepel_0.9.3      
[41] Formula_1.2-4       listenv_0.9.0       jpeg_0.1-10         import_1.3.0        units_0.8-2        
[46] parallelly_1.35.0   glue_1.6.2          ps_1.6.0            codetools_0.2-19    stringi_1.7.5      
[51] munsell_0.5.0       pillar_1.9.0        rappdirs_0.3.3      R6_2.5.1            lattice_0.20-45    
[56] png_0.1-8           segmented_1.6-4     renv_0.14.0         class_7.3-21        Rcpp_1.0.10        
[61] uuid_0.1-4          gridExtra_2.3       nlme_3.1-161        rgeos_0.6-1         pkgconfig_2.0.3    
```
