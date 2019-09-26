---
title: "Bayesian Skyline Plot tests"
author: "Louis du Plessis"
date: 'Last modified: 26 Sep 2019'
output:
  html_document:
    toc: yes
    toc_float: yes
    theme: paper
    self_contained: no
    lib_dir: lib
    keep_md: yes
  md_document:
    toc: yes
    variant: markdown_github
layout: page
editor_options: 
  chunk_output_type: inline
---

# Summary

Tests to check that the new Bayesian Skyline Plot implementation (`bsp.distributions.BSP`) returns the same results as the old Bayesian Skyline Plot implementation in the BEAST2 core (`beast.evolution.tree.coalescent.BayesianSkyline`).
    





---
  
# Homochronous dataset (HCV)
  
## Likelihood comparison

Run two separate analyses with the same data, operators and seeds, using the old and new Bayesian skyline models. Check that the maximum absolute difference between the likelihoods of the old and new models are within the tolerance.


```r
  lf1 <- readLog("../../examples/BSP/output/hcv_bspold_mcpgamma_127.log",  burnin=0)
  lf2 <- readLog("../../examples/BSP/output/hcv_bsp_mcpgamma_127.log", burnin=0)
  maxdiff <- max(abs(lf1[, "BayesianSkyline"] - lf2[, "BayesianSkyline"]))
```

Maximum likelihood difference: 1.955152e-08

Pass test: **TRUE**


## Direct likelihood comparison

Run one analysis, with both the old and new models as treepriors, sharing the same parameters. Check that the maximum absolute difference between the likelihoods of the old and new models are within the tolerance.


```r
  lf <- readLog("../../examples/BSP/output/hcv_bsp_mcpgamma_skylinecomparison_127.log", burnin=0)
  maxdiff <- max(abs(lf[ ,"BayesianSkylineNew"] - lf[, "BayesianSkyline"]))
```

Maximum likelihood difference: 5.684342e-14

Pass test: **TRUE**




## Change time logger comparison

Check that the change times logged by `BayesianSkylineChangeTimesLogger` are identical to the change times calculated from the logged trees. 


```r
    # Change point times logged by change time logger
    lf <- readLog("../../examples/BSP/output/hcv_bsp_mcpgamma_127.log", burnin=0)
    changeTimes <- getLogFileSubset(lf, "bChangeTimes")
    
    # Calculate change point times from the logged trees and logged group sizes
    tree <- readTreeLog("../../examples/BSP/output/hcv_bsp_mcpgamma_127.trees", burnin=0)
    groupSizes <- getLogFileSubset(lf,"bGroupSizes")
    treeTimes  <- getChangeTimes(tree, groupSizes)
    
    maxdiff <- max(abs(changeTimes - treeTimes))
```

Maximum change time difference: 2.700062e-13

Pass test: **TRUE**




---

# Heterochronous dataset (bison)
  
## Likelihood comparison

Run two separate analyses with the same data, operators and seeds, using the old and new Bayesian skyline models. Check that the maximum absolute difference between the likelihoods of the old and new models are within the tolerance.


```r
  lf1 <- readLog("../../examples/BSP/output/bison_bspold_mcpgamma_127.log",  burnin=0)
  lf2 <- readLog("../../examples/BSP/output/bison_bsp_mcpgamma_127.log", burnin=0)
  maxdiff <- max(abs(lf1[, "BayesianSkyline"] - lf2[, "BayesianSkyline"]))
```

Maximum likelihood difference: 1.336753e-08

Pass test: **TRUE**


## Direct likelihood comparison

Run one analysis, with both the old and new models as treepriors, sharing the same parameters. Check that the maximum absolute difference between the likelihoods of the old and new models are within the tolerance.


```r
  lf <- readLog("../../examples/BSP/output/bison_bsp_mcpgamma_skylinecomparison_127.log", burnin=0)
  maxdiff <- max(abs(lf[ ,"BayesianSkylineNew"] - lf[, "BayesianSkyline"]))
```

Maximum likelihood difference: 4.547474e-13

Pass test: **TRUE**



## Change time logger comparison

Check that the change times logged by `BayesianSkylineChangeTimesLogger` are identical to the change times calculated from the logged trees. 




```r
    # Change point times logged by change time logger
    lf <- readLog("../../examples/BSP/output/bison_bsp_mcpgamma_127.log", burnin=0)
    changeTimes <- getLogFileSubset(lf, "bPopSizeChangeTimes")
    
    # Calculate change point times from the logged trees and logged group sizes
    tree <- readTreeLog("../../examples/BSP/output/bison_bsp_mcpgamma_127.trees", burnin=0)
    groupSizes <- getLogFileSubset(lf,"bGroupSizes")
    treeTimes  <- getChangeTimes(tree, groupSizes)
    
    maxdiff <- max(abs(changeTimes - treeTimes))
```

Maximum change time difference: 1.009539e-10

Pass test: **TRUE**


---

# Session info


```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS Sierra 10.12.6
## 
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] rskylinetools_0.1.1 beastio_0.2.2      
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.2      ape_5.3         lattice_0.20-38 digest_0.6.21  
##  [5] grid_3.5.1      nlme_3.1-139    magrittr_1.5    evaluate_0.13  
##  [9] coda_0.19-3     stringi_1.4.3   rmarkdown_1.12  tools_3.5.1    
## [13] stringr_1.4.0   parallel_3.5.1  xfun_0.6        yaml_2.2.0     
## [17] compiler_3.5.1  htmltools_0.3.6 knitr_1.22
```

