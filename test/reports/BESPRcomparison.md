---
title: "Bayesian Epoch Sampling Skyline Plot R-comparisonn"
author: "Louis du Plessis"
date: 'Last modified: 09 Oct 2019'
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

Tests to check that the Bayesian Epoch Sampling Skyline Plot implementation (`bsp.distributions.BESP`) returns the same likelihood as the `rskylinetools` `esp_lk` function.
Keep in mind that the likelihoods computed in R need to be calculated with `full = FALSE`, which will remove the constant terms from the likelihood (as is done in the
BEAST2 implementation).

    



---

# Heterochronous dataset (bison)

## Density-defined sampling model

Check that the maximum absolute difference between the likelihoods of the skyline plot calculated in BEAST2 and R are within the tolerance.


```r
  lf    <- readLog("../../examples/BESP/output/bison_besp_mcpgamma_limits_densdefined_127.log", burnin=0)
  trees <- readTreeLog("../../examples/BESP/output/bison_besp_mcpgamma_limits_densdefined_127.trees", burnin=0)

  popSizes      <- getLogFileSubset(lf, "bPopSizes")
  groupSizes    <- getLogFileSubset(lf, "bGroupSizes")
  sampIntensity <- lf[, "bSamplingIntensity"]
  beastlk       <- lf[, "BayesianSkyline"]

  rlk <- sapply(1:length(beastlk), function(i) esp_lk(trees[[i]], popSizes = popSizes[i,], segmentSizes = groupSizes[i,], sampIntensities = sampIntensity[i], full=FALSE))
  
  maxdiff <- max(abs(beastlk - rlk))
```

Maximum likelihood difference: 1.117587e-08

Pass test: **TRUE**

---

## Frequency-defined sampling model

Check that the maximum absolute difference between the likelihoods of the skyline plot calculated in BEAST2 and R are within the tolerance.


```r
  lf    <- readLog("../../examples/BESP/output/bison_besp_mcpgamma_limits_freqdefined_127.log", burnin=0)
  trees <- readTreeLog("../../examples/BESP/output/bison_besp_mcpgamma_limits_freqdefined_127.trees", burnin=0)

  popSizes        <- getLogFileSubset(lf, "bPopSizes")
  groupSizes      <- getLogFileSubset(lf, "bGroupSizes")
  sampIntensities <- getLogFileSubset(lf, "bSamplingIntensity")
  epochTimes      <- getLogFileSubset(lf, "bSamplingIntensityChangeTimes")[1,]
  epochEvents     <- sapply(epochTimes, function(i) which(abs(i - getSamplingTimes(trees[[1]])) < tol))
  epochSizes      <- c(epochEvents[1], diff(epochEvents))
  beastlk         <- lf[, "BayesianSkyline"]

  rlk <- sapply(1:length(beastlk), function(i) esp_lk(trees[[i]], popSizes = popSizes[i,], segmentSizes = groupSizes[i,], sampIntensities = sampIntensities[i,], epochSizes = epochSizes, full=FALSE))
  
  maxdiff <- max(abs(beastlk - rlk))
```

Maximum likelihood difference: 1.117587e-08

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
## [1] rskylinetools_0.2.1 beastio_0.2.4      
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.2      ape_5.3         lattice_0.20-38 digest_0.6.21  
##  [5] grid_3.5.1      nlme_3.1-139    magrittr_1.5    evaluate_0.13  
##  [9] coda_0.19-3     stringi_1.4.3   rmarkdown_1.12  tools_3.5.1    
## [13] stringr_1.4.0   parallel_3.5.1  xfun_0.6        yaml_2.2.0     
## [17] compiler_3.5.1  htmltools_0.3.6 knitr_1.22
```

