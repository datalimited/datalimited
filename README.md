## Stock-assessment methods for data-limited fisheries

[![Build Status](https://magnum.travis-ci.com/datalimited/datalimited.svg?token=QExyQi6ySw3SZD4gggYN&branch=master)](https://magnum.travis-ci.com/datalimited/datalimited)

The R package datalimited can be installed from GitHub with:

```R
# install.packages("devtools")
devtools::install_github("hadley/devtools") # temporarily needed due to a bug 
devtools::install_github("datalimited/datalimited") 
```

Because the package includes C++ code, you will need a C++ compiler to install the package from source code. RStudio has a [good article](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites) on setting this up.

The package implements the methods used in Rosenberg et al. (2014) including the following:

- Catch-MSY based on Martell and Froese (2013); see `?cmsy`
- Catch-only method with sample importance resampling based on asconcellos and Cochrane (2005); see `?comsir`
- Panel regression in the style of Costello et al. (2012); see `?prm`
- State-space catch-only model; see `?sscom`

```R
library("datalimited")
help(package = "datalimited")
```

### References

Costello, C., D. Ovando, R. Hilborn, S. D. Gaines, O. Deschenes, and S. E. Lester. 2012. Status and Solutions for the World’s Unassessed Fisheries. Science 338:517-520.

Martell, S., and R. Froese. 2013. A simple method for estimating MSY from catch and resilience. Fish and Fisheries 14:504-514.

Rosenberg, A. A., M. J. Fogarty, A. B. Cooper, M. Dickey-Collas, E. A. Fulton, N. L. Gutiérrez, K. J. W. Hyde, K. M. Kleisner, C. Longo, C. V. Minte-Vera, C. Minto, I. Mosqueira, G. C. Osio, D. Ovando, E. R. Selig, J. T. Thorson, and Y. Ye. 2014. Developing new approaches to global stock status assessment and fishery production potential of the seas. FAO Fisheries and Aquaculture Circular, Rome, Italy.

Vasconcellos, M., and K. Cochrane. 2005. Overview of World Status of Data-Limited Fisheries: Inferences from Landings Statistics. Pages 1-20 in G. H. Kruse, V. F. Gallucci, D. E. Hay, R. I. Perry, R. M. Peterman, T. C. Shirley, P. D. Spencer, B. Wilson, and D. Woodby, editors. Fisheries Assessment and Management in Data-Limited Situations. Alaska Sea Grant, University of Alaska Fairbanks.
