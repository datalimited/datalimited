## Stock assessment methods for data-limited fisheries

The R package datalimited can be installed from GitHub with:

```R
# install.packages("devtools")
devtools::install_github("datalimited/datalimited", 
  auth_token = "625bc066956a5f285e4fe45f4084bb6a26afb9ed")
```

The `auth_token` argument is only required as long as this repository is private.

Because the package includes C++ code, you will need a C++ compiler to install from the source code. RStudio has a [good article](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites) on setting this up. We will also release compiled versions of the package via GitHub and CRAN.

The package currently implements the following methods:

- Catch-only method with sample importance resampling (Vasconcellos and Cochrane 2005) via the function `comsir()`
- Catch-MSY (Martell and Froese 2013) via the function `cmsy()`

### References

Costello, C., D. Ovando, R. Hilborn, S. D. Gaines, O. Deschenes, and S. E. Lester. 2012. Status and Solutions for the World’s Unassessed Fisheries. Science 338:517-520.

Martell, S., and R. Froese. 2013. A simple method for estimating MSY from catch and resilience. Fish and Fisheries 14:504-514.

Rosenberg, A. A., M. J. Fogarty, A. B. Cooper, M. Dickey-Collas, E. A. Fulton, N. L. Gutiérrez, K. J. W. Hyde, K. M. Kleisner, C. Longo, C. V. Minte-Vera, C. Minto, I. Mosqueira, G. C. Osio, D. Ovando, E. R. Selig, J. T. Thorson, and Y. Ye. 2014. Developing new approaches to global stock status assessment and fishery production potential of the seas. FAO Fisheries and Aquaculture Circular, Rome, Italy.

Vasconcellos, M., and K. Cochrane. 2005. Overview of World Status of Data-Limited Fisheries: Inferences from Landings Statistics. Pages 1-20 in G. H. Kruse, V. F. Gallucci, D. E. Hay, R. I. Perry, R. M. Peterman, T. C. Shirley, P. D. Spencer, B. Wilson, and D. Woodby, editors. Fisheries Assessment and Management in Data-Limited Situations. Alaska Sea Grant, University of Alaska Fairbanks.
