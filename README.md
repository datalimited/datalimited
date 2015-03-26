## Stock assessment methods for data-limited fisheries

The R package datalimited can be installed with:

```R
# install.packages("devtools")
devtools::install_github("data-limited-wg/datalimited", 
  auth_token = "625bc066956a5f285e4fe45f4084bb6a26afb9ed")
```

The package implements the following methods:

### Catch-only method with sample importance resampling: `comsir()`

Vasconcellos and Cochrane (2005)

### Catch-MSY: `cmsy()`

Martell and Froese (2013)

### References

Costello, C., D. Ovando, R. Hilborn, S. D. Gaines, O. Deschenes, and S. E. Lester. 2012. Status and Solutions for the World’s Unassessed Fisheries. Science 338:517-520.

Martell, S., and R. Froese. 2013. A simple method for estimating MSY from catch and resilience. Fish and Fisheries 14:504-514.

Rosenberg, A. A., M. J. Fogarty, A. B. Cooper, M. Dickey-Collas, E. A. Fulton, N. L. Gutiérrez, K. J. W. Hyde, K. M. Kleisner, C. Longo, C. V. Minte-Vera, C. Minto, I. Mosqueira, G. C. Osio, D. Ovando, E. R. Selig, J. T. Thorson, and Y. Ye. 2014. Developing new approaches to global stock status assessment and fishery production potential of the seas. FAO Fisheries and Aquaculture Circular, Rome, Italy.

Vasconcellos, M., and K. Cochrane. 2005. Overview of World Status of Data-Limited Fisheries: Inferences from Landings Statistics. Pages 1-20 in G. H. Kruse, V. F. Gallucci, D. E. Hay, R. I. Perry, R. M. Peterman, T. C. Shirley, P. D. Spencer, B. Wilson, and D. Woodby, editors. Fisheries Assessment and Management in Data-Limited Situations. Alaska Sea Grant, University of Alaska Fairbanks.
