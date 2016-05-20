### Bringing Elton and Grinnell together: a quantitative framework to represent the biogeography of ecological interaction networks

**Analysis built under R 3.1.2 (2014-10-31)**

## External libraries

We installed several external libraries needed to built the R packages ```ncdf4```, ```rgdal``` with Homebrew (Multiplateform packages manager).

```code
brew install gdal
brew install netcdf
```

## R Packages dependencies

- doParallel
- dismo 
- raster
- igraph
- RColorBrewer
- sp
- rgdal
- viridis
- ncdf4

#### Install packages in one step

```r
ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
}

# Apply ipak function on packages
packages <-c('doParallel','dismo','raster','igraph','RColorBrewer','sp','rgdal','viridis','ncdf4')
ipak(packages)
```
## Getting started

1. Set ```TRUE``` or ```FALSE``` the ```parallel``` option in ```source_all.r```  
2. Run ``` Rscript source_all.r ``` in your terminal

