
## BEDASSLE 2.0-a1 ReadMe

**NOTE THIS PACKAGE IS BEING ACTIVELY DEVELOPED AND SHOULD ONLY BE USED FOR TESTING PUPOSES**

This repo contains the code for the method **BEDASSLE** - a statistical tool 
for modeling pairwise genetic distance as a function of geographic and ecological 
or environmental distance, for estimating the relative contributions of these 
distances to genetic differentiation, and for statistically comparing models with 
different distance predictors (e.g., geographic distance alone vs. geographic AND 
ecological distance). BEDASSLE stands for: 
**B**ayesian **E**stimation of 
**D**ifferentiation in **A**lleles by **S**patial **S**tructure and **L**ocal **E**cology.

This package replaces (**and is not backwards-compatible with**), BEDASSLE 1.0. The 
paper describing the original method can be found here:

 * [paper](https://doi.org/10.1111/evo.12193)

## Installation

Upon installation, the different **BEDASSLE** models will be compiled, which may 
spit lots of text, and possibly some warnings, to your screen. This is 
totally normal, and you should only be concerned if you get errors 
and the installation fails.

To install from github:

```r
	library(devtools)
	install_github("gbradburd/bedassle")
```

Note that Windows users may have to download Rtools as a 
standalone executable before trying to install the **BEDASSLE** R package.


## Contact

After referring to the manual, 
please direct all queries to bradburd (at) msu.edu, 
or post as issues on the git repo.