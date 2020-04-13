# PhyDOSE R Package

## Installation guide

1. First, install and load the devtools package
```
install.packages("devtools")
libarary(devtools)
```
2. Run the following
```
install_github("elkebir-group/PhyDOSE", subdir="pkg")
```

## How to run PhyDOSE R package

Example execution:
```
PhyDOSE("./data/example.txt", conf_level=0.95, fn_rate=0, f_mult=1)
```
The first argument is the input filename of trees and their respective frequency matrices. Refer to ./data/example.txt for formatting.
The other arguments are optional.
