# Extended Joint Models for Longitudinal and Time-to-Event Data: A Tutorial

R code accompanying the tutorial *“Extended Joint Models for Longitudinal and Time-to-Event Data: A Tutorial.”* The repository contains the code used to generate the simulated datasets and reproduce all code listings presented in the tutorial.

## Repository structure

The code is organized into three R scripts:

- `functions.R` contains the functions used for data simulation.
- `data_simulation.R` generates the simulated longitudinal and event-time datasets used throughout the tutorial and saves them in the `Data/` directory.
- `tutorial_analyses.R` contains the analyses and code listings presented in the tutorial.

## How to use the files

To reproduce the examples, run `tutorial_analyses.R`. This script automatically sources `data_simulation.R`, which in turn sources `functions.R`, to generate the datasets and reproduce the models and analyses presented in the tutorial.

The script requires R and the packages used throughout the tutorial, including **JMbayes2**. The current CRAN version of **JMbayes2** can be installed with:

```r
install.packages("JMbayes2")
```
