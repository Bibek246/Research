# Memory Dynamics Shiny App

Interactive Shiny app for exploring a two-state microbial memory model based on current Aim 3 modeling work.

## Features
- No stress, periodic stress, continuous Hill stress
- Deterministic and stochastic simulations
- Compare runs
- Memory fraction, population trajectories, phase portrait

## Run locally

```r
install.packages(c("shiny", "bslib", "ggplot2", "dplyr", "tidyr", "scales", "htmltools"))
shiny::runApp("app.R")
```
