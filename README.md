# Memory Dynamics Shiny App

An interactive **R Shiny** app for exploring a two-state microbial memory model based on current Aim 3 modeling work.

This app is designed to help visualize how **unexposed cells** and **memory cells** change over time under different stress conditions, and how parameters such as switching, decay, growth cost, and noise affect the dynamics.

## Features

- Interactive simulation of a **two-state memory model**
- Stress scenarios:
  - No stress
  - Periodic stress
  - Continuous Hill-type stress
- Deterministic and stochastic simulation modes
- Adjustable biological and simulation parameters
- Population trajectories
- Memory fraction over time
- Stress forcing visualization
- Phase portrait
- Side-by-side comparison mode
- Preset scenarios based on the current research Rmd

## Model Overview

The app tracks two cell populations:

- **U** = unexposed cells
- **M** = memory cells

The model includes:

- Logistic growth
- Stress-induced switching from **U** to **M**
- Baseline symmetric switching between states
- Memory decay from **M** back to **U**
- Reduced growth of memory cells
- Optional stochastic noise

This app is an interactive implementation of a developing research model and may change as the underlying Rmd/model evolves.

## Main Parameters

- `r0` — baseline growth rate
- `c_cost` — growth multiplier for memory cells
- `alpha` — stress penalty on unexposed-cell growth
- `mu` — baseline switching rate
- `eps` — memory decay rate
- `beta` — stress-driven switching rate
- `sigma` — stochastic noise strength

## Run Locally

Install the required packages in R:

```r
install.packages(c(
  "shiny",
  "bslib",
  "ggplot2",
  "dplyr",
  "tidyr",
  "scales",
  "htmltools"
))

Then run the app:

shiny::runApp("app.R")
```
Deployment

This app can be deployed to a Shiny hosting service such as shinyapps.io.

Typical deployment workflow:

Keep the code in a public GitHub repository
Deploy the live app from the same project folder
Share:
GitHub repository link for code access
Shiny deployment link for browser access
Repository Goals

This repository is intended to allow others to:

View the code
Run the app locally
Modify and extend the model
Use the app as a starting point for future Aim 3 modeling work
Status

This is an active research development project.
The current version reflects the present modeling direction and interface, but the equations, assumptions, and implementation may continue to evolve.

License

This project is released under the MIT License.
See the LICENSE file for details.
Use this link to access the website:

https://research12.shinyapps.io/research/
