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
License

MIT

GitHub recommends having a README because it tells others what the project does and how to use it. :contentReference[oaicite:5]{index=5}

---

## 4. Add a license

If you want others to modify it, add a license.

Good choices:

- **MIT** = very permissive
- **GPL-3.0** = modifications must stay open
- **Apache-2.0** = also permissive with patent language

If you want:
> “others can use it and modify if necessary”

then **MIT** is usually the simplest. GitHub supports adding a license file directly. :contentReference[oaicite:6]{index=6}

---

## 5. Create a new GitHub repository

On GitHub:
- click **New repository**
- name it something like `memory-dynamics-app`
- make it public
- create it empty, or initialize it with README if you want

GitHub’s docs cover creating a new repository and adding README, `.gitignore`, and license. :contentReference[oaicite:7]{index=7}

---

## 6. Upload the files

### Easiest way
Use GitHub web upload:
- open repo
- **Add file**
- **Upload files**
- drag in your files
- commit

### Better way with Git
From terminal inside your project folder:

```bash
git init
git add .
git commit -m "Initial commit of memory dynamics Shiny app"
git branch -M main
git remote add origin https://github.com/YOUR-USERNAME/memory-dynamics-app.git
git push -u origin main

```
