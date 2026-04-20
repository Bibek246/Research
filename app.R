library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(htmltools)

clip01 <- function(x) pmax(0, pmin(1, x))
`%||%` <- function(x, y) if (is.null(x)) y else x

get_display_labels <- function(label_U = "U", label_M = "M") {
  c(
    U = ifelse(nchar(trimws(label_U)) == 0, "U", trimws(label_U)),
    M = ifelse(nchar(trimws(label_M)) == 0, "M", trimws(label_M)),
    total = "total"
  )
}

stress_signal <- function(times, scenario, beta, period_hours, duty_cycle,
                          s_base, s_amp, forcing_period, rho50, hill_n) {
  if (scenario == "No stress") {
    tibble(time = times, stress_raw = 0, beta_eff = 0, rho = 0, iota = 0)
  } else if (scenario == "Periodic stress") {
    iota <- (((times) %/% period_hours) %% 2)
    tibble(time = times, stress_raw = iota, beta_eff = beta * iota, rho = iota, iota = iota)
  } else {
    stress <- sin((times) / (forcing_period / (2 * pi))) + 0.5
    stress <- pmax(stress, 0)
    rho <- ifelse(stress <= 0, 0, 1 / (1 + (rho50 / stress)^hill_n))
    tibble(time = times, stress_raw = stress, beta_eff = beta * rho, rho = rho, iota = NA_real_)
  }
}

simulate_memory <- function(scenario = "No stress",
                            tmax = 120,
                            dt = 0.1,
                            U0 = 0.98,
                            M0 = 0.02,
                            r0 = 0.8,
                            c_cost = 0.8,
                            alpha = 0.4,
                            mu = 0.002,
                            eps = 0.03,
                            beta = 0.10,
                            sigma = 0.02,
                            stochastic = TRUE,
                            seed = 123,
                            period_hours = 8,
                            duty_cycle = 0.5,
                            s_base = 0.6,
                            s_amp = 0.4,
                            forcing_period = 24,
                            rho50 = 0.5,
                            hill_n = 4) {
  
  times <- seq(0, tmax, by = dt)
  n <- length(times)
  sig <- stress_signal(times, scenario, beta, period_hours, duty_cycle,
                       s_base, s_amp, forcing_period, rho50, hill_n)
  
  U <- numeric(n)
  M <- numeric(n)
  U[1] <- U0
  M[1] <- M0
  
  if (!is.null(seed) && stochastic) set.seed(seed)
  
  for (i in 2:n) {
    b <- sig$beta_eff[i - 1]
    iota <- sig$iota[i - 1]
    rho <- sig$rho[i - 1]
    Kterm <- (1 - U[i - 1] - M[i - 1])
    
    rU <- if (scenario == "Periodic stress") {
      (1 - iota * alpha) * r0
    } else if (scenario == "Continuous stress") {
      (1 - alpha * rho) * r0
    } else {
      r0
    }
    rM <- c_cost * r0
    
    dU_det <- (rU * Kterm * U[i - 1]) - (b * U[i - 1]) - mu * (U[i - 1] - M[i - 1]) + eps * M[i - 1]
    dM_det <- (rM * Kterm * M[i - 1]) + (b * U[i - 1]) - mu * (M[i - 1] - U[i - 1]) - eps * M[i - 1]
    
    if (stochastic) {
      dW_U <- rnorm(1, mean = 0, sd = sqrt(dt))
      dW_M <- rnorm(1, mean = 0, sd = sqrt(dt))
      U_new <- U[i - 1] + dU_det * dt + sigma * U[i - 1] * dW_U
      M_new <- M[i - 1] + dM_det * dt + sigma * M[i - 1] * dW_M
    } else {
      U_new <- U[i - 1] + dU_det * dt
      M_new <- M[i - 1] + dM_det * dt
    }
    
    U[i] <- max(0, U_new)
    M[i] <- max(0, M_new)
  }
  
  tibble(
    time = times,
    U = U,
    M = M,
    total = U + M,
    memory_fraction = ifelse((U + M) > 0, M / (U + M), 0),
    stress_raw = sig$stress_raw,
    beta_eff = sig$beta_eff,
    rho = sig$rho,
    iota = sig$iota
  )
}

simulate_replicates <- function(reps = 20, seed_base = 123, ...) {
  bind_rows(lapply(seq_len(reps), function(i) {
    dat <- simulate_memory(seed = seed_base + i - 1, ...)
    dat$replicate <- i
    dat
  }))
}

summarize_replicates <- function(rep_dat) {
  rep_dat %>%
    group_by(time) %>%
    summarize(
      U_med = median(U),
      U_lo = quantile(U, 0.1),
      U_hi = quantile(U, 0.9),
      M_med = median(M),
      M_lo = quantile(M, 0.1),
      M_hi = quantile(M, 0.9),
      mem_med = median(memory_fraction),
      mem_lo = quantile(memory_fraction, 0.1),
      mem_hi = quantile(memory_fraction, 0.9),
      total_med = median(total),
      total_lo = quantile(total, 0.1),
      total_hi = quantile(total, 0.9),
      stress_raw = first(stress_raw),
      beta_eff = first(beta_eff),
      .groups = "drop"
    )
}

calc_equilibrium <- function(scenario, mu, eps, beta) {
  if (scenario == "No stress") {
    mu / (eps + 2 * mu)
  } else if (scenario == "Periodic stress") {
    (beta + mu) / (beta + 2 * mu + eps)
  } else {
    NA_real_
  }
}

build_params <- function(input, suffix = "") {
  getv <- function(id) input[[paste0(id, suffix)]]
  list(
    scenario = getv("scenario"),
    tmax = getv("tmax"),
    dt = getv("dt"),
    U0 = getv("U0"),
    M0 = getv("M0"),
    r0 = getv("r0"),
    c_cost = getv("c_cost"),
    alpha = getv("alpha"),
    mu = getv("mu"),
    eps = getv("eps"),
    beta = getv("beta"),
    sigma = if (isTRUE(getv("stochastic"))) getv("sigma") else 0,
    stochastic = getv("stochastic"),
    seed = getv("seed"),
    period_hours = getv("period_hours") %||% 8,
    duty_cycle = getv("duty_cycle") %||% 0.5,
    s_base = getv("s_base") %||% 0.6,
    s_amp = getv("s_amp") %||% 0.4,
    forcing_period = getv("forcing_period") %||% 24,
    rho50 = getv("rho50") %||% 0.5,
    hill_n = getv("hill_n") %||% 4
  )
}

preset_specs <- list(
  no_stress = list(
    scenario = "No stress", stochastic = FALSE, tmax = 120, dt = 0.1,
    U0 = 0.98, M0 = 0.02, r0 = 0.8, c_cost = 0.8, alpha = 0.4, mu = 0.002, eps = 0.03,
    beta = 0.10, sigma = 0.02, period_hours = 8, duty_cycle = 0.5,
    s_base = 0.6, s_amp = 0.4, forcing_period = 24, rho50 = 0.5, hill_n = 4
  ),
  periodic_8h = list(
    scenario = "Periodic stress", stochastic = TRUE, tmax = 120, dt = 0.1,
    U0 = 0.98, M0 = 0.02, r0 = 0.8, c_cost = 0.75, alpha = 0.4, mu = 0.002, eps = 0.03,
    beta = 0.10, sigma = 0.02, period_hours = 8, duty_cycle = 0.5,
    s_base = 0.6, s_amp = 0.4, forcing_period = 24, rho50 = 0.5, hill_n = 4
  ),
  periodic_24h = list(
    scenario = "Periodic stress", stochastic = TRUE, tmax = 168, dt = 0.1,
    U0 = 0.98, M0 = 0.02, r0 = 0.8, c_cost = 0.75, alpha = 0.4, mu = 0.002, eps = 0.03,
    beta = 0.10, sigma = 0.02, period_hours = 24, duty_cycle = 0.5,
    s_base = 0.6, s_amp = 0.4, forcing_period = 24, rho50 = 0.5, hill_n = 4
  ),
  continuous_hill = list(
    scenario = "Continuous stress", stochastic = TRUE, tmax = 168, dt = 0.1,
    U0 = 0.98, M0 = 0.02, r0 = 0.8, c_cost = 0.75, alpha = 0.4, mu = 0.002, eps = 0.03,
    beta = 0.10, sigma = 0.02, period_hours = 8, duty_cycle = 0.5,
    s_base = 0.6, s_amp = 0.4, forcing_period = 24, rho50 = 0.5, hill_n = 4
  )
)

apply_preset <- function(session, spec, suffix = "") {
  upd_sel <- function(id, value) updateSelectInput(session, paste0(id, suffix), selected = value)
  upd_chk <- function(id, value) updateCheckboxInput(session, paste0(id, suffix), value = value)
  upd_sld <- function(id, value) updateSliderInput(session, paste0(id, suffix), value = value)
  upd_num <- function(id, value) updateNumericInput(session, paste0(id, suffix), value = value)
  upd_sel("scenario", spec$scenario)
  upd_chk("stochastic", spec$stochastic)
  upd_sld("tmax", spec$tmax)
  upd_sld("dt", spec$dt)
  upd_sld("U0", spec$U0)
  upd_sld("M0", spec$M0)
  upd_sld("r0", spec$r0)
  upd_sld("c_cost", spec$c_cost)
  upd_sld("alpha", spec$alpha)
  upd_sld("mu", spec$mu)
  upd_sld("eps", spec$eps)
  upd_sld("beta", spec$beta)
  upd_sld("sigma", spec$sigma)
  upd_sld("period_hours", spec$period_hours)
  upd_sld("duty_cycle", spec$duty_cycle)
  upd_sld("s_base", spec$s_base)
  upd_sld("s_amp", spec$s_amp)
  upd_sld("forcing_period", spec$forcing_period)
  upd_sld("rho50", spec$rho50)
  upd_sld("hill_n", spec$hill_n)
  upd_num("seed", 123)
}

kpi_card <- function(title, value_output, subtitle = NULL) {
  card(
    class = "border-0 shadow-sm bg-white bg-opacity-10",
    card_body(
      div(style = "font-size: .8rem; text-transform: uppercase; letter-spacing: .05em; opacity: .8;", title),
      div(textOutput(value_output), style = "font-size: 1.5rem; font-weight: 700; line-height: 1.1;"),
      if (!is.null(subtitle)) div(subtitle, style = "font-size: .78rem; opacity: .8;")
    )
  )
}

scenario_controls <- function(prefix = "") {
  tagList(
    conditionalPanel(
      sprintf("input['%sscenario'] == 'Periodic stress'", prefix),
      div(
        class = "control-block",
        h5("Periodic forcing"),
        sliderInput(paste0("period_hours", prefix), "Half-cycle length (hours)", min = 2, max = 48, value = 8, step = 1),
        helpText("Rmd logic alternates 0/1 every half-cycle, so 8 means 8h off then 8h on."),
        sliderInput(paste0("duty_cycle", prefix), "Duty cycle (kept for compatibility; ignored in Rmd mode)", min = 0.1, max = 0.9, value = 0.5, step = 0.05)
      )
    ),
    conditionalPanel(
      sprintf("input['%sscenario'] == 'Continuous stress'", prefix),
      div(
        class = "control-block",
        h5("Continuous forcing"),
        sliderInput(paste0("s_base", prefix), HTML("Baseline stress (<i>s</i><sub>0</sub>)"), min = 0, max = 2, value = 0.6, step = 0.05),
        sliderInput(paste0("s_amp", prefix), "Stress amplitude", min = 0, max = 2, value = 0.4, step = 0.05),
        sliderInput(paste0("forcing_period", prefix), "Forcing period (hours)", min = 4, max = 72, value = 24, step = 1),
        sliderInput(paste0("rho50", prefix), HTML("Hill midpoint (<i>&rho;</i><sub>50</sub>)"), min = 0.05, max = 2, value = 0.5, step = 0.05),
        sliderInput(paste0("hill_n", prefix), HTML("Hill coefficient (<i>n</i>)"), min = 1, max = 10, value = 4, step = 1)
      )
    )
  )
}

parameter_panel <- function(prefix = "", title = "Primary run", compare = FALSE) {
  suffix <- if (prefix == "") "" else paste0("_", prefix)
  tagList(
    div(
      class = "control-block",
      h5(title),
      if (compare) div(class = "preset-row",
                       actionButton(paste0("preset_no", suffix), "No stress", class = "btn btn-outline-light btn-sm"),
                       actionButton(paste0("preset_p8", suffix), "8h", class = "btn btn-outline-light btn-sm"),
                       actionButton(paste0("preset_p24", suffix), "24h", class = "btn btn-outline-light btn-sm"),
                       actionButton(paste0("preset_cont", suffix), "Hill", class = "btn btn-outline-light btn-sm")),
      selectInput(paste0("scenario", suffix), "Stress regime", c("No stress", "Periodic stress", "Continuous stress")),
      checkboxInput(paste0("stochastic", suffix), "Include stochastic noise", TRUE),
      sliderInput(paste0("tmax", suffix), "Simulation horizon (hours)", min = 24, max = 240, value = 120, step = 12),
      sliderInput(paste0("dt", suffix), "Time step", min = 0.01, max = 0.5, value = 0.10, step = 0.01),
      numericInput(paste0("seed", suffix), "Random seed", value = 123, min = 1, step = 1)
    ),
    div(
      class = "control-block",
      h5("Initial conditions"),
      sliderInput(paste0("U0", suffix), "Initial U", min = 0, max = 1, value = 0.98, step = 0.01),
      sliderInput(paste0("M0", suffix), "Initial M", min = 0, max = 1, value = 0.02, step = 0.01)
    ),
    div(
      class = "control-block",
      h5("Core parameters"),
      sliderInput(paste0("r0", suffix), HTML("Growth rate of U (<i>r</i><sub>0</sub>)"), min = 0.05, max = 2.0, value = 0.80, step = 0.01),
      sliderInput(paste0("c_cost", suffix), HTML("Relative growth of M (<i>c</i>)"), min = 0.1, max = 1.2, value = 0.80, step = 0.01),
      sliderInput(paste0("alpha", suffix), HTML("Stress penalty on U growth (<i>&alpha;</i>)"), min = 0, max = 1.0, value = 0.40, step = 0.01),
      sliderInput(paste0("mu", suffix), HTML("Symmetric baseline switching (<i>&mu;</i>)"), min = 0, max = 0.05, value = 0.002, step = 0.001),
      sliderInput(paste0("eps", suffix), HTML("Memory decay (<i>&epsilon;</i>)"), min = 0, max = 0.20, value = 0.03, step = 0.001),
      sliderInput(paste0("beta", suffix), HTML("Stress-driven switching (<i>&beta;</i>)"), min = 0, max = 0.50, value = 0.10, step = 0.005),
      conditionalPanel(sprintf("input['%sstochastic'] == true", suffix),
                       sliderInput(paste0("sigma", suffix), HTML("Noise amplitude (<i>&sigma;</i>)"), min = 0, max = 0.30, value = 0.02, step = 0.005))
    ),
    scenario_controls(prefix = suffix)
  )
}

ui <- page_navbar(
  title = div(style = "font-weight:700;", "Memory Dynamics Explorer"),
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#6d28d9",
    secondary = "#0f766e",
    bg = "#f8fafc",
    fg = "#0f172a",
    base_font = font_google("Inter"),
    heading_font = font_google("Space Grotesk")
  ),
  header = tags$head(tags$style(HTML("
    body { background: #f8fafc; }
    .app-hero{
      padding: 1rem 1.25rem;
      border-radius: 18px;
      background: linear-gradient(135deg, rgba(109,40,217,.08), rgba(14,165,233,.08));
      border: 1px solid rgba(15,23,42,.08);
      margin: .75rem .75rem 0 .75rem
    }
    .control-block{
      padding: .9rem 1rem;
      border-radius: 16px;
      background: white;
      border: 1px solid #e2e8f0;
      margin-bottom: .8rem;
      box-shadow: 0 2px 10px rgba(15,23,42,.04);
    }
    .preset-row .btn{
      margin-right: .4rem;
      margin-bottom: .4rem;
      border-radius: 999px;
    }
    .metric-card{
      background: white;
      border: 1px solid #e2e8f0;
      border-radius: 18px;
      box-shadow: 0 2px 10px rgba(15,23,42,.04);
      min-height: 130px;
    }
    .metric-title{
      font-size: .78rem;
      text-transform: uppercase;
      letter-spacing: .06em;
      color: #475569;
      margin-bottom: .35rem;
    }
    .metric-value{
      font-size: 2rem;
      font-weight: 800;
      color: #0f172a;
      line-height: 1.1;
    }
    .metric-sub{
      font-size: .82rem;
      color: #64748b;
      margin-top: .35rem;
    }
    .card-header{
      font-weight: 700;
      font-size: 1.05rem;
    }
    .small-muted{font-size:.84rem;color:#64748b;}
    .nav-link{ font-weight: 600; }
    .bslib-sidebar-layout > .main {
      padding-top: .5rem;
    }
  "))),
  
  nav_panel(
    "Simulator",
    div(
      class = "app-hero",
      h3("Prototype"),
      p("Use the tabs below to focus on one plot at a time.")
    ),
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        open = "desktop",
        div(
          class = "control-block",
          h5("Quick presets"),
          div(
            class = "preset-row",
            actionButton("preset_no", "No stress analytic", class = "btn btn-outline-primary btn-sm"),
            actionButton("preset_p8", "Periodic 8h", class = "btn btn-outline-primary btn-sm"),
            actionButton("preset_p24", "Periodic 24h", class = "btn btn-outline-primary btn-sm"),
            actionButton("preset_cont", "Continuous Hill", class = "btn btn-outline-primary btn-sm")
          )
        ),
        div(
          class = "control-block",
          h5("State names"),
          textInput("label_U", "Name for the unexposed state", value = "U"),
          textInput("label_M", "Name for the memory state", value = "M")
        ),
        div(
          class = "control-block",
          h5("Playback"),
          sliderInput("view_time", "Displayed time", min = 0, max = 120, value = 120, step = 0.5,
                      animate = animationOptions(interval = 350, loop = FALSE)),
          p(class = "small-muted", "This slider trims the plotted trajectory up to the selected time point."),
          checkboxInput("follow_end", "Auto-jump displayed time to the simulation end when parameters change", FALSE)
        ),
        div(
          class = "control-block",
          h5("Stochastic replicates"),
          checkboxInput("show_bands", "Show replicate uncertainty bands", TRUE),
          sliderInput("replicates", "Number of replicate runs", min = 5, max = 100, value = 30, step = 5),
          p(class = "small-muted", "Bands show the 10th–90th percentile envelope across stochastic runs.")
        ),
        div(
          class = "control-block",
          h5("Export"),
          downloadButton("download_data", "Simulated data (.csv)", class = "btn btn-primary btn-sm w-100 mb-2"),
          downloadButton("download_reps", "Replicate summary (.csv)", class = "btn btn-outline-primary btn-sm w-100 mb-2"),
          downloadButton("download_traj", "Trajectory plot (.png)", class = "btn btn-outline-primary btn-sm w-100 mb-2"),
          downloadButton("download_frac", "Memory plot (.png)", class = "btn btn-outline-primary btn-sm w-100")
        ),
        accordion(
          accordion_panel("Primary run parameters", parameter_panel(title = "Primary run"))
        )
      ),
      div(
        layout_column_wrap(
          width = 1/4,
          class = "mb-3",
          card(class = "metric-card", card_body(div(class="metric-title", "Final U"), div(class="metric-value", textOutput("kpi_u")), div(class="metric-sub", "Unexposed fraction at displayed time"))),
          card(class = "metric-card", card_body(div(class="metric-title", "Final M"), div(class="metric-value", textOutput("kpi_m")), div(class="metric-sub", "Memory fraction at displayed time"))),
          card(class = "metric-card", card_body(div(class="metric-title", "Memory fraction"), div(class="metric-value", textOutput("kpi_frac")), div(class="metric-sub", "M / (U + M)"))),
          card(class = "metric-card", card_body(div(class="metric-title", "Analytic M*"), div(class="metric-value", textOutput("kpi_eq")), div(class="metric-sub", "Shown only for the no-stress case")))
        ),
        navset_card_tab(
          height = "auto",
          nav_panel(
            "Overview",
            layout_column_wrap(
              width = 1/2,
              card(full_screen = TRUE, card_header("Population trajectories"), card_body(plotOutput("traj_plot_overview", height = "500px"))),
              card(full_screen = TRUE, card_header("Memory fraction over time"), card_body(plotOutput("frac_plot_overview", height = "500px")))
            ),
            card(full_screen = TRUE, class = "mt-3", card_header("Stress forcing"), card_body(plotOutput("stress_plot_overview", height = "300px")))
          ),
          nav_panel(
            "Population trajectories",
            card(full_screen = TRUE, card_header("Population trajectories"), card_body(plotOutput("traj_plot", height = "700px")))
          ),
          nav_panel(
            "Memory fraction",
            card(full_screen = TRUE, card_header("Memory fraction over time"), card_body(plotOutput("frac_plot", height = "700px")))
          ),
          nav_panel(
            "Stress forcing",
            card(full_screen = TRUE, card_header("Stress forcing"), card_body(plotOutput("stress_plot", height = "700px")))
          ),
          nav_panel(
            "Phase portrait",
            card(full_screen = TRUE, card_header("Phase portrait"), card_body(plotOutput("phase_plot", height = "700px")))
          )
        )
      )
    )
  ),
  
  nav_panel(
    "Compare runs",
    div(class = "app-hero", h3("Side-by-side comparison"), p("Run A and Run B stay available, but the plots are separated into tabs so each one has room to breathe.")),
    layout_sidebar(
      sidebar = sidebar(
        width = 420,
        open = "desktop",
        accordion(
          accordion_panel("Run A", parameter_panel(prefix = "a", title = "Run A", compare = TRUE)),
          accordion_panel("Run B", parameter_panel(prefix = "b", title = "Run B", compare = TRUE))
        )
      ),
      navset_card_tab(
        nav_panel(
          "Memory fraction",
          card(full_screen = TRUE, card_header("Memory fraction comparison"), card_body(plotOutput("compare_frac_plot", height = "720px")))
        ),
        nav_panel(
          "Population comparison",
          card(full_screen = TRUE, card_header("Population comparison"), card_body(plotOutput("compare_state_plot", height = "720px")))
        ),
        nav_panel(
          "Stress forcing",
          card(full_screen = TRUE, card_header("Stress forcing comparison"), card_body(plotOutput("compare_stress_plot", height = "720px")))
        )
      )
    )
  ),
  
  nav_panel(
    "Inference placeholder",
    card(
      class = "m-3 border-0 shadow-sm",
      card_body(
        h3("Future pomp fitting workflow"),
        p("This tab is still a placeholder for the next phase: attaching data and fitting the mechanistic model with pomp-style inference."),
        tags$ul(
          tags$li("Observed data input: population counts, memory fractions, or single-cell summary statistics."),
          tags$li("Measurement model: how noisy observations relate to latent U and M states."),
          tags$li("Parameter estimation: iterated filtering, particle filtering, PMCMC, or likelihood profiling."),
          tags$li("Model comparison: periodic forcing vs continuous Hill forcing vs extended models."),
          tags$li("Posterior or profile-based uncertainty on key biological parameters like epsilon, mu, and beta.")
        ),
        p(class = "small-muted", "Right now this app is for simulation, intuition, and design discussion. The next serious step is attaching real data and fitting."),
        h4("Suggested next coding targets"),
        tags$ol(
          tags$li("Create a synthetic observation layer with binomial or Gaussian noise."),
          tags$li("Wrap the model in pomp components: rprocess, rmeasure, dmeasure, initializer."),
          tags$li("Add a data upload widget and start with one fitting objective."),
          tags$li("Use this simulator tab to decide which parameters are practically identifiable.")
        )
      )
    )
  ),
  
  nav_panel(
    "Equations & notes",
    card(
      class = "m-3 border-0 shadow-sm",
      card_body(
        h3("Model structure"),
        p(HTML("This app implements the current draft model more directly: <b>U</b> (unexposed) and <b>M</b> (memory).")),
        tags$ul(
          tags$li(HTML("U grows at rate <i>r</i><sub>0</sub> with logistic crowding.")),
          tags$li(HTML("M grows at rate <i>c r</i><sub>0</sub>, allowing a memory cost or benefit.")),
          tags$li(HTML("Baseline switching is symmetric with terms <i>-&mu;(U-M)</i> and <i>-&mu;(M-U)</i>.")),
          tags$li(HTML("Stress-driven switching adds <i>&beta;U</i> or <i>&beta;&rho;(s)U</i>, while stress also penalizes U growth through <i>&alpha;</i>.")),
          tags$li(HTML("Memory decays back to U at rate <i>&epsilon;</i>.")),
          tags$li("Optional stochasticity is added with Euler-style noise, as in the Rmd prototype.")
        ),
        h4("Why this layout changed"),
        p("The previous version compressed too many outputs into one screen. This version separates the plots into tabs and gives each one a much larger drawing area so it is easier to read curves, legends, and axes."),
        h4("Suggested uses"),
        tags$ul(
          tags$li("Use Overview for fast scanning."),
          tags$li("Switch to a dedicated plot tab when you want to inspect one graph carefully."),
          tags$li("Use full-screen mode on each card during meetings."),
          tags$li("Use the compare tab to isolate how one parameter changes the model outcome.")
        )
      )
    )
  )
)

make_theme <- function() {
  theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#e2e8f0"),
      panel.grid.major.y = element_line(color = "#e2e8f0"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "top"
    )
}

plot_traj_gg <- function(dat, rs = NULL, view_time = NULL,
                         labels = c(U = "U", M = "M", total = "total")) {
  
  long <- dat %>%
    select(time, U, M, total) %>%
    pivot_longer(cols = c(U, M, total), names_to = "state", values_to = "value") %>%
    mutate(
      state = factor(state, levels = c("U", "M", "total")),
      state_label = recode(as.character(state),
                           U = labels["U"],
                           M = labels["M"],
                           total = labels["total"])
    )
  
  g <- ggplot() +
    make_theme() +
    labs(x = "Time (hours)", y = "Population fraction", color = NULL, fill = NULL)
  
  if (!is.null(rs)) {
    band_long <- bind_rows(
      rs %>% transmute(time, state = "U", lo = U_lo, hi = U_hi),
      rs %>% transmute(time, state = "M", lo = M_lo, hi = M_hi),
      rs %>% transmute(time, state = "total", lo = total_lo, hi = total_hi)
    ) %>%
      mutate(
        state = factor(state, levels = c("U", "M", "total")),
        state_label = recode(as.character(state),
                             U = labels["U"],
                             M = labels["M"],
                             total = labels["total"])
      )
    
    g <- g +
      geom_ribbon(
        data = band_long,
        aes(time, ymin = lo, ymax = hi, fill = state_label),
        alpha = 0.12
      ) +
      scale_fill_manual(values = setNames(c("#2563eb", "#db2777", "#ca8a04"), labels))
  }
  
  g +
    geom_line(data = long, aes(time, value, color = state_label), linewidth = 1.15) +
    scale_color_manual(values = setNames(c("#2563eb", "#db2777", "#ca8a04"), labels)) +
    {if (!is.null(view_time)) geom_vline(xintercept = view_time, linetype = "dashed", color = "#475569")} +
    coord_cartesian(xlim = c(0, max(dat$time)), ylim = c(0, 1))
}

plot_frac_gg <- function(dat, rs = NULL, view_time = NULL, eq = NA_real_,
                         labels = c(U = "U", M = "M")) {
  g <- ggplot() +
    make_theme() +
    labs(x = "Time (hours)", y = paste0(labels["M"], " / (", labels["U"], " + ", labels["M"], ")"))
  
  if (!is.null(rs)) {
    g <- g + geom_ribbon(data = rs, aes(time, ymin = mem_lo, ymax = mem_hi),
                         fill = "#10b981", alpha = 0.16)
  }
  
  g +
    geom_line(data = dat, aes(time, memory_fraction), color = "#10b981", linewidth = 1.2) +
    geom_point(data = dat %>% slice_tail(n = 1), aes(time, memory_fraction), color = "#06b6d4", size = 2.8) +
    {if (!is.na(eq)) geom_hline(yintercept = eq, linetype = "dotted", color = "#f59e0b", linewidth = 1)} +
    {if (!is.null(view_time)) geom_vline(xintercept = view_time, linetype = "dashed", color = "#475569")} +
    coord_cartesian(xlim = c(0, max(dat$time)), ylim = c(0, 1))
}

plot_stress_gg <- function(dat, scenario, view_time = NULL) {
  yvar <- if (scenario == "Continuous stress") "stress_raw" else "beta_eff"
  lab <- if (scenario == "Continuous stress") "Stress level s(t)" else expression(beta[eff](t))
  ggplot(dat, aes(time, .data[[yvar]])) +
    geom_line(color = "#8b5cf6", linewidth = 1.2) +
    geom_point(data = dat %>% slice_tail(n = 1), color = "#06b6d4", size = 2.8) +
    {if (!is.null(view_time)) geom_vline(xintercept = view_time, linetype = "dashed", color = "#475569")} +
    make_theme() +
    coord_cartesian(xlim = c(0, max(dat$time))) +
    labs(x = "Time (hours)", y = lab)
}

plot_phase_gg <- function(dat, labels = c(U = "U", M = "M")) {
  ggplot(dat, aes(U, M)) +
    geom_path(color = "#ef4444", linewidth = 1.2) +
    geom_point(data = dat %>% slice_head(n = 1), color = "#334155", size = 2.6) +
    geom_point(data = dat %>% slice_tail(n = 1), color = "#06b6d4", size = 3.2) +
    make_theme() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = labels["U"], y = labels["M"])
}

server <- function(input, output, session) {
  observeEvent(input$preset_no, apply_preset(session, preset_specs$no_stress))
  observeEvent(input$preset_p8, apply_preset(session, preset_specs$periodic_8h))
  observeEvent(input$preset_p24, apply_preset(session, preset_specs$periodic_24h))
  observeEvent(input$preset_cont, apply_preset(session, preset_specs$continuous_hill))
  
  observeEvent(input$preset_no_a, apply_preset(session, preset_specs$no_stress, "_a"))
  observeEvent(input$preset_p8_a, apply_preset(session, preset_specs$periodic_8h, "_a"))
  observeEvent(input$preset_p24_a, apply_preset(session, preset_specs$periodic_24h, "_a"))
  observeEvent(input$preset_cont_a, apply_preset(session, preset_specs$continuous_hill, "_a"))
  
  observeEvent(input$preset_no_b, apply_preset(session, preset_specs$no_stress, "_b"))
  observeEvent(input$preset_p8_b, apply_preset(session, preset_specs$periodic_8h, "_b"))
  observeEvent(input$preset_p24_b, apply_preset(session, preset_specs$periodic_24h, "_b"))
  observeEvent(input$preset_cont_b, apply_preset(session, preset_specs$continuous_hill, "_b"))
  
  for (id in c("", "_a", "_b")) {
    local({
      suffix <- id
      observe({
        total_init <- input[[paste0("U0", suffix)]] + input[[paste0("M0", suffix)]]
        if (!is.na(total_init) && total_init > 1) {
          showNotification("Initial U + M exceeds 1. The simulator rescales internally to keep the system biologically valid.", type = "message", duration = 4)
        }
      })
    })
  }
  
  display_labels <- reactive({
    get_display_labels(input$label_U, input$label_M)
  })
  
  primary_params <- reactive(build_params(input))
  compare_params_a <- reactive(build_params(input, "_a"))
  compare_params_b <- reactive(build_params(input, "_b"))
  
  sim_data <- reactive({
    do.call(simulate_memory, primary_params())
  })
  
  rep_data <- reactive({
    req(input$show_bands)
    pars <- primary_params()
    if (!isTRUE(pars$stochastic)) return(NULL)
    do.call(simulate_replicates, c(list(reps = input$replicates, seed_base = pars$seed), pars[names(pars) != "seed"]))
  })
  
  rep_summary <- reactive({
    dat <- rep_data()
    if (is.null(dat)) return(NULL)
    summarize_replicates(dat)
  })
  
  sim_a <- reactive(do.call(simulate_memory, compare_params_a()))
  sim_b <- reactive(do.call(simulate_memory, compare_params_b()))
  
  observeEvent(input$tmax, {
    current <- isolate(input$view_time %||% input$tmax)
    updateSliderInput(session, "view_time", max = input$tmax, value = min(current, input$tmax))
  }, ignoreInit = TRUE)
  
  observeEvent(sim_data(), {
    if (isTRUE(input$follow_end)) updateSliderInput(session, "view_time", value = max(sim_data()$time))
  }, ignoreInit = TRUE)
  
  current_row <- reactive(sim_data() %>% filter(time <= input$view_time) %>% slice_tail(n = 1))
  
  output$kpi_u <- renderText(percent(current_row()$U, accuracy = 0.1))
  output$kpi_m <- renderText(percent(current_row()$M, accuracy = 0.1))
  output$kpi_frac <- renderText(percent(current_row()$memory_fraction, accuracy = 0.1))
  output$kpi_eq <- renderText({
    eq <- calc_equilibrium(primary_params()$scenario, primary_params()$mu, primary_params()$eps, primary_params()$beta)
    if (is.na(eq)) "N/A" else percent(eq, accuracy = 0.1)
  })
  
  output$traj_plot <- renderPlot({
    dat <- sim_data() %>% filter(time <= input$view_time)
    rs <- rep_summary()
    print(plot_traj_gg(dat, rs, input$view_time, labels = display_labels()))
  }, res = 110)
  
  output$traj_plot_overview <- renderPlot({
    dat <- sim_data() %>% filter(time <= input$view_time)
    rs <- rep_summary()
    print(plot_traj_gg(dat, rs, input$view_time, labels = display_labels()))
  }, res = 100)
  
  output$frac_plot <- renderPlot({
    dat <- sim_data() %>% filter(time <= input$view_time)
    rs <- rep_summary()
    eq <- calc_equilibrium(primary_params()$scenario, primary_params()$mu, primary_params()$eps, primary_params()$beta)
    print(plot_frac_gg(dat, rs, input$view_time, eq, labels = display_labels()))
  }, res = 110)
  
  output$frac_plot_overview <- renderPlot({
    dat <- sim_data() %>% filter(time <= input$view_time)
    rs <- rep_summary()
    eq <- calc_equilibrium(primary_params()$scenario, primary_params()$mu, primary_params()$eps, primary_params()$beta)
    print(plot_frac_gg(dat, rs, input$view_time, eq, labels = display_labels()))
  }, res = 100)
  
  output$stress_plot <- renderPlot({
    dat <- sim_data() %>% filter(time <= input$view_time)
    print(plot_stress_gg(dat, primary_params()$scenario, input$view_time))
  }, res = 110)
  
  output$stress_plot_overview <- renderPlot({
    dat <- sim_data() %>% filter(time <= input$view_time)
    print(plot_stress_gg(dat, primary_params()$scenario, input$view_time))
  }, res = 100)
  
  output$phase_plot <- renderPlot({
    dat <- sim_data() %>% filter(time <= input$view_time)
    print(plot_phase_gg(dat, labels = display_labels()))
  }, res = 110)
  
  output$compare_frac_plot <- renderPlot({
    dat <- bind_rows(
      sim_a() %>% mutate(run = "Run A"),
      sim_b() %>% mutate(run = "Run B")
    )
    g <- ggplot(dat, aes(time, memory_fraction, color = run)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("Run A" = "#2563eb", "Run B" = "#f59e0b")) +
      make_theme() + 
      labs(
        x = "Time (hours)",
        y = paste0(display_labels()["M"], " / (", display_labels()["U"], " + ", display_labels()["M"], ")"),
        color = NULL
      )
    print(g)
  }, res = 110)
  
  output$compare_state_plot <- renderPlot({
    labels <- display_labels()
    
    dat <- bind_rows(
      sim_a() %>% transmute(time, run = "Run A", U, M),
      sim_b() %>% transmute(time, run = "Run B", U, M)
    ) %>%
      pivot_longer(cols = c(U, M), names_to = "state", values_to = "value") %>%
      mutate(
        state = recode(state,
                       U = labels["U"],
                       M = labels["M"])
      )
    
    g <- ggplot(dat, aes(time, value, color = run, linetype = state)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("Run A" = "#2563eb", "Run B" = "#f59e0b")) +
      make_theme() + 
      labs(x = "Time (hours)", y = "Population fraction", color = NULL, linetype = NULL)
    print(g)
  }, res = 110)
  
  output$compare_stress_plot <- renderPlot({
    dat <- bind_rows(
      sim_a() %>% mutate(run = "Run A"),
      sim_b() %>% mutate(run = "Run B")
    )
    g <- ggplot(dat, aes(time, beta_eff, color = run)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("Run A" = "#2563eb", "Run B" = "#f59e0b")) +
      make_theme() + labs(x = "Time (hours)", y = expression(beta[eff](t)), color = NULL)
    print(g)
  }, res = 110)
  
  output$download_data <- downloadHandler(
    filename = function() paste0("memory_simulation_", Sys.Date(), ".csv"),
    content = function(file) write.csv(sim_data(), file, row.names = FALSE)
  )
  
  output$download_reps <- downloadHandler(
    filename = function() paste0("memory_replicate_summary_", Sys.Date(), ".csv"),
    content = function(file) {
      rs <- rep_summary()
      if (is.null(rs)) rs <- tibble(message = "Replicate bands unavailable unless stochastic mode and show_bands are enabled.")
      write.csv(rs, file, row.names = FALSE)
    }
  )
  
  output$download_traj <- downloadHandler(
    filename = function() paste0("memory_trajectory_", Sys.Date(), ".png"),
    content = function(file) {
      labels <- display_labels()
      dat <- sim_data()
      rs <- rep_summary()
      
      long <- dat %>%
        select(time, U, M, total) %>%
        pivot_longer(cols = c(U, M, total), names_to = "state", values_to = "value") %>%
        mutate(
          state = factor(state, levels = c("U", "M", "total")),
          state_label = recode(as.character(state),
                               U = labels["U"],
                               M = labels["M"],
                               total = labels["total"])
        )
      
      g <- ggplot() +
        theme_minimal(base_size = 14) +
        theme(
          plot.background = element_rect(fill = "#111827", color = NA),
          panel.background = element_rect(fill = "#111827", color = NA),
          legend.position = "top",
          text = element_text(color = "white"),
          axis.text = element_text(color = "white"),
          axis.title = element_text(color = "white"),
          legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white")
        )
      
      if (!is.null(rs)) {
        band_long <- bind_rows(
          rs %>% transmute(time, state = "U", lo = U_lo, hi = U_hi),
          rs %>% transmute(time, state = "M", lo = M_lo, hi = M_hi),
          rs %>% transmute(time, state = "total", lo = total_lo, hi = total_hi)
        ) %>%
          mutate(
            state = factor(state, levels = c("U", "M", "total")),
            state_label = recode(as.character(state),
                                 U = labels["U"],
                                 M = labels["M"],
                                 total = labels["total"])
          )
        
        g <- g +
          geom_ribbon(data = band_long, aes(time, ymin = lo, ymax = hi, fill = state_label), alpha = 0.12) +
          scale_fill_manual(values = setNames(c("#60a5fa", "#f472b6", "#facc15"), labels))
      }
      
      g <- g +
        geom_line(data = long, aes(time, value, color = state_label), linewidth = 1.15) +
        scale_color_manual(values = setNames(c("#60a5fa", "#f472b6", "#facc15"), labels)) +
        coord_cartesian(xlim = c(0, max(dat$time)), ylim = c(0, 1)) +
        labs(x = "Time (hours)", y = "Population fraction", color = NULL, fill = NULL)
      
      ggsave(file, g, width = 10, height = 6, dpi = 300, bg = "#111827")
    }
  )
  
  output$download_frac <- downloadHandler(
    filename = function() paste0("memory_fraction_", Sys.Date(), ".png"),
    content = function(file) {
      labels <- display_labels()
      dat <- sim_data()
      rs <- rep_summary()
      g <- ggplot() +
        theme_minimal(base_size = 14) +
        theme(
          plot.background = element_rect(fill = "#111827", color = NA),
          panel.background = element_rect(fill = "#111827", color = NA),
          text = element_text(color = "white"),
          axis.text = element_text(color = "white"),
          axis.title = element_text(color = "white")
        )
      if (!is.null(rs)) {
        g <- g + geom_ribbon(data = rs, aes(time, ymin = mem_lo, ymax = mem_hi), fill = "#34d399", alpha = 0.16)
      }
      g <- g +
        geom_line(data = dat, aes(time, memory_fraction), color = "#34d399", linewidth = 1.2) +
        coord_cartesian(xlim = c(0, max(dat$time)), ylim = c(0, 1)) +
        labs(x = "Time (hours)", y = paste0(labels["M"], " / (", labels["U"], " + ", labels["M"], ")"))
      
      ggsave(file, g, width = 9, height = 5.5, dpi = 300, bg = "#111827")
    }
  )
}

shinyApp(ui, server)