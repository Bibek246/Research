library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(htmltools)

clip01 <- function(x) pmax(0, pmin(1, x))
`%||%` <- function(x, y) if (is.null(x)) y else x
clamp <- function(x, lo, hi) max(lo, min(hi, x))

normalize_range <- function(minv, maxv, default_min, default_max, tiny = 1e-6) {
  minv <- suppressWarnings(as.numeric(minv))
  maxv <- suppressWarnings(as.numeric(maxv))
  if (is.na(minv)) minv <- default_min
  if (is.na(maxv)) maxv <- default_max
  if (minv >= maxv) {
    mid <- mean(c(minv, maxv))
    span <- max(abs(default_max - default_min), tiny)
    minv <- mid - span / 2
    maxv <- mid + span / 2
  }
  list(min = minv, max = maxv)
}

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
    stress <- pmax(0, s_base + s_amp * sin((times) / (forcing_period / (2 * pi))))
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
  
  total0 <- U0 + M0
  if (total0 > 1 && total0 > 0) {
    U0 <- U0 / total0
    M0 <- M0 / total0
  }
  
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


# -----------------------------
# Formaldehyde resistance model
# -----------------------------
# This model follows the slide-deck structure:
# da/dt = f(d) * H_A(d) * (1 + eta * xi) - L(d) * a
# dd/dt = Ki * (D_D(t) - d) - (Ka * a * d)/(ka + d) - L(d) * d
# where a = resistance/enzyme level and d = intracellular formaldehyde/damage.

hill_response <- function(x, k, n) {
  x <- pmax(0, x)
  k <- max(k, 1e-12)
  n <- max(n, 1e-12)
  (x^n) / (k^n + x^n + 1e-12)
}

formaldehyde_external <- function(times,
                                  D1 = 0.50,
                                  pulse_start = 150,
                                  pulse_duration = 150,
                                  include_retest = TRUE,
                                  D2 = 0.50,
                                  retest_start = 300,
                                  retest_duration = 40) {
  DD <- ifelse(times >= pulse_start & times <= (pulse_start + pulse_duration), D1, 0)
  if (isTRUE(include_retest)) {
    DD2 <- ifelse(times >= retest_start & times <= (retest_start + retest_duration), D2, 0)
    DD <- pmax(DD, DD2)
  }
  DD
}

formaldehyde_growth <- function(d, Lmax = 0.012, L50 = 0.40, L_n = 3) {
  Lmax / (1 + (pmax(0, d) / max(L50, 1e-12))^max(L_n, 1e-12))
}

formaldehyde_expression <- function(d,
                                    basal = 0.0002,
                                    hmax = 0.010,
                                    K_expr = 0.30,
                                    h_n = 4,
                                    f_amp = 0.50,
                                    f50 = 0.40,
                                    f_n = 2) {
  H_A <- basal + hmax * hill_response(d, K_expr, h_n)
  f_mod <- 1 + f_amp * hill_response(d, f50, f_n)
  H_A * f_mod
}

simulate_formaldehyde <- function(tmax = 500,
                                  dt = 0.5,
                                  a0 = 0.25,
                                  d0 = 0,
                                  D1 = 0.50,
                                  pulse_start = 150,
                                  pulse_duration = 150,
                                  include_retest = TRUE,
                                  D2 = 0.50,
                                  retest_start = 300,
                                  retest_duration = 40,
                                  Ki = 0.030,
                                  Ka = 0.070,
                                  ka = 0.20,
                                  Lmax = 0.012,
                                  L50 = 0.40,
                                  L_n = 3,
                                  basal = 0.0002,
                                  hmax = 0.010,
                                  K_expr = 0.30,
                                  h_n = 4,
                                  f_amp = 0.50,
                                  f50 = 0.40,
                                  f_n = 2,
                                  noise_amp = 0.08,
                                  stochastic = TRUE,
                                  seed = 123,
                                  resistance_threshold = 0.50) {
  dt <- max(dt, 1e-6)
  tmax <- max(tmax, dt)
  times <- seq(0, tmax, by = dt)
  n <- length(times)
  DD <- formaldehyde_external(times, D1, pulse_start, pulse_duration, include_retest, D2, retest_start, retest_duration)
  
  a <- numeric(n)
  d <- numeric(n)
  growth <- numeric(n)
  expression <- numeric(n)
  processing <- numeric(n)
  diffusion <- numeric(n)
  
  a[1] <- max(0, a0)
  d[1] <- max(0, d0)
  
  if (!is.null(seed) && stochastic) set.seed(seed)
  
  for (i in 2:n) {
    a_prev <- max(0, a[i - 1])
    d_prev <- max(0, d[i - 1])
    DD_prev <- DD[i - 1]
    
    L_prev <- formaldehyde_growth(d_prev, Lmax, L50, L_n)
    expr_prev <- formaldehyde_expression(d_prev, basal, hmax, K_expr, h_n, f_amp, f50, f_n)
    noise_mult <- if (isTRUE(stochastic)) max(0, 1 + noise_amp * rnorm(1)) else 1
    proc_prev <- (Ka * a_prev * d_prev) / (ka + d_prev + 1e-12)
    diff_prev <- Ki * (DD_prev - d_prev)
    
    da_det <- expr_prev * noise_mult - L_prev * a_prev
    dd_det <- diff_prev - proc_prev - L_prev * d_prev
    
    a[i] <- max(0, a_prev + da_det * dt)
    d[i] <- max(0, d_prev + dd_det * dt)
  }
  
  growth <- formaldehyde_growth(d, Lmax, L50, L_n)
  expression <- formaldehyde_expression(d, basal, hmax, K_expr, h_n, f_amp, f50, f_n)
  processing <- (Ka * a * d) / (ka + d + 1e-12)
  diffusion <- Ki * (DD - d)
  
  tibble(
    time = times,
    a = a,
    d = d,
    DD = DD,
    growth = growth,
    relative_growth = ifelse(Lmax > 0, growth / Lmax, NA_real_),
    expression = expression,
    processing = processing,
    diffusion = diffusion,
    activated = a >= resistance_threshold,
    pulse_start = pulse_start,
    resistance_threshold = resistance_threshold
  )
}

calc_activation_lag <- function(dat) {
  if (is.null(dat) || nrow(dat) == 0) return(NA_real_)
  start_time <- dat$pulse_start[1]
  threshold <- dat$resistance_threshold[1]
  idx <- which(dat$time >= start_time & dat$a >= threshold)
  if (length(idx) == 0) NA_real_ else dat$time[min(idx)] - start_time
}

simulate_formaldehyde_replicates <- function(reps = 25, seed_base = 123, ...) {
  bind_rows(lapply(seq_len(reps), function(i) {
    dat <- simulate_formaldehyde(seed = seed_base + i - 1, ...)
    dat$replicate <- i
    dat
  }))
}

summarize_formaldehyde_replicates <- function(rep_dat) {
  rep_dat %>%
    group_by(time) %>%
    summarize(
      a_med = median(a),
      a_lo = quantile(a, 0.1),
      a_hi = quantile(a, 0.9),
      d_med = median(d),
      d_lo = quantile(d, 0.1),
      d_hi = quantile(d, 0.9),
      growth_med = median(relative_growth, na.rm = TRUE),
      growth_lo = quantile(relative_growth, 0.1, na.rm = TRUE),
      growth_hi = quantile(relative_growth, 0.9, na.rm = TRUE),
      expression_med = median(expression, na.rm = TRUE),
      expression_lo = quantile(expression, 0.1, na.rm = TRUE),
      expression_hi = quantile(expression, 0.9, na.rm = TRUE),
      processing_med = median(processing, na.rm = TRUE),
      processing_lo = quantile(processing, 0.1, na.rm = TRUE),
      processing_hi = quantile(processing, 0.9, na.rm = TRUE),
      diffusion_med = median(diffusion, na.rm = TRUE),
      diffusion_lo = quantile(diffusion, 0.1, na.rm = TRUE),
      diffusion_hi = quantile(diffusion, 0.9, na.rm = TRUE),
      DD = first(DD),
      .groups = "drop"
    )
}

summarize_formaldehyde_lags <- function(rep_dat) {
  if (is.null(rep_dat) || nrow(rep_dat) == 0) return(tibble())
  rep_dat %>%
    group_by(replicate) %>%
    group_modify(~ tibble(
      activation_lag = calc_activation_lag(.x),
      final_resistance = last(.x$a),
      max_intracellular_formaldehyde = max(.x$d),
      min_relative_growth = min(.x$relative_growth, na.rm = TRUE)
    )) %>%
    ungroup()
}

build_formaldehyde_params <- function(input) {
  list(
    tmax = input$form_tmax,
    dt = input$form_dt,
    a0 = input$form_a0,
    d0 = input$form_d0,
    D1 = input$form_D1,
    pulse_start = input$form_pulse_start,
    pulse_duration = input$form_pulse_duration,
    include_retest = input$form_include_retest,
    D2 = input$form_D2 %||% input$form_D1,
    retest_start = input$form_retest_start %||% 300,
    retest_duration = input$form_retest_duration %||% 40,
    Ki = input$form_Ki,
    Ka = input$form_Ka,
    ka = input$form_ka,
    Lmax = input$form_Lmax,
    L50 = input$form_L50,
    L_n = input$form_Ln,
    basal = input$form_basal,
    hmax = input$form_hmax,
    K_expr = input$form_Kexpr,
    h_n = input$form_Hn,
    f_amp = input$form_famp,
    f50 = input$form_f50,
    f_n = input$form_Fn,
    noise_amp = if (isTRUE(input$form_stochastic)) input$form_noise else 0,
    stochastic = input$form_stochastic,
    seed = input$form_seed,
    resistance_threshold = input$form_threshold
  )
}

formaldehyde_parameter_panel <- function() {
  tagList(
    div(
      class = "control-block",
      h5("Simulation"),
      sliderInput("form_tmax", "Simulation horizon (minutes)", min = 60, max = 1200, value = 500, step = 10),
      sliderInput("form_dt", "Time step (minutes)", min = 0.05, max = 5, value = 0.50, step = 0.05),
      sliderInput("form_view_time", "Displayed time", min = 0, max = 500, value = 500, step = 1,
                  animate = animationOptions(interval = 350, loop = FALSE)),
      numericInput("form_seed", "Random seed", value = 123, min = 1, step = 1),
      checkboxInput("form_stochastic", "Include expression noise", TRUE),
      conditionalPanel("input.form_stochastic == true",
                       sliderInput("form_noise", HTML("Expression noise amplitude (&eta;)"), min = 0, max = 0.50, value = 0.08, step = 0.01)),
      checkboxInput("form_show_bands", "Show replicate uncertainty bands", TRUE),
      sliderInput("form_replicates", "Number of replicate runs", min = 5, max = 100, value = 30, step = 5)
    ),
    div(
      class = "control-block",
      h5("External formaldehyde input"),
      sliderInput("form_D1", HTML("Primary extracellular formaldehyde D<sub>D</sub>(t)"), min = 0, max = 2, value = 0.50, step = 0.05),
      sliderInput("form_pulse_start", "Primary exposure start (min)", min = 0, max = 900, value = 150, step = 5),
      sliderInput("form_pulse_duration", "Primary exposure duration (min)", min = 1, max = 500, value = 150, step = 5),
      checkboxInput("form_include_retest", "Include retest pulse", TRUE),
      conditionalPanel(
        "input.form_include_retest == true",
        sliderInput("form_D2", HTML("Retest extracellular formaldehyde D<sub>D</sub>(t)"), min = 0, max = 2, value = 0.50, step = 0.05),
        sliderInput("form_retest_start", "Retest start (min)", min = 0, max = 1000, value = 300, step = 5),
        sliderInput("form_retest_duration", "Retest duration (min)", min = 1, max = 300, value = 40, step = 5)
      )
    ),
    div(
      class = "control-block",
      h5("Initial conditions and activation"),
      sliderInput("form_a0", "Initial resistance/enzyme level (a0)", min = 0, max = 1.5, value = 0.25, step = 0.01),
      sliderInput("form_d0", "Initial intracellular formaldehyde/damage (d0)", min = 0, max = 1.5, value = 0, step = 0.01),
      sliderInput("form_threshold", "Resistance activation threshold", min = 0.05, max = 1.5, value = 0.50, step = 0.05)
    ),
    div(
      class = "control-block",
      h5("Diffusion and processing"),
      sliderInput("form_Ki", HTML("Diffusion rate K<sub>i</sub>"), min = 0, max = 0.20, value = 0.030, step = 0.005),
      sliderInput("form_Ka", HTML("Enzyme Vmax K<sub>a</sub>"), min = 0, max = 0.30, value = 0.070, step = 0.005),
      sliderInput("form_ka", HTML("Michaelis constant k<sub>a</sub>"), min = 0.01, max = 2, value = 0.20, step = 0.01)
    ),
    div(
      class = "control-block",
      h5("Growth and dilution"),
      sliderInput("form_Lmax", HTML("Maximum growth/dilution L<sub>max</sub>"), min = 0.001, max = 0.05, value = 0.012, step = 0.001),
      sliderInput("form_L50", HTML("Growth inhibition midpoint L<sub>50</sub>"), min = 0.05, max = 2, value = 0.40, step = 0.05),
      sliderInput("form_Ln", "Growth inhibition Hill coefficient", min = 1, max = 10, value = 3, step = 1)
    ),
    div(
      class = "control-block",
      h5("Enzyme expression"),
      sliderInput("form_basal", "Basal expression", min = 0, max = 0.005, value = 0.0002, step = 0.0001),
      sliderInput("form_hmax", HTML("Maximum induced expression H<sub>A,max</sub>"), min = 0, max = 0.05, value = 0.010, step = 0.001),
      sliderInput("form_Kexpr", HTML("Expression midpoint K<sub>expr</sub>"), min = 0.01, max = 2, value = 0.30, step = 0.01),
      sliderInput("form_Hn", "Expression Hill coefficient", min = 1, max = 10, value = 4, step = 1),
      sliderInput("form_famp", HTML("Expression modulation strength f<sub>amp</sub>"), min = 0, max = 5, value = 0.50, step = 0.05),
      sliderInput("form_f50", HTML("Modulation midpoint f<sub>50</sub>"), min = 0.01, max = 2, value = 0.40, step = 0.01),
      sliderInput("form_Fn", "Modulation Hill coefficient", min = 1, max = 10, value = 2, step = 1)
    ),
    div(
      class = "control-block",
      h5("Export"),
      downloadButton("download_form_data", "Formaldehyde simulation (.csv)", class = "btn btn-primary btn-sm w-100 mb-2"),
      downloadButton("download_form_reps", "Replicate lag summary (.csv)", class = "btn btn-outline-primary btn-sm w-100 mb-2"),
      downloadButton("download_form_plot", "Concentration plot (.png)", class = "btn btn-outline-primary btn-sm w-100")
    )
  )
}

plot_formaldehyde_conc_gg <- function(dat, rs = NULL, view_time = NULL) {
  g <- ggplot() +
    make_theme() +
    labs(x = "Time (minutes)", y = "Concentration / model units", color = NULL, fill = NULL)
  
  if (!is.null(rs)) {
    band_long <- bind_rows(
      rs %>% transmute(time, quantity = "Resistance/enzyme a", lo = a_lo, hi = a_hi, med = a_med),
      rs %>% transmute(time, quantity = "Intracellular formaldehyde d", lo = d_lo, hi = d_hi, med = d_med)
    )
    g <- g +
      geom_ribbon(data = band_long, aes(time, ymin = lo, ymax = hi, fill = quantity), alpha = 0.14) +
      geom_line(data = band_long, aes(time, med, color = quantity), linewidth = 1.15)
  } else {
    long <- dat %>%
      select(time, a, d) %>%
      pivot_longer(cols = c(a, d), names_to = "quantity", values_to = "value") %>%
      mutate(quantity = recode(quantity, a = "Resistance/enzyme a", d = "Intracellular formaldehyde d"))
    g <- g + geom_line(data = long, aes(time, value, color = quantity), linewidth = 1.15)
  }
  
  g +
    geom_line(data = dat, aes(time, DD, color = "Extracellular formaldehyde DD(t)"), linewidth = 1.0, linetype = "dashed") +
    geom_hline(aes(yintercept = first(dat$resistance_threshold), color = "Resistance threshold"), linetype = "dotted", linewidth = 0.9) +
    scale_color_manual(values = c(
      "Resistance/enzyme a" = "#dc2626",
      "Intracellular formaldehyde d" = "#2563eb",
      "Extracellular formaldehyde DD(t)" = "#475569",
      "Resistance threshold" = "#f59e0b"
    )) +
    scale_fill_manual(values = c(
      "Resistance/enzyme a" = "#dc2626",
      "Intracellular formaldehyde d" = "#2563eb"
    )) +
    {if (!is.null(view_time)) geom_vline(xintercept = view_time, linetype = "dashed", color = "#475569")} +
    coord_cartesian(xlim = c(0, max(dat$time)))
}

plot_formaldehyde_process_gg <- function(dat, rs = NULL, view_time = NULL) {
  if (!is.null(rs)) {
    gdat <- rs %>%
      transmute(time,
                `Relative growth L(d) / Lmax` = growth_med,
                `Enzyme expression f(d)H_A(d)` = expression_med,
                `Enzyme processing` = processing_med,
                `Diffusion flux` = diffusion_med,
                `Extracellular formaldehyde DD(t)` = DD) %>%
      pivot_longer(cols = -time, names_to = "quantity", values_to = "value")
  } else {
    gdat <- dat %>%
      transmute(time,
                `Relative growth L(d) / Lmax` = relative_growth,
                `Enzyme expression f(d)H_A(d)` = expression,
                `Enzyme processing` = processing,
                `Diffusion flux` = diffusion) %>%
      pivot_longer(cols = -time, names_to = "quantity", values_to = "value")
  }
  
  ggplot(gdat, aes(time, value)) +
    geom_line(color = "#0f766e", linewidth = 1.1) +
    {if (!is.null(view_time)) geom_vline(xintercept = view_time, linetype = "dashed", color = "#475569")} +
    facet_wrap(~quantity, scales = "free_y", ncol = 1) +
    make_theme() +
    theme(legend.position = "none") +
    coord_cartesian(xlim = c(0, max(dat$time))) +
    labs(x = "Time (minutes)", y = NULL)
}

plot_formaldehyde_lag_gg <- function(lag_dat) {
  lag_dat <- lag_dat %>% filter(!is.na(activation_lag))
  ggplot(lag_dat, aes(activation_lag)) +
    geom_histogram(bins = 20, fill = "#8b5cf6", color = "white", alpha = 0.9) +
    make_theme() +
    theme(legend.position = "none") +
    labs(x = "Activation lag after primary exposure starts (minutes)", y = "Replicate count")
}

plot_formaldehyde_phase_gg <- function(dat, view_time = NULL) {
  ggplot(dat, aes(d, a)) +
    geom_path(color = "#ef4444", linewidth = 1.15) +
    geom_point(data = dat %>% slice_head(n = 1), color = "#334155", size = 2.6) +
    geom_point(data = dat %>% slice_tail(n = 1), color = "#06b6d4", size = 3.2) +
    geom_hline(aes(yintercept = first(resistance_threshold)), linetype = "dotted", color = "#f59e0b") +
    make_theme() +
    labs(x = "Intracellular formaldehyde/damage d", y = "Resistance/enzyme level a")
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

range_editor_ui <- function(suffix = "") {
  ns <- function(id) paste0(id, suffix)
  div(
    class = "control-block",
    h5("Custom slider ranges"),
    p(class = "small-muted", "Research mode lets you redefine the range of the main sliders. Turning it off returns them to safe defaults."),
    checkboxInput(ns("advanced_ranges"), "Research mode: customize slider ranges", FALSE),
    conditionalPanel(
      sprintf("input['%s'] == true", ns("advanced_ranges")),
      tags$div(class = "range-grid",
               tags$div(tags$strong("Parameter"), tags$strong("Min"), tags$strong("Max")),
               tags$div(tags$span("Simulation horizon"), numericInput(ns("tmax_min"), NULL, 24, min = 1), numericInput(ns("tmax_max"), NULL, 240, min = 2)),
               tags$div(tags$span("Time step"), numericInput(ns("dt_min"), NULL, 0.01, min = 0.001, step = 0.01), numericInput(ns("dt_max"), NULL, 0.50, min = 0.01, step = 0.01)),
               tags$div(tags$span("Initial U"), numericInput(ns("U0_min"), NULL, 0, min = 0, max = 1, step = 0.01), numericInput(ns("U0_max"), NULL, 1, min = 0, max = 1, step = 0.01)),
               tags$div(tags$span("Initial M"), numericInput(ns("M0_min"), NULL, 0, min = 0, max = 1, step = 0.01), numericInput(ns("M0_max"), NULL, 1, min = 0, max = 1, step = 0.01)),
               tags$div(tags$span("r0"), numericInput(ns("r0_min"), NULL, 0.05, min = 0, step = 0.01), numericInput(ns("r0_max"), NULL, 2.0, min = 0.01, step = 0.01)),
               tags$div(tags$span("c_cost"), numericInput(ns("c_cost_min"), NULL, 0.10, step = 0.01), numericInput(ns("c_cost_max"), NULL, 1.20, step = 0.01)),
               tags$div(tags$span("alpha"), numericInput(ns("alpha_min"), NULL, 0, step = 0.01), numericInput(ns("alpha_max"), NULL, 1.0, step = 0.01)),
               tags$div(tags$span("mu"), numericInput(ns("mu_min"), NULL, 0, step = 0.001), numericInput(ns("mu_max"), NULL, 0.05, step = 0.001)),
               tags$div(tags$span("eps"), numericInput(ns("eps_min"), NULL, 0, step = 0.001), numericInput(ns("eps_max"), NULL, 0.20, step = 0.001)),
               tags$div(tags$span("beta"), numericInput(ns("beta_min"), NULL, 0, step = 0.005), numericInput(ns("beta_max"), NULL, 0.50, step = 0.005)),
               tags$div(tags$span("sigma"), numericInput(ns("sigma_min"), NULL, 0, step = 0.005), numericInput(ns("sigma_max"), NULL, 0.30, step = 0.005))
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
    range_editor_ui(suffix),
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
    .bslib-sidebar-layout > .main { padding-top: .5rem; }
    .equation-box{
      padding: .75rem .9rem;
      border-radius: 14px;
      background: #f8fafc;
      border: 1px solid #e2e8f0;
      font-family: 'Times New Roman', serif;
      font-size: 1.1rem;
      margin-bottom: .6rem;
    }

    .note-grid{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
      gap: .9rem;
      margin-top: .75rem;
      margin-bottom: .75rem;
    }
    .note-card{
      padding: .9rem 1rem;
      border-radius: 16px;
      background: #ffffff;
      border: 1px solid #e2e8f0;
      box-shadow: 0 2px 10px rgba(15,23,42,.04);
    }
    .note-card h5{ margin-bottom: .45rem; font-weight: 800; }
    .note-badge{
      display: inline-block;
      border-radius: 999px;
      padding: .15rem .55rem;
      font-size: .78rem;
      font-weight: 700;
      background: rgba(109,40,217,.10);
      color: #5b21b6;
      margin-bottom: .45rem;
    }
    .equation-box .eq-label{
      display: block;
      font-family: Inter, sans-serif;
      font-size: .78rem;
      text-transform: uppercase;
      letter-spacing: .06em;
      color: #64748b;
      margin-bottom: .25rem;
    }
    .equation-box .eq-main{
      display: block;
      color: #0f172a;
    }
    .doc-table th{ background: #f1f5f9; }
    .doc-table td, .doc-table th{ vertical-align: top; }
    .range-grid > div {
      display: grid;
      grid-template-columns: 1.2fr .9fr .9fr;
      gap: .4rem;
      align-items: center;
      margin-bottom: .35rem;
    }
    .range-grid .form-group { margin-bottom: 0; }
  "))),
  nav_panel(
    "Simulator",
    div(class = "app-hero", h3("Prototype"), p("Use the tabs below to focus on one plot at a time.")),
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
        accordion(accordion_panel("Primary run parameters", parameter_panel(title = "Primary run")))
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
          nav_panel("Population trajectories", card(full_screen = TRUE, card_header("Population trajectories"), card_body(plotOutput("traj_plot", height = "700px")))),
          nav_panel("Memory fraction", card(full_screen = TRUE, card_header("Memory fraction over time"), card_body(plotOutput("frac_plot", height = "700px")))),
          nav_panel("Stress forcing", card(full_screen = TRUE, card_header("Stress forcing"), card_body(plotOutput("stress_plot", height = "700px")))),
          nav_panel("Phase portrait", card(full_screen = TRUE, card_header("Phase portrait"), card_body(plotOutput("phase_plot", height = "700px"))))
        )
      )
    )
  ),
  nav_panel(
    "Formaldehyde model",
    div(
      class = "app-hero",
      h3("Mechanistic formaldehyde resistance model"),
      p("This page adds the intracellular resistance model from the formaldehyde slide deck: extracellular formaldehyde enters the cell, resistance/enzyme expression rises, enzyme processing removes intracellular formaldehyde, and growth dilution lowers both a and d.")
    ),
    layout_sidebar(
      sidebar = sidebar(
        width = 390,
        open = "desktop",
        formaldehyde_parameter_panel()
      ),
      div(
        layout_column_wrap(
          width = 1/4,
          class = "mb-3",
          card(class = "metric-card", card_body(div(class="metric-title", "Resistance a"), div(class="metric-value", textOutput("form_kpi_a")), div(class="metric-sub", "Resistance/enzyme level at displayed time"))),
          card(class = "metric-card", card_body(div(class="metric-title", "Intracellular d"), div(class="metric-value", textOutput("form_kpi_d")), div(class="metric-sub", "Formaldehyde/damage at displayed time"))),
          card(class = "metric-card", card_body(div(class="metric-title", "Relative growth"), div(class="metric-value", textOutput("form_kpi_growth")), div(class="metric-sub", "L(d) / Lmax at displayed time"))),
          card(class = "metric-card", card_body(div(class="metric-title", "Activation lag"), div(class="metric-value", textOutput("form_kpi_lag")), div(class="metric-sub", "First time a crosses threshold after exposure")))
        ),
        navset_card_tab(
          height = "auto",
          nav_panel(
            "Overview",
            layout_column_wrap(
              width = 1/2,
              card(full_screen = TRUE, card_header("Formaldehyde and resistance trajectories"), card_body(plotOutput("form_conc_plot_overview", height = "520px"))),
              card(full_screen = TRUE, card_header("Growth, expression, processing"), card_body(plotOutput("form_process_plot_overview", height = "520px")))
            )
          ),
          nav_panel("Concentrations", card(full_screen = TRUE, card_header("Resistance a, intracellular formaldehyde d, and extracellular input"), card_body(plotOutput("form_conc_plot", height = "720px")))),
          nav_panel("Model processes", card(full_screen = TRUE, card_header("Growth, expression, enzyme processing, and diffusion"), card_body(plotOutput("form_process_plot", height = "720px")))),
          nav_panel("Lag distribution", card(full_screen = TRUE, card_header("Activation lag across stochastic replicates"), card_body(plotOutput("form_lag_plot", height = "720px")))),
          nav_panel("Phase portrait", card(full_screen = TRUE, card_header("Phase portrait: intracellular formaldehyde d vs resistance a"), card_body(plotOutput("form_phase_plot", height = "720px")))),
          nav_panel(
            "Equations",
            card(
              class = "border-0 shadow-sm",
              card_body(
                h4("Implemented ODE structure"),
                p(HTML("<b>a</b> is the resistance/enzyme level and <b>d</b> is intracellular formaldehyde/damage.")),
                tags$div(class = "equation-box", HTML("da/dt = f(d)H<sub>A</sub>(d)(1 + &eta;&xi;) - L(d)a")),
                tags$div(class = "equation-box", HTML("dd/dt = K<sub>i</sub>(D<sub>D</sub>(t) - d) - K<sub>a</sub>ad/(k<sub>a</sub> + d) - L(d)d")),
                tags$ul(
                  tags$li(HTML("<b>D<sub>D</sub>(t)</b>: extracellular formaldehyde input pulse.")),
                  tags$li(HTML("<b>K<sub>i</sub>(D<sub>D</sub>(t)-d)</b>: diffusion into/out of the cell.")),
                  tags$li(HTML("<b>K<sub>a</sub>ad/(k<sub>a</sub>+d)</b>: enzyme processing/removal of intracellular formaldehyde.")),
                  tags$li(HTML("<b>L(d)</b>: growth-dependent dilution, reduced when d is high.")),
                  tags$li(HTML("<b>f(d)H<sub>A</sub>(d)</b>: damage-induced resistance/enzyme expression.")),
                  tags$li(HTML("<b>&eta;&xi;</b>: optional expression noise to create cell-to-cell lag variation."))
                ),
                p(class = "small-muted", "The parameter values here are adjustable working defaults for simulation and discussion, not final fitted biological estimates.")
              )
            )
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
        nav_panel("Memory fraction", card(full_screen = TRUE, card_header("Memory fraction comparison"), card_body(plotOutput("compare_frac_plot", height = "720px")))),
        nav_panel("Population comparison", card(full_screen = TRUE, card_header("Population comparison"), card_body(plotOutput("compare_state_plot", height = "720px")))),
        nav_panel("Stress forcing", card(full_screen = TRUE, card_header("Stress forcing comparison"), card_body(plotOutput("compare_stress_plot", height = "720px"))))
      )
    )
  ),
  nav_panel(
    "Inference placeholder",
    card(class = "m-3 border-0 shadow-sm", card_body(
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
    ))
  ),
  nav_panel(
    "Equations & notes",
    div(
      class = "app-hero",
      h3("Model documentation and research notes"),
      p("This tab now documents both model layers in the app: the U/M population-memory model and the mechanistic formaldehyde-resistance model. Use this page during meetings to explain what each parameter means, how the equations connect to the biology, and how to check whether the implementation matches the research model.")
    ),
    navset_card_tab(
      height = "auto",
      nav_panel(
        "Research overview",
        div(
          class = "m-3",
          div(
            class = "note-grid",
            div(
              class = "note-card",
              span(class = "note-badge", "Population memory model"),
              h5("What it answers"),
              p("How does a population move between a normal/unexposed state and a memory state after environmental stress?"),
              tags$ul(
                tags$li(HTML("Tracks <b>U</b>: normal or unexposed cells.")),
                tags$li(HTML("Tracks <b>M</b>: cells carrying the memory phenotype.")),
                tags$li("Useful for parameter-space exploration: switching, memory decay, stress intensity, and fitness cost.")
              )
            ),
            div(
              class = "note-card",
              span(class = "note-badge", "Formaldehyde model"),
              h5("What it answers"),
              p("How can intracellular formaldehyde, resistance/enzyme expression, growth dilution, and stochasticity create delayed resistance activation?"),
              tags$ul(
                tags$li(HTML("Tracks <b>a</b>: resistance/enzyme level.")),
                tags$li(HTML("Tracks <b>d</b>: intracellular formaldehyde or damage.")),
                tags$li("Useful for studying lag time, enzyme processing, and retest behavior.")
              )
            ),
            div(
              class = "note-card",
              span(class = "note-badge", "Long-run goal"),
              h5("How the two pieces fit together"),
              p("The population model is good for bulk memory dynamics. The formaldehyde model explains one possible intracellular mechanism that can generate resistant or memory-like cells. Later, the two can be connected by making the transition into M depend on whether a crosses a resistance threshold.")
            )
          ),
          card(
            class = "border-0 shadow-sm mt-3",
            card_body(
              h4("Main biological questions represented in the app"),
              tags$ol(
                tags$li("Why does a stress-induced phenotype persist for multiple generations instead of disappearing immediately through dilution?"),
                tags$li("Can the observed behavior be explained by discrete state switching, continuous dilution, delayed expression, or a combination of these?"),
                tags$li("For formaldehyde resistance, how do intracellular formaldehyde, enzyme expression, and growth-dependent dilution affect activation lag?"),
                tags$li("Which outputs should be compared to experiments: mean behavior, single-cell variance, lag-time distributions, and recovery after retest?")
              )
            )
          )
        )
      ),
      nav_panel(
        "Memory model equations",
        div(
          class = "m-3",
          card(
            class = "border-0 shadow-sm",
            card_body(
              h4("U/M population-memory model"),
              p(HTML("This is the model used in the <b>Simulator</b> and <b>Compare runs</b> tabs. It tracks a normal/unexposed state <b>U</b> and a memory state <b>M</b>.")),
              tags$div(class = "equation-box", HTML("<span class='eq-label'>Unexposed state</span><span class='eq-main'>dU/dt = r<sub>U</sub>(t)(1 - U - M)U - &beta;<sub>eff</sub>(t)U - &mu;(U - M) + &epsilon;M</span>")),
              tags$div(class = "equation-box", HTML("<span class='eq-label'>Memory state</span><span class='eq-main'>dM/dt = r<sub>M</sub>(1 - U - M)M + &beta;<sub>eff</sub>(t)U - &mu;(M - U) - &epsilon;M</span>")),
              tags$div(class = "equation-box", HTML("<span class='eq-label'>Memory fraction</span><span class='eq-main'>Memory fraction = M / (U + M)</span>")),
              h5("Parameter meanings"),
              tags$table(
                class = "table table-sm table-striped doc-table",
                tags$thead(tags$tr(tags$th("Term"), tags$th("Meaning"), tags$th("Where it appears in the app"))),
                tags$tbody(
                  tags$tr(tags$td(HTML("<b>U</b>")), tags$td("Normal/unexposed population state."), tags$td("Population trajectories, phase portrait, compare tab.")),
                  tags$tr(tags$td(HTML("<b>M</b>")), tags$td("Memory phenotype state."), tags$td("Population trajectories, memory fraction plot.")),
                  tags$tr(tags$td(HTML("r<sub>0</sub>")), tags$td("Base growth rate of U."), tags$td("Growth rate of U slider.")),
                  tags$tr(tags$td(HTML("c")), tags$td("Relative growth of M compared with U. Values below 1 mean memory has a fitness cost."), tags$td("Relative growth of M slider.")),
                  tags$tr(tags$td(HTML("&alpha;")), tags$td("Stress penalty on U growth."), tags$td("Stress penalty on U growth slider.")),
                  tags$tr(tags$td(HTML("&mu;")), tags$td("Baseline symmetric switching between U and M."), tags$td("Symmetric baseline switching slider.")),
                  tags$tr(tags$td(HTML("&epsilon;")), tags$td("Memory decay from M back toward U."), tags$td("Memory decay slider.")),
                  tags$tr(tags$td(HTML("&beta;<sub>eff</sub>(t)")), tags$td("Stress-driven movement from U into M."), tags$td("Stress-driven switching plus stress-regime controls.")),
                  tags$tr(tags$td(HTML("&sigma;")), tags$td("Optional stochastic noise amplitude."), tags$td("Noise amplitude slider when stochastic mode is enabled."))
                )
              ),
              h5("Stress regimes"),
              tags$ul(
                tags$li(HTML("<b>No stress:</b> &beta;<sub>eff</sub>(t) = 0.")),
                tags$li(HTML("<b>Periodic stress:</b> stress alternates off/on every chosen half-cycle; &beta;<sub>eff</sub>(t) = &beta;iota(t).")),
                tags$li(HTML("<b>Continuous stress:</b> stress follows a smooth forcing curve and is transformed by a Hill response; &beta;<sub>eff</sub>(t) = &beta;&rho;(s)."))
              ),
              p(class = "small-muted", "Interpretation: the U/M model is a population-level model. It is useful for asking which combinations of switching, decay, stress, and fitness cost can reproduce the observed persistence of memory.")
            )
          )
        )
      ),
      nav_panel(
        "Formaldehyde equations",
        div(
          class = "m-3",
          card(
            class = "border-0 shadow-sm",
            card_body(
              h4("Mechanistic formaldehyde-resistance model"),
              p(HTML("This is the model used in the <b>Formaldehyde model</b> tab. It follows the slide-deck structure: formaldehyde enters the cell, resistance/enzyme expression rises, enzyme processing removes intracellular formaldehyde, and growth dilution reduces both resistance and intracellular formaldehyde.")),
              tags$div(class = "equation-box", HTML("<span class='eq-label'>Resistance/enzyme dynamics</span><span class='eq-main'>da/dt = f(d)H<sub>A</sub>(d)(1 + &eta;&xi;) - L(d)a</span>")),
              tags$div(class = "equation-box", HTML("<span class='eq-label'>Intracellular formaldehyde/damage dynamics</span><span class='eq-main'>dd/dt = K<sub>i</sub>(D<sub>D</sub>(t) - d) - K<sub>a</sub>ad/(k<sub>a</sub> + d) - L(d)d</span>")),
              h5("Functional forms implemented in this app"),
              tags$div(class = "equation-box", HTML("<span class='eq-label'>Growth / dilution</span><span class='eq-main'>L(d) = L<sub>max</sub> / (1 + (d/L<sub>50</sub>)<sup>n<sub>L</sub></sup>)</span>")),
              tags$div(class = "equation-box", HTML("<span class='eq-label'>Damage-induced expression</span><span class='eq-main'>H<sub>A</sub>(d) = basal + h<sub>max</sub>d<sup>n<sub>H</sub></sup> / (K<sub>expr</sub><sup>n<sub>H</sub></sup> + d<sup>n<sub>H</sub></sup>)</span>")),
              tags$div(class = "equation-box", HTML("<span class='eq-label'>Expression modulation</span><span class='eq-main'>f(d) = 1 + f<sub>amp</sub>d<sup>n<sub>f</sub></sup> / (f<sub>50</sub><sup>n<sub>f</sub></sup> + d<sup>n<sub>f</sub></sup>)</span>")),
              h5("Term-by-term mapping"),
              tags$table(
                class = "table table-sm table-striped doc-table",
                tags$thead(tags$tr(tags$th("Paper/model term"), tags$th("Meaning"), tags$th("Code/app meaning"))),
                tags$tbody(
                  tags$tr(tags$td(HTML("<b>a</b>")), tags$td("Resistance or enzyme level."), tags$td("Red resistance/enzyme trajectory and activation threshold.")),
                  tags$tr(tags$td(HTML("<b>d</b>")), tags$td("Intracellular formaldehyde or damage."), tags$td("Blue intracellular formaldehyde/damage trajectory.")),
                  tags$tr(tags$td(HTML("D<sub>D</sub>(t)")), tags$td("Extracellular formaldehyde concentration."), tags$td("Primary pulse and optional retest pulse.")),
                  tags$tr(tags$td(HTML("K<sub>i</sub>(D<sub>D</sub>(t)-d)")), tags$td("Diffusion into or out of the cell depending on the concentration difference."), tags$td("Diffusion curve in the model-process plot.")),
                  tags$tr(tags$td(HTML("K<sub>a</sub>ad/(k<sub>a</sub>+d)")), tags$td("Michaelis-Menten-like enzyme processing/removal of intracellular formaldehyde."), tags$td("Processing curve; increasing K_a should clear d faster.")),
                  tags$tr(tags$td(HTML("L(d)a and L(d)d")), tags$td("Dilution caused by cell growth."), tags$td("Growth/dilution curve; high d reduces growth.")),
                  tags$tr(tags$td(HTML("f(d)H<sub>A</sub>(d)")), tags$td("Damage-induced enzyme expression."), tags$td("Expression curve; increases as d activates resistance expression.")),
                  tags$tr(tags$td(HTML("1 + &eta;&xi;")), tags$td("Expression noise."), tags$td("Stochastic toggle and noise amplitude; creates lag variation across replicates."))
                )
              ),
              p(class = "small-muted", "Important: the equation structure is implemented directly, but exact numerical reproduction of the slide-deck figure requires the original parameter values. The app defaults are working simulation values for exploration, not fitted estimates.")
            )
          )
        )
      ),
      nav_panel(
        "Validation checklist",
        div(
          class = "m-3",
          card(
            class = "border-0 shadow-sm",
            card_body(
              h4("How to check that the website model matches the research model"),
              p("Use these checks after changing the code or parameters."),
              tags$table(
                class = "table table-sm table-striped doc-table",
                tags$thead(tags$tr(tags$th("Check"), tags$th("What to do"), tags$th("Expected result"))),
                tags$tbody(
                  tags$tr(tags$td("Equation audit"), tags$td("Download the formaldehyde simulation CSV and inspect expression, diffusion, processing, dilution, da_dt, and dd_dt."), tags$td(HTML("da_dt should equal expression - dilution_a. dd_dt should equal diffusion - processing - dilution_d."))),
                  tags$tr(tags$td("No formaldehyde"), tags$td("Set primary and retest extracellular formaldehyde to 0."), tags$td("d should stay low; a should not strongly activate except for basal expression.")),
                  tags$tr(tags$td("Pulse exposure"), tags$td("Apply a formaldehyde pulse."), tags$td("d rises first, then a rises, then d declines when processing becomes strong enough.")),
                  tags$tr(tags$td("No enzyme processing"), tags$td(HTML("Set K<sub>a</sub> = 0.")), tags$td("d should remain higher because the enzyme-removal term is disabled.")),
                  tags$tr(tags$td("Strong enzyme processing"), tags$td(HTML("Increase K<sub>a</sub>.")), tags$td("d should clear faster and activation lag should usually shorten.")),
                  tags$tr(tags$td("No stochasticity"), tags$td("Turn off stochastic expression noise."), tags$td("Replicate trajectories should become identical or nearly identical.")),
                  tags$tr(tags$td("With stochasticity"), tags$td("Turn stochastic expression noise back on."), tags$td("Replicates should show different lag times, matching the idea of lag-time variation."))
                )
              ),
              h5("Recommended code-level checks"),
              tags$ul(
                tags$li(HTML("Inside <code>simulate_formaldehyde()</code>, confirm that <code>da_det &lt;- expression - dilution_a</code>.")),
                tags$li(HTML("Inside <code>simulate_formaldehyde()</code>, confirm that <code>dd_det &lt;- diffusion - processing - dilution_d</code>.")),
                tags$li(HTML("Confirm that the output tibble includes <code>expression</code>, <code>diffusion</code>, <code>processing</code>, <code>dilution_a</code>, <code>dilution_d</code>, <code>da_dt</code>, and <code>dd_dt</code>.")),
                tags$li("Confirm that the plotted curves use the same variables that are exported in the CSV.")
              )
            )
          )
        )
      ),
      nav_panel(
        "Outputs & interpretation",
        div(
          class = "m-3",
          card(
            class = "border-0 shadow-sm",
            card_body(
              h4("How to interpret the website outputs"),
              tags$table(
                class = "table table-sm table-striped doc-table",
                tags$thead(tags$tr(tags$th("Output"), tags$th("What it means"), tags$th("Why it matters"))),
                tags$tbody(
                  tags$tr(tags$td("Population trajectories"), tags$td("Shows U, M, and total population over time."), tags$td("Useful for bulk memory dynamics and fitness-cost effects.")),
                  tags$tr(tags$td("Memory fraction"), tags$td(HTML("Shows M / (U + M).")), tags$td("Useful for measuring persistence or loss of memory.")),
                  tags$tr(tags$td("Stress forcing"), tags$td(HTML("Shows &beta;<sub>eff</sub>(t) or continuous stress intensity.")), tags$td("Useful for checking that the intended stress schedule is actually being applied.")),
                  tags$tr(tags$td("Formaldehyde concentration plot"), tags$td(HTML("Shows resistance/enzyme <b>a</b>, intracellular formaldehyde/damage <b>d</b>, and extracellular input D<sub>D</sub>(t).")), tags$td("This is the closest visual comparison to the formaldehyde slide-deck behavior.")),
                  tags$tr(tags$td("Model-process plot"), tags$td("Shows growth, expression, processing, diffusion, and dilution terms."), tags$td("Useful for debugging and explaining which mechanism is driving the trajectory.")),
                  tags$tr(tags$td("Lag distribution"), tags$td("Shows when stochastic replicates cross the resistance threshold."), tags$td("Useful for comparing the model to single-cell or microcolony lag-time data.")),
                  tags$tr(tags$td("Phase portraits"), tags$td("Shows the path through state space."), tags$td("Useful for seeing whether trajectories converge, cycle, or move through different regimes."))
                )
              ),
              h5("Best practice for research use"),
              tags$ol(
                tags$li("Start with deterministic settings to confirm the expected mechanism."),
                tags$li("Turn on stochasticity only after the deterministic behavior makes sense."),
                tags$li("Use the exported CSV to verify equations numerically."),
                tags$li("Use the compare tab for one-parameter-at-a-time sensitivity checks."),
                tags$li("Do not interpret default parameter values as biological estimates until they are fit to data.")
              )
            )
          )
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
  g <- ggplot() +
    make_theme() +
    labs(x = "Time (hours)", y = "Population fraction", color = NULL, fill = NULL)
  
  if (!is.null(rs)) {
    band_long <- bind_rows(
      rs %>% transmute(time, state = "U", lo = U_lo, hi = U_hi, med = U_med),
      rs %>% transmute(time, state = "M", lo = M_lo, hi = M_hi, med = M_med),
      rs %>% transmute(time, state = "total", lo = total_lo, hi = total_hi, med = total_med)
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
      geom_line(data = band_long, aes(time, med, color = state_label), linewidth = 1.15) +
      scale_fill_manual(values = setNames(c("#2563eb", "#db2777", "#ca8a04"), labels)) +
      scale_color_manual(values = setNames(c("#2563eb", "#db2777", "#ca8a04"), labels))
  } else {
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
    
    g <- g +
      geom_line(data = long, aes(time, value, color = state_label), linewidth = 1.15) +
      scale_color_manual(values = setNames(c("#2563eb", "#db2777", "#ca8a04"), labels))
  }
  
  g +
    {if (!is.null(view_time)) geom_vline(xintercept = view_time, linetype = "dashed", color = "#475569")} +
    coord_cartesian(xlim = c(0, max(dat$time)), ylim = c(0, 1))
}

plot_frac_gg <- function(dat, rs = NULL, view_time = NULL, eq = NA_real_,
                         labels = c(U = "U", M = "M")) {
  g <- ggplot() +
    make_theme() +
    labs(x = "Time (hours)", y = paste0(labels["M"], " / (", labels["U"], " + ", labels["M"], ")"))
  
  if (!is.null(rs)) {
    g <- g +
      geom_ribbon(data = rs, aes(time, ymin = mem_lo, ymax = mem_hi), fill = "#10b981", alpha = 0.16) +
      geom_line(data = rs, aes(time, mem_med), color = "#10b981", linewidth = 1.2) +
      geom_point(data = rs %>% slice_tail(n = 1), aes(time, mem_med), color = "#06b6d4", size = 2.8)
  } else {
    g <- g +
      geom_line(data = dat, aes(time, memory_fraction), color = "#10b981", linewidth = 1.2) +
      geom_point(data = dat %>% slice_tail(n = 1), aes(time, memory_fraction), color = "#06b6d4", size = 2.8)
  }
  
  g +
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

sync_slider_range <- function(session, input, suffix, slider_id,
                              default_min, default_max, step = NULL) {
  full_id <- paste0(slider_id, suffix)
  adv_id <- paste0("advanced_ranges", suffix)
  min_id <- paste0(slider_id, "_min", suffix)
  max_id <- paste0(slider_id, "_max", suffix)
  
  observe({
    adv <- isTRUE(input[[adv_id]])
    rr <- if (adv) {
      normalize_range(input[[min_id]], input[[max_id]], default_min, default_max)
    } else {
      list(min = default_min, max = default_max)
    }
    current <- input[[full_id]]
    if (is.null(current)) current <- rr$min
    current <- clamp(as.numeric(current), rr$min, rr$max)
    updateSliderInput(session, full_id, min = rr$min, max = rr$max, value = current, step = step)
  })
}

register_dynamic_ranges <- function(session, input, suffix = "") {
  sync_slider_range(session, input, suffix, "tmax", 24, 240, step = 12)
  sync_slider_range(session, input, suffix, "dt", 0.01, 0.50, step = 0.01)
  sync_slider_range(session, input, suffix, "U0", 0, 1, step = 0.01)
  sync_slider_range(session, input, suffix, "M0", 0, 1, step = 0.01)
  sync_slider_range(session, input, suffix, "r0", 0.05, 2.0, step = 0.01)
  sync_slider_range(session, input, suffix, "c_cost", 0.10, 1.20, step = 0.01)
  sync_slider_range(session, input, suffix, "alpha", 0, 1.0, step = 0.01)
  sync_slider_range(session, input, suffix, "mu", 0, 0.05, step = 0.001)
  sync_slider_range(session, input, suffix, "eps", 0, 0.20, step = 0.001)
  sync_slider_range(session, input, suffix, "beta", 0, 0.50, step = 0.005)
  sync_slider_range(session, input, suffix, "sigma", 0, 0.30, step = 0.005)
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
  
  register_dynamic_ranges(session, input, "")
  register_dynamic_ranges(session, input, "_a")
  register_dynamic_ranges(session, input, "_b")
  
  observeEvent(input$form_tmax, {
    current <- isolate(input$form_view_time %||% input$form_tmax)
    updateSliderInput(session, "form_view_time", max = input$form_tmax, value = min(current, input$form_tmax))
  }, ignoreInit = TRUE)
  
  observeEvent(input$form_pulse_start, {
    if (!is.null(input$form_pulse_start) && !is.null(input$form_retest_start) && input$form_retest_start < input$form_pulse_start) {
      updateSliderInput(session, "form_retest_start", value = input$form_pulse_start + input$form_pulse_duration)
    }
  }, ignoreInit = TRUE)
  
  form_params <- reactive(build_formaldehyde_params(input))
  
  form_data <- reactive({
    do.call(simulate_formaldehyde, form_params())
  })
  
  form_rep_data <- reactive({
    req(input$form_show_bands)
    pars <- form_params()
    if (!isTRUE(pars$stochastic)) return(NULL)
    do.call(simulate_formaldehyde_replicates, c(list(reps = input$form_replicates, seed_base = pars$seed), pars[names(pars) != "seed"]))
  })
  
  form_rep_summary <- reactive({
    dat <- form_rep_data()
    if (is.null(dat)) return(NULL)
    summarize_formaldehyde_replicates(dat)
  })
  
  form_lag_summary <- reactive({
    dat <- form_rep_data()
    if (is.null(dat)) return(NULL)
    summarize_formaldehyde_lags(dat)
  })
  
  form_current_row <- reactive({
    form_data() %>% filter(time <= input$form_view_time) %>% slice_tail(n = 1)
  })
  
  output$form_kpi_a <- renderText(sprintf("%.3f", form_current_row()$a))
  output$form_kpi_d <- renderText(sprintf("%.3f", form_current_row()$d))
  output$form_kpi_growth <- renderText(percent(form_current_row()$relative_growth, accuracy = 0.1))
  output$form_kpi_lag <- renderText({
    lg <- calc_activation_lag(form_data())
    if (is.na(lg)) "Not reached" else paste0(round(lg, 1), " min")
  })
  
  output$form_conc_plot <- renderPlot({
    dat <- form_data() %>% filter(time <= input$form_view_time)
    rs <- form_rep_summary()
    if (!is.null(rs)) rs <- rs %>% filter(time <= input$form_view_time)
    print(plot_formaldehyde_conc_gg(dat, rs, input$form_view_time))
  }, res = 110)
  
  output$form_conc_plot_overview <- renderPlot({
    dat <- form_data() %>% filter(time <= input$form_view_time)
    rs <- form_rep_summary()
    if (!is.null(rs)) rs <- rs %>% filter(time <= input$form_view_time)
    print(plot_formaldehyde_conc_gg(dat, rs, input$form_view_time))
  }, res = 100)
  
  output$form_process_plot <- renderPlot({
    dat <- form_data() %>% filter(time <= input$form_view_time)
    rs <- form_rep_summary()
    if (!is.null(rs)) rs <- rs %>% filter(time <= input$form_view_time)
    print(plot_formaldehyde_process_gg(dat, rs, input$form_view_time))
  }, res = 110)
  
  output$form_process_plot_overview <- renderPlot({
    dat <- form_data() %>% filter(time <= input$form_view_time)
    rs <- form_rep_summary()
    if (!is.null(rs)) rs <- rs %>% filter(time <= input$form_view_time)
    print(plot_formaldehyde_process_gg(dat, rs, input$form_view_time))
  }, res = 100)
  
  output$form_lag_plot <- renderPlot({
    lag_dat <- form_lag_summary()
    validate(need(!is.null(lag_dat) && nrow(lag_dat) > 0, "Turn on stochastic mode and replicate bands to view the lag distribution."))
    validate(need(any(!is.na(lag_dat$activation_lag)), "No replicate reached the resistance threshold. Try lowering the threshold or increasing expression."))
    print(plot_formaldehyde_lag_gg(lag_dat))
  }, res = 110)
  
  output$form_phase_plot <- renderPlot({
    dat <- form_data() %>% filter(time <= input$form_view_time)
    print(plot_formaldehyde_phase_gg(dat, input$form_view_time))
  }, res = 110)
  
  output$download_form_data <- downloadHandler(
    filename = function() paste0("formaldehyde_simulation_", Sys.Date(), ".csv"),
    content = function(file) write.csv(form_data(), file, row.names = FALSE)
  )
  
  output$download_form_reps <- downloadHandler(
    filename = function() paste0("formaldehyde_replicate_lags_", Sys.Date(), ".csv"),
    content = function(file) {
      lag_dat <- form_lag_summary()
      if (is.null(lag_dat)) lag_dat <- tibble(message = "Replicate lag summary unavailable unless stochastic mode and replicate bands are enabled.")
      write.csv(lag_dat, file, row.names = FALSE)
    }
  )
  
  output$download_form_plot <- downloadHandler(
    filename = function() paste0("formaldehyde_concentration_plot_", Sys.Date(), ".png"),
    content = function(file) {
      dat <- form_data()
      rs <- form_rep_summary()
      g <- plot_formaldehyde_conc_gg(dat, rs, max(dat$time))
      ggsave(file, g, width = 10, height = 6, dpi = 300, bg = "white")
    }
  )
  
  
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
    if (!is.null(rs)) rs <- rs %>% filter(time <= input$view_time)
    print(plot_traj_gg(dat, rs, input$view_time, labels = display_labels()))
  }, res = 110)
  
  output$traj_plot_overview <- renderPlot({
    dat <- sim_data() %>% filter(time <= input$view_time)
    rs <- rep_summary()
    if (!is.null(rs)) rs <- rs %>% filter(time <= input$view_time)
    print(plot_traj_gg(dat, rs, input$view_time, labels = display_labels()))
  }, res = 100)
  
  output$frac_plot <- renderPlot({
    dat <- sim_data() %>% filter(time <= input$view_time)
    rs <- rep_summary()
    if (!is.null(rs)) rs <- rs %>% filter(time <= input$view_time)
    eq <- calc_equilibrium(primary_params()$scenario, primary_params()$mu, primary_params()$eps, primary_params()$beta)
    print(plot_frac_gg(dat, rs, input$view_time, eq, labels = display_labels()))
  }, res = 110)
  
  output$frac_plot_overview <- renderPlot({
    dat <- sim_data() %>% filter(time <= input$view_time)
    rs <- rep_summary()
    if (!is.null(rs)) rs <- rs %>% filter(time <= input$view_time)
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
    dat <- bind_rows(sim_a() %>% mutate(run = "Run A"), sim_b() %>% mutate(run = "Run B"))
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
      mutate(state = recode(state, U = labels["U"], M = labels["M"]))
    
    g <- ggplot(dat, aes(time, value, color = run, linetype = state)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("Run A" = "#2563eb", "Run B" = "#f59e0b")) +
      make_theme() +
      labs(x = "Time (hours)", y = "Population fraction", color = NULL, linetype = NULL)
    print(g)
  }, res = 110)
  
  output$compare_stress_plot <- renderPlot({
    dat <- bind_rows(sim_a() %>% mutate(run = "Run A"), sim_b() %>% mutate(run = "Run B"))
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
          rs %>% transmute(time, state = "U", lo = U_lo, hi = U_hi, med = U_med),
          rs %>% transmute(time, state = "M", lo = M_lo, hi = M_hi, med = M_med),
          rs %>% transmute(time, state = "total", lo = total_lo, hi = total_hi, med = total_med)
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
          geom_line(data = band_long, aes(time, med, color = state_label), linewidth = 1.15) +
          scale_fill_manual(values = setNames(c("#60a5fa", "#f472b6", "#facc15"), labels)) +
          scale_color_manual(values = setNames(c("#60a5fa", "#f472b6", "#facc15"), labels))
      } else {
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
        
        g <- g +
          geom_line(data = long, aes(time, value, color = state_label), linewidth = 1.15) +
          scale_color_manual(values = setNames(c("#60a5fa", "#f472b6", "#facc15"), labels))
      }
      
      g <- g +
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
        g <- g +
          geom_ribbon(data = rs, aes(time, ymin = mem_lo, ymax = mem_hi), fill = "#34d399", alpha = 0.16) +
          geom_line(data = rs, aes(time, mem_med), color = "#34d399", linewidth = 1.2)
      } else {
        g <- g +
          geom_line(data = dat, aes(time, memory_fraction), color = "#34d399", linewidth = 1.2)
      }
      
      g <- g +
        coord_cartesian(xlim = c(0, max(dat$time)), ylim = c(0, 1)) +
        labs(x = "Time (hours)", y = paste0(labels["M"], " / (", labels["U"], " + ", labels["M"], ")"))
      
      ggsave(file, g, width = 9, height = 5.5, dpi = 300, bg = "#111827")
    }
  )
}

shinyApp(ui, server)
