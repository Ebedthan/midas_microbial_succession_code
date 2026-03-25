# =============================================================================
# Microbial secondary succession West African post-forest landscapes
# Author: Anicet Ebou, ediman dot ebou at inphb dot ci
# Tested on Ubuntu Linux 22.04
# Alpha diversity models: richness, Shannon, Simpson
# Confounders included: crop_type, crop_year, topography, soil_type
# =============================================================================

source(here::here("02_R/mylib.R"))

dir.create("03_models", showWarnings = FALSE)
dir.create("04_figures", showWarnings = FALSE)
dir.create("05_results", showWarnings = FALSE)

# Step 0. Load data
data <- read.csv(here::here("01_data/data.csv"))

# Step 1. Confounder balance check
# Before modelling, verify that confounders are not perfectly collinear with age

cat("\n=== Confounder balance across age classes ===\n")
cat("\n-- crop_type × age category --\n")
print(table(data$crop_type, data$categorie, useNA = "ifany"))

cat("\n-- soil_type × age category --\n")
print(table(data$soil_type, data$categorie, useNA = "ifany"))

cat("\n-- topography × age category --\n")
print(table(data$topography, data$categorie, useNA = "ifany"))

cat("\n-- crop_year × age category (mean per class) --\n")
data |>
  dplyr::filter(categorie != "OGF") |>
  dplyr::group_by(categorie) |>
  dplyr::summarise(
    mean_crop_year = mean(crop_year, na.rm = TRUE),
    sd_crop_year = sd(crop_year,   na.rm = TRUE),
    n = dplyr::n()
  ) |>
  print()

# Step 2. Exploratory plots
data.div <- data |>
  dplyr::select(
    categorie, crop_type, crop_year, topography, soil_type,
    dplyr::starts_with("amf"),
    dplyr::starts_with("bacterial")
  ) |>
  tidyr::pivot_longer(
    cols = amf_richness:bacterial_simpson,
    names_to = "estimate",
    values_to = "value"
  )

facet.labels <- c(
  amf_richness = "AMF richness",
  amf_shannon = "AMF Shannon",
  amf_simpson = "AMF Simpson",
  bacterial_richness = "Bacterial richness",
  bacterial_shannon = "Bacterial Shannon",
  bacterial_simpson = "Bacterial Simpson"
)

p1 <- ggplot2::ggplot(data.div, ggplot2::aes(x = categorie, y = value)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.15, seed = 42),
    alpha = 0.6,
    size = 1.8
  ) +
  ggplot2::facet_wrap(
    ~ estimate, scales = "free",
    labeller = ggplot2::labeller(estimate = facet.labels)
  ) +
  ggplot2::labs(x = "Age category", y = NULL) +
  theme_Publication()

export_plot(p1, here::here("04_figures/Fig_exploratory.png"), width = 13, height = 8)

# Step 3. Build long modelling dataset
make_long <- function(data, metric) {
  data |>
    dplyr::select(
      site = forest,
      sample, categorie, age,
      crop_type, crop_year,
      topography, soil_type,
      bact_val = paste0("bacterial_", metric),
      amf_val = paste0("amf_", metric)
    ) |>
    tidyr::pivot_longer(
      cols = c(bact_val, amf_val),
      names_to = "guild",
      values_to = "diversity"
    ) |>
    dplyr::mutate(
      guild = factor(
        guild,
        levels = c("bact_val", "amf_val"),
        labels = c("Bacteria", "AMF")
      ),
      is_ogf = ifelse(categorie == "OGF", 1, 0),
      age_model = ifelse(is_ogf == 1, NA, age),
      
      # crop_year: use 0 for OGF (no agricultural legacy)
      crop_year_model = ifelse(is_ogf == 1, 0, crop_year),
      
      # crop_type: recode OGF as explicit "none" level
      crop_type = dplyr::case_when(
        is_ogf == 1 ~ "cacao", # used as reference level as most common in data
        TRUE ~ crop_type
      ),
      crop_type = factor(crop_type, levels = c("cacao", "cacao_cafe", "cafe", "coton", "igname", "mais", "riz")),
      topography = factor(topography),
      soil_type = factor(soil_type),
      metric = metric
    )
}

long.richness <- make_long(data, "richness")
long.shannon <- make_long(data, "shannon")
long.simpson <- make_long(data, "simpson")

# Step 4. Fit function
fit_model <- function(long_data, metric, family) {
  dir.create(here::here("03_models"), showWarnings = FALSE)
  d <- long_data |>
    dplyr::filter(!is.na(age_model) | is_ogf == 1)
  
  if (family == "negbinomial") {
    priors <- c(
      brms::prior(normal(0, 1), class = "b"),
      brms::prior(lognormal(-1, 0.5), class = "sds"),
      brms::prior(lognormal(-1, 0.5), class = "sd"),
      brms::prior(gamma(0.01, 0.01), class = "shape")
    )
    fam <- brms::negbinomial()
  } else {
    priors <- c(
      brms::prior(normal(0, 1), class = "b"),
      brms::prior(lognormal(-1, 0.5), class = "sds"),
      brms::prior(lognormal(-1, 0.5), class = "sd"),
      brms::prior(lognormal(-1, 0.5), class = "sigma")
    )
    fam <- brms::lognormal()
  }
  
  brms::brm(
    diversity ~ s(age_model, by = guild) +
      guild +
      is_ogf +
      crop_year_model +
      crop_type +
      topography +
      soil_type +
      (1 | site),
    data = d,
    family = fam,
    prior = priors,
    chains = 4,
    cores = 4,
    iter = 4000,
    warmup = 2000,
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    seed = 42,
    file = file.path(here::here("03_models"), paste0("fit_confounders_", metric))
  )
}

cat("\n>>> Fitting richness model (negative binomial)...\n")
fit_richness <- fit_model(long.richness, "richness", "negbinomial")

cat("\n>>> Fitting Shannon model (lognormal)...\n")
fit_shannon  <- fit_model(long.shannon,  "shannon",  "lognormal")

cat("\n>>> Fitting Simpson model (lognormal)...\n")
fit_simpson  <- fit_model(long.simpson,  "simpson",  "lognormal")

# Step 5. Parameter extraction
extract_params <- function(fit, metric_label, family) {
  post <- brms::as_draws_df(fit)
  # Fixed effects
  fixed <- post |>
    dplyr::select(
      dplyr::starts_with("b_"),
      dplyr::any_of(c("sigma", "shape"))
    ) |>
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to = "parameter",
      values_to = "value"
    ) |>
    dplyr::group_by(parameter) |>
    dplyr::summarise(
      mean = mean(value),
      median = median(value),
      sd = sd(value),
      lci = quantile(value, 0.025),
      uci = quantile(value, 0.975),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      component = dplyr::case_when(
        grepl("^b_crop_type", parameter) ~ "Confounder: crop type",
        grepl("^b_topography", parameter) ~ "Confounder: topography",
        grepl("^b_soil_type", parameter) ~ "Confounder: soil type",
        grepl("crop_year", parameter) ~ "Confounder: crop year",
        TRUE ~ "Fixed effects"
      ),
      parameter = dplyr::case_when(
        parameter == "b_Intercept"  ~ "Intercept (Bacteria, log-scale)",
        parameter == "b_guildAMF" ~ "Guild: AMF vs Bacteria",
        parameter == "b_is_ogf" ~ "Old-growth forest effect",
        parameter == "b_crop_year_model"~ "Crop year (continuous)",
        parameter == "sigma" ~ "Residual SD (σ)",
        parameter == "shape" ~ "Neg. binomial shape",
        # crop_type levels
        grepl("^b_crop_type", parameter) ~ paste0("Crop: ", gsub("b_crop_type", "", parameter)),
        # topography levels
        grepl("^b_topography", parameter) ~ paste0("Topography: ", gsub("b_topography", "", parameter)),
        # soil_type levels
        grepl("^b_soil_type", parameter)~ paste0("Soil: ", gsub("b_soil_type", "", parameter)),
        TRUE ~ parameter
      )
    )
  
  # Smooth term SDs
  smooth <- post |>
    dplyr::select(dplyr::starts_with("sds_")) |>
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to = "parameter",
      values_to = "value"
    ) |>
    dplyr::group_by(parameter) |>
    dplyr::summarise(
      mean = mean(value),
      median = median(value),
      sd = sd(value),
      lci = quantile(value, 0.025),
      uci = quantile(value, 0.975),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      component = "Smooth terms",
      parameter = dplyr::case_when(
        grepl("Bacteria", parameter) ~ "s(age) - Bacteria",
        grepl("AMF", parameter) ~ "s(age) - AMF",
        TRUE ~ parameter
      )
    )
  
  # Random effects
  random <- post |>
    dplyr::select(dplyr::starts_with("sd_")) |>
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to = "parameter",
      values_to = "value"
    ) |>
    dplyr::group_by(parameter) |>
    dplyr::summarise(
      mean = mean(value),
      median = median(value),
      sd = sd(value),
      lci = quantile(value, 0.025),
      uci = quantile(value, 0.975),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      component = "Random effects",
      parameter = "Site intercept SD"
    )
  
  dplyr::bind_rows(fixed, smooth, random) |>
    dplyr::mutate(metric = metric_label)
}

# Step 6. Posterior probabilities

compute_posteriors <- function(fit, long_data, metric_label) {
  sf_data <- long_data |> dplyr::filter(is_ogf == 0)
  modal_crop <- names(sort(table(sf_data$crop_type), decreasing = TRUE))[1]
  modal_topo <- names(sort(table(sf_data$topography), decreasing = TRUE))[1]
  modal_soil <- names(sort(table(sf_data$soil_type),  decreasing = TRUE))[1]
  mean_cyear <- mean(sf_data$crop_year_model, na.rm = TRUE)
  
  cat(sprintf(
    "\n[%s] Marginalising over: crop_type=%s, topography=%s, soil_type=%s, crop_year=%.1f\n",
    metric_label, modal_crop, modal_topo, modal_soil, mean_cyear
  ))
  
  # Secondary forest newdata: all confounders at modal/mean values
  make_nd_sf <- function(age_val, guild_val) {
    data.frame(
      age_model = age_val,
      guild = guild_val,
      is_ogf = 0,
      crop_year_model = mean_cyear,
      crop_type = modal_crop,
      topography = modal_topo,
      soil_type = modal_soil,
      site = NA
    )
  }
  
  # OGF newdata: crop_type and crop_year_model set to their OGF training values
  # Use the most common topography/soil_type seen in OGF plots specifically
  ogf_data <- long_data |> dplyr::filter(is_ogf == 1)
  modal_topo_ogf <- names(sort(table(ogf_data$topography), decreasing = TRUE))[1]
  modal_soil_ogf <- names(sort(table(ogf_data$soil_type),  decreasing = TRUE))[1]
  
  make_nd_ogf <- function(age_val, guild_val) {
    data.frame(
      age_model = age_val,
      guild = guild_val,
      is_ogf = 1,
      crop_year_model = 0, # OGF: no agricultural legacy
      crop_type = "cacao", # valid level seen during training
      topography = modal_topo_ogf,
      soil_type = modal_soil_ogf,
      site = NA
    )
  }
  
  mean_age <- long_data |>
    dplyr::filter(!is.na(age_model)) |>
    dplyr::pull(age_model) |>
    mean()
  
  p_yb <- brms::posterior_epred(fit, newdata = make_nd_sf(5, "Bacteria"), re_formula = NA)
  p_ob <- brms::posterior_epred(fit, newdata = make_nd_sf(40, "Bacteria"), re_formula = NA)
  p_ya <- brms::posterior_epred(fit, newdata = make_nd_sf(5, "AMF"), re_formula = NA)
  p_oa <- brms::posterior_epred(fit, newdata = make_nd_sf(40, "AMF"), re_formula = NA)
  p_ogf_b <- brms::posterior_epred(fit, newdata = make_nd_ogf(mean_age, "Bacteria"), re_formula = NA)
  p_ogf_a <- brms::posterior_epred(fit, newdata = make_nd_ogf(mean_age, "AMF"), re_formula = NA)
  
  tibble::tibble(
    metric = metric_label,
    quantity = c(
      "P(Bacteria increases, age 5 => 40)",
      "P(AMF increases, age 5 => 40)",
      "P(AMF trend > Bacteria trend)",
      "P(Bacteria at 40 yr ≥ OGF level)",
      "P(AMF at 40 yr ≥ OGF level)"
    ),
    estimate = c(
      mean(p_ob > p_yb),
      mean(p_oa > p_ya),
      mean((p_oa - p_ya) > (p_ob - p_yb)),
      mean(p_ob >= p_ogf_b),
      mean(p_oa >= p_ogf_a)
    )
  )
}

params_richness <- extract_params(fit_richness, "Richness", "negbinomial")
params_shannon <- extract_params(fit_shannon, "Shannon", "lognormal")
params_simpson <- extract_params(fit_simpson, "Simpson", "lognormal")

post_richness <- compute_posteriors(fit_richness, long.richness, "Richness")
post_shannon <- compute_posteriors(fit_shannon, long.shannon, "Shannon")
post_simpson <- compute_posteriors(fit_simpson, long.simpson, "Simpson")

all_params <- dplyr::bind_rows(params_richness, params_shannon, params_simpson) |>
  dplyr::mutate(
    dplyr::across(c(mean, median, sd, lci, uci), ~ round(.x, 3)),
    ci95 = paste0("[", lci, ", ", uci, "]"),
    metric = factor(metric, levels = c("Richness", "Shannon", "Simpson"))
  ) |>
  dplyr::select(metric, component, parameter, mean, median, sd, ci95)

all_posteriors <- dplyr::bind_rows(post_richness, post_shannon, post_simpson) |>
  dplyr::mutate(
    estimate = round(estimate, 3),
    metric = factor(metric, levels = c("Richness", "Shannon", "Simpson"))
  )

print(all_params, n = Inf)
print(all_posteriors, n = Inf)

write.csv(all_params, here::here("05_results/params_confounders.csv"), row.names = FALSE)
write.csv(all_posteriors, here::here("05_results/posteriors_confounders.csv"), row.names = FALSE)

sf_data <- long.richness |> dplyr::filter(is_ogf == 0)
modal_crop <- names(sort(table(sf_data$crop_type), decreasing = TRUE))[1]
modal_topo <- names(sort(table(sf_data$topography), decreasing = TRUE))[1]
modal_soil <- names(sort(table(sf_data$soil_type),  decreasing = TRUE))[1]
mean_cyear <- mean(sf_data$crop_year_model, na.rm = TRUE)

newdat.base <- expand.grid(
  age_model = seq(
    min(long.richness$age_model, na.rm = TRUE),
    max(long.richness$age_model, na.rm = TRUE),
    length.out = 100
  ),
  guild = levels(long.richness$guild),
  is_ogf = 0,
  crop_type = modal_crop,
  crop_year_model = mean_cyear,
  topography = modal_topo,
  soil_type = modal_soil,
  site = NA
)

# Get posterior predictions from each model, tag with metric
get_preds <- function(fit, newdat, metric_label) {
  tidybayes::add_epred_draws(
    fit,
    newdata = newdat,
    re_formula = NA
  ) |>
    dplyr::mutate(metric = metric_label)
}

pred.richness <- get_preds(fit_richness, newdat.base, "Richness")
pred.shannon <- get_preds(fit_shannon, newdat.base, "Shannon")
pred.simpson <- get_preds(fit_simpson, newdat.base, "Simpson")

pred.all <- dplyr::bind_rows(pred.richness, pred.shannon, pred.simpson) |>
  dplyr::mutate(
    metric = factor(metric, levels = c("Richness", "Shannon", "Simpson")),
    guild = factor(guild,  levels = c("Bacteria", "AMF"))
  )

# Observed data: bind all three long datasets
obs.all <- dplyr::bind_rows(
  long.richness |> dplyr::mutate(metric = "Richness"),
  long.shannon  |> dplyr::mutate(metric = "Shannon"),
  long.simpson  |> dplyr::mutate(metric = "Simpson")
) |>
  dplyr::filter(is_ogf == 0) |>   # secondary forest only on the curve panel
  dplyr::mutate(
    metric = factor(metric, levels = c("Richness", "Shannon", "Simpson")),
    guild = factor(guild,  levels = c("Bacteria", "AMF"))
  )

# OGF reference points
ogf.all <- dplyr::bind_rows(
  long.richness |> dplyr::mutate(metric = "Richness"),
  long.shannon |> dplyr::mutate(metric = "Shannon"),
  long.simpson |> dplyr::mutate(metric = "Simpson")
) |>
  dplyr::filter(is_ogf == 1) |>
  dplyr::mutate(
    metric = factor(metric, levels = c("Richness", "Shannon", "Simpson")),
    guild = factor(guild, levels = c("Bacteria", "AMF")),
    # Place OGF points just beyond the age axis for visual separation
    age_plot = max(long.richness$age_model, na.rm = TRUE) * 1.12
  )
guild_colours <- c("Bacteria" = "#2166ac", "AMF" = "#d6604d")
guild_fills <- c("Bacteria" = "#2166ac", "AMF" = "#d6604d")

p_div <- ggplot2::ggplot(
  pred.all,
  ggplot2::aes(x = age_model, y = .epred, colour = guild, fill = guild)
) +
  tidybayes::stat_lineribbon(
    .width = 0.95,
    alpha = 0.20,
    linewidth = 0.9
  ) +
  # Observed secondary forest points
  ggplot2::geom_point(
    data = obs.all,
    mapping = ggplot2::aes(x = age_model, y = diversity, colour = guild),
    shape = 16,
    size = 1.8,
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  # OGF observed points as triangles at right margin
  ggplot2::geom_point(
    data = ogf.all,
    mapping = ggplot2::aes(x = age_plot, y = diversity, colour = guild),
    shape = 17,
    size = 2.2,
    alpha = 0.8,
    inherit.aes = FALSE
  ) +
  # Vertical dashed separator before OGF strip
  ggplot2::geom_vline(
    xintercept = max(long.richness$age_model, na.rm = TRUE) * 1.06,
    linetype = "dashed",
    colour = "grey60",
    linewidth = 0.4
  ) +
  # OGF label on top of first facet only
  ggplot2::annotate(
    "text",
    x = max(long.richness$age_model, na.rm = TRUE) * 1.12,
    y = Inf,
    label = "OGF",
    colour = "grey40",
    size = 3,
    vjust = 1.5,
    hjust = 0.5
  ) +
  ggplot2::scale_colour_manual(
    values = guild_colours,
    labels = c("Bacteria", "AM fungi")
  ) +
  ggplot2::scale_fill_manual(
    values = guild_fills,
    labels = c("Bacteria", "AM fungi")
  ) +
  ggplot2::facet_wrap(
    guild ~ metric,
    scales = "free_y",
    ncol = 3
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(0, 10, 20, 30, 40),
    expand = ggplot2::expansion(mult = c(0.02, 0.18))
  ) +
  ggplot2::labs(
    x = "Time since abandonment (years)",
    y = "Alpha diversity",
    colour = "Guild",
    fill = "Guild",
    caption = "Shaded band: 95% credible interval | Circles: secondary forest plots | Triangles: old-growth reference plots"
  ) +
  theme_Publication() +
  ggplot2::theme(
    legend.position = "bottom",
    strip.text.y = ggplot2::element_text(face = "bold", size = 11),
    plot.caption = ggplot2::element_text(size = 8, colour = "grey50"),
    panel.spacing = ggplot2::unit(0.8, "lines")
  )
export_plot(p_div, here::here("04_figures/Fig_alpha_diversity.png"), width = 13, height = 8)
