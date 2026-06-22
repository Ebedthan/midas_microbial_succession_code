# =============================================================================
# Microbial secondary succession West African post-forest landscapes
# Author: Anicet Ebou, ediman dot ebou at inphb dot ci
# Tested on Ubuntu Linux 22.04
# My personal lib for this project
# =============================================================================

# loading my custom fns
source(here::here("02_R/mylib.R"))

dir.create("03_models", showWarnings = FALSE)
dir.create("04_figures", showWarnings = FALSE)
dir.create("05_results", showWarnings = FALSE)

# loading diversity data
data <- read.csv(here::here("01_data/div.csv"), row.names = 1) |>
  dplyr::mutate(
    is_ogf = ifelse(categorie == "OGF", 1, 0),
    age = ifelse(is_ogf == 1, 0, age)
    )

# a simple fn to fit models
fit_div_model <- function(guild_name, df) {
  d <- df |> dplyr::filter(guild == guild_name)

  bf_joint <- brms::bf(
    richness ~ is_ogf * thinf +
      (1 - is_ogf) * (th0 + (thinf - th0) * (1 - exp(-lam * age))),
    th0   ~ 1,
    thinf ~ 1,
    lam ~ 1,
    nl = TRUE
  )

  priors <- c(
    brms::prior(normal(600, 200), nlpar = "th0", lb = 0),
    brms::prior(normal(800, 200), nlpar = "thinf", lb = 0),
    brms::prior(normal(0.1, 0.1), nlpar = "lam", lb = 0),
    brms::prior(normal(0, 1), class = "sigma")
  )

  brms::brm(
    formula = bf_joint,
    data = d,
    family = lognormal(),
    prior = priors,
    chains  = 4, iter = 4000, warmup = 2000, cores = 4, seed = 42,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    file = paste0("03_models/fit_joint_", guild_name)
  )
}

fit_joint_bact <- fit_div_model("Bacteria", data)
fit_joint_amf <- fit_div_model("AMF", data)

# extract raw posterior draws
extract_posteriors <- function(fit, guild_name) {
  brms::as_draws_df(fit) |>
    dplyr::transmute(
      guild = guild_name,
      # Parameters are on log scale inside lognormal
      th0 = exp(b_th0_Intercept),
      thinf = exp(b_thinf_Intercept),
      lam = exp(b_lam_Intercept),
      sigma = sigma,
      # Derived: recovery half-life (time to reach 50% of the way to thinf)
      half_life = log(2) / lam
    )
}

post_bact <- extract_posteriors(fit_joint_bact, "Bacteria")
post_amf <- extract_posteriors(fit_joint_amf, "AMF")
post_all <- dplyr::bind_rows(post_bact, post_amf)


# Summary table
param_summary <- post_all |>
  dplyr::group_by(guild) |>
  dplyr::summarise(
    dplyr::across(
      c(th0, thinf, lam, sigma, half_life),
      list(
        mean = mean,
        sd = sd,
        lo = ~ quantile(.x, 0.025),
        hi = ~ quantile(.x, 0.975)
      ),
      .names = "{.col}__{.fn}"
    )
  ) |>
  tidyr::pivot_longer(
    cols = -guild,
    names_to = c("parameter", "stat"),
    names_sep = "__"
  ) |>
  tidyr::pivot_wider(names_from = stat, values_from = value) |>
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

param_summary

# Formatted summary
param_table <- param_summary |>
  dplyr::mutate(
    label = paste0(round(mean, 2), " [", round(lo, 2), ", ", round(hi, 2), "]")
  ) |>
  dplyr::select(guild, parameter, label) |>
  tidyr::pivot_wider(names_from = guild, values_from = label)

param_table

# Posterior overlap between guilds
compare_guilds <- post_all |>
  dplyr::select(guild, th0, thinf, lam, half_life) |>
  tidyr::pivot_longer(cols = -guild, names_to = "parameter") |>
  dplyr::group_by(parameter) |>
  dplyr::summarise(
    # Probability that Bacteria value > AMF value
    p_bact_greater = mean(
      value[guild == "Bacteria"] > value[guild == "AMF"]
    ),
    # Difference (Bacteria - AMF): mean and 95% CI
    diff_mean = mean(value[guild == "Bacteria"] - value[guild == "AMF"]),
    diff_lo = quantile(value[guild == "Bacteria"] - value[guild == "AMF"], 0.025),
    diff_hi = quantile(value[guild == "Bacteria"] - value[guild == "AMF"], 0.975)
  ) |>
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

compare_guilds
