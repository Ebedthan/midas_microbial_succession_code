# =============================================================================
# Microbial secondary succession West African post-forest landscapes
# Author: Anicet Ebou, ediman dot ebou at inphb dot ci
# Tested on Ubuntu Linux 22.04
# Beta diversity analysis:
#   - 1. Pairwise dissimilarity-based beta regression (turnover)
#   - 2. Within-site dispersion along chronosequence (stability)
#   - 3. Distance to OGF (compositional recovery) with confounders
# =============================================================================

source(here::here("02_R/mylib.R"))

dir.create("03_models", showWarnings = FALSE)
dir.create("04_figures", showWarnings = FALSE)
dir.create("05_results", showWarnings = FALSE)

# Step 0. Load data and ASV tables
data <- read.csv("data.csv")
asv.bact <- read.csv("asv_table.csv", row.names = 1, check.names = FALSE)
asv.amf <- read.csv("amf_table.csv", row.names = 1, check.names = FALSE)

asv.bact <- t(asv.bact)
asv.amf <- t(asv.amf)

asv.bact.rel <- asv.bact / rowSums(asv.bact)
asv.amf.rel <- asv.amf  / rowSums(asv.amf)
asv.bact.pa <- asv.bact > 0
asv.amf.pa <- asv.amf  > 0

# Step 1. Compute all dissimilarity matrices
dist.bact.sorensen <- vegan::vegdist(asv.bact.pa, method = "bray")
dist.amf.sorensen <- vegan::vegdist(asv.amf.pa, method = "bray")
dist.bact.horn <- vegan::vegdist(asv.bact.rel, method = "horn")
dist.amf.horn <- vegan::vegdist(asv.amf.rel, method = "horn")
dist.bact.mh <- vegan::vegdist(asv.bact, method = "morisita")
dist.amf.mh <- vegan::vegdist(asv.amf, method = "morisita")
dist.bact.bray <- vegan::vegdist(asv.bact.rel, method = "bray")
dist.amf.bray <- vegan::vegdist(asv.amf.rel, method = "bray")

metrics <- c("Sorensen", "Horn", "Morisita-Horn", "Bray-Curtis")
guilds <- c("Bacteria", "AMF")

guild_colours <- c("Bacteria" = "#2166ac", "AMF" = "#d6604d")
guild_fills <- c("Bacteria" = "#2166ac", "AMF" = "#d6604d")

# Helpers 
## retrieve the right distance matrix
get_dist <- function(guild, metric) {
  key <- paste(guild, metric, sep = "_")
  list(
    "Bacteria_Sorensen" = dist.bact.sorensen,
    "Bacteria_Horn" = dist.bact.horn,
    "Bacteria_Morisita-Horn" = dist.bact.mh,
    "Bacteria_Bray-Curtis" = dist.bact.bray,
    "AMF_Sorensen" = dist.amf.sorensen,
    "AMF_Horn" = dist.amf.horn,
    "AMF_Morisita-Horn" = dist.amf.mh,
    "AMF_Bray-Curtis" = dist.amf.bray
  )[[key]]
}

## safe model file key
model_key <- function(metric, guild) {
  paste0(tolower(gsub("[- ]", "", metric)), "_", tolower(guild))
}

# Step 2. NMDS ordination
nmds.bact <- vegan::metaMDS(dist.bact.horn, k = 2, trymax = 100)
nmds.amf  <- vegan::metaMDS(dist.amf.horn, k = 2, trymax = 100)

cat("Bacteria NMDS stress:", nmds.bact$stress, "\n")
cat("AMF NMDS stress:", nmds.amf$stress, "\n")

rownames(data) <- data$sample
bact_scores <- as.data.frame(nmds.bact$points)
amf_scores  <- as.data.frame(nmds.amf$points)

stopifnot(all(rownames(bact_scores) %in% rownames(data)))
stopifnot(all(rownames(amf_scores)  %in% rownames(data)))

## Align metadata to score row order
meta_bact <- data[rownames(bact_scores), ]
meta_amf  <- data[rownames(amf_scores), ]

max_sf_age <- 100
display_age_bact <- ifelse(is.na(meta_bact$age), 100, meta_bact$age)
display_age_amf <- ifelse(is.na(meta_amf$age), 100, meta_amf$age)

# flag OGF plots for optional shape distinction
is_ogf_bact <- is.na(meta_bact$age)
is_ogf_amf <- is.na(meta_amf$age)

nmds.df <- dplyr::bind_rows(
  bact_scores |>
    dplyr::mutate(
      guild = "Bacteria",
      forest = meta_bact$forest,
      age = display_age_bact,
      is_ogf = is_ogf_bact,
      sample = rownames(bact_scores)
    ),
  amf_scores |>
    dplyr::mutate(
      guild = "AMF",
      forest = meta_amf$forest,
      age = display_age_amf,
      is_ogf = is_ogf_amf,
      sample = rownames(amf_scores)
    )
)

p_nmds <- ggplot2::ggplot(
  nmds.df,
  ggplot2::aes(
    x = MDS1,
    y = MDS2,
    colour = forest,
    size = age,
    shape = is_ogf
  )
) +
  ggplot2::geom_point(alpha = 0.8) +
  ggplot2::scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("FALSE" = "Secondary forest", "TRUE" = "Old-growth forest"),
    name = NULL
  ) +
  ggsci::scale_color_aaas() +
  ggplot2::scale_size_continuous(
    name = "Successional age",
    breaks = c(10, 20, 30, 40, max_sf_age * 1.2),
    labels = c("10", "20", "30", "40", "OGF")
  ) +
  ggplot2::facet_wrap(~guild) +
  ggplot2::labs(
    x = "NMDS1",
    y = "NMDS2",
    colour = "Forest"
  ) +
  theme_Publication()
export_plot(p_nmds, here::here("04_figures/Fig_NMDS.png"), width = 13, height = 8)

# Variance partitioning
sf_idx <- which(data$categorie != "OGF")
sf_data <- data[sf_idx, ]

# Hellinger-transformed community matrices for SF plots only
bact_hell_sf <- vegan::decostand(asv.bact[sf_idx, ], method = "hellinger")
amf_hell_sf  <- vegan::decostand(asv.amf[sf_idx, ], method = "hellinger")

# Predictor matrices
X_age <- data.frame(age = sf_data$age)
X_zone <- data.frame(zone = sf_data$zone)
X_edaph <- model.matrix(~ soil_type + topography - 1, data = sf_data)
X_site <- model.matrix(~ forest - 1, data = sf_data)

# bacteria VP
vp_bact <- vegan::varpart(bact_hell_sf, X_age, X_zone, X_edaph, X_site)

# AMF VP
vp_amf <- vegan::varpart(amf_hell_sf, X_age, X_zone, X_edaph, X_site)

# Test significance of each fraction via RDA
rda_age_bact <- vegan::rda(
  bact_hell_sf ~ age + Condition(as.matrix(X_zone)) +
    Condition(X_edaph) + Condition(X_site), data = sf_data
  )
rda_age_amf <- vegan::rda(
  amf_hell_sf  ~ age + Condition(as.matrix(X_zone)) +
    Condition(X_edaph) + Condition(X_site), data = sf_data
  )

anova_age_bact <- vegan::anova.cca(rda_age_bact, permutations = 9999)
anova_age_amf <- vegan::anova.cca(rda_age_amf, permutations = 9999)

# Save results
saveRDS(list(
  vp_bact = vp_bact, 
  vp_amf = vp_amf,
  anova_age_bact = anova_age_bact,
  anova_age_amf  = anova_age_amf
  ), here::here("03_models/varpart_results.rds"))

# Step 3. Build pairwise dissimilarity dataframe
make_pairs <- function(dist_mat, data, guild_name, metric_name) {
  dist_mat <- as.matrix(dist_mat)
  sf_idx <- which(data$categorie != "OGF")
  sf_data <- data[sf_idx, ]
  sf_mat <- dist_mat[sf_idx, sf_idx]
  pairs <- which(lower.tri(sf_mat), arr.ind = TRUE)
  
  tibble::tibble(
    guild = guild_name,
    metric = metric_name,
    forest_i = sf_data$forest[pairs[,1]],
    forest_j = sf_data$forest[pairs[,2]],
    age_i = sf_data$age[pairs[,1]],
    age_j = sf_data$age[pairs[,2]],
    crop_type_i = sf_data$crop_type[pairs[,1]],
    crop_type_j = sf_data$crop_type[pairs[,2]],
    soil_type_i = sf_data$soil_type[pairs[,1]],
    soil_type_j = sf_data$soil_type[pairs[,2]],
    topo_i = sf_data$topography[pairs[,1]],
    topo_j = sf_data$topography[pairs[,2]],
    dissimilarity = sf_mat[pairs]
  ) |>
    dplyr::mutate(
      age_diff = abs(age_i - age_j),
      same_crop = crop_type_i == crop_type_j,
      same_soil = soil_type_i == soil_type_j,
      same_topo = topo_i == topo_j,
      # scale age_diff to [0,1] for numerical stability
      age_diff_s = age_diff / max(age_diff, na.rm = TRUE),
      dissimilarity = pmin(pmax(dissimilarity, 1e-4), 1 - 1e-4)
    )
}

pairs.all <- purrr::map_dfr(guilds, function(g) {
  purrr::map_dfr(metrics, function(m) {
    make_pairs(get_dist(g, m), data, g, m)
  })
}) |>
  dplyr::mutate(
    guild  = factor(guild,  levels = guilds),
    metric = factor(metric, levels = metrics)
  )

# Step 4. Fit pairwise turnover models
fit_turnover <- function(metric_name, guild_name) {
  d <- pairs.all |> dplyr::filter(metric == metric_name, guild == guild_name)
  brms::brm(
    dissimilarity ~ age_diff_s + same_crop + same_soil + same_topo +
      (1 | forest_i) + (1 | forest_j),
    data = d,
    family = brms::Beta(),
    prior = c(
      brms::prior(normal(0, 1), class = "b"),
      brms::prior(normal(0, 1), class = "Intercept"),
      brms::prior(lognormal(0, 1), class = "phi"),
      brms::prior(lognormal(-1, 0.5), class = "sd")
    ),
    chains = 4, cores = 4, iter = 4000, warmup = 2000, seed = 42,
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    file = file.path(here::here("03_models"), paste0("fit_turnover_", model_key(metric_name, guild_name)))
  )
}

turnover_fits <- purrr::map(
  setNames(
    purrr::cross2(metrics, guilds),
    purrr::map_chr(purrr::cross2(metrics, guilds), ~ model_key(.x[[1]], .x[[2]]))
  ),
  function(x) {
    cat(sprintf("  %s — %s\n", x[[1]], x[[2]]))
    fit_turnover(x[[1]], x[[2]])
  }
)

# Step 5. Turnover posterior probabilities
turnover_probs <- purrr::map_dfr(
  purrr::cross2(metrics, guilds),
  function(x) {
    m <- x[[1]]; g <- x[[2]]
    fit <- turnover_fits[[model_key(m, g)]]
    post <- brms::as_draws_df(fit)
    tibble::tibble(
      guild = g,
      metric = m,
      P_increase_with_age = mean(post$b_age_diff_s > 0),
      age_effect_mean = round(mean(post$b_age_diff_s), 3),
      age_effect_lci = round(quantile(post$b_age_diff_s, 0.025), 3),
      age_effect_uci = round(quantile(post$b_age_diff_s, 0.975), 3)
    )
  }
) |>
  dplyr::mutate(
    guild = factor(guild, levels = guilds),
    metric = factor(metric, levels = metrics)
  )

# Step 6. Compute beta dispersion per sample
disp.all <- purrr::map_dfr(guilds, function(g) {
  purrr::map_dfr(metrics, function(m) {
    bd <- vegan::betadisper(get_dist(g, m), data$forest)
    tibble::tibble(
      sample = data$sample,
      forest = data$forest,
      age = data$age,
      categorie = data$categorie,
      crop_type = data$crop_type,
      crop_year = data$crop_year,
      soil_type = data$soil_type,
      topography = data$topography,
      dispersion = bd$distances,
      guild = g,
      metric = m
    )
  })
}) |>
  dplyr::mutate(
    is_ogf = ifelse(categorie == "OGF", 1, 0),
    age_model = ifelse(is_ogf == 1, NA, age),
    age_s = age_model / max(age_model, na.rm = TRUE),
    crop_year_model = ifelse(is_ogf == 1, 0, crop_year),
    crop_type = factor(
      dplyr::case_when(is_ogf == 1 ~ "cacao", TRUE ~ crop_type),
      levels = c("cacao", "cacao_cafe", "cafe", "coton", "igname", "mais", "riz")
    ),
    topography = factor(topography),
    soil_type  = factor(soil_type),
    guild = factor(guild,  levels = guilds),
    metric = factor(metric, levels = metrics),
    dispersion = pmax(dispersion, 1e-6)
  )

# Step 7. Fit dispersion models
fit_dispersion <- function(metric_name, guild_name) {
  d <- disp.all |>
    dplyr::filter(metric == metric_name, guild == guild_name, is_ogf == 0)
  brms::brm(
    dispersion ~ s(age_s) + crop_year_model + crop_type + topography + soil_type +
      (1 | forest),
    data = d,
    family = lognormal(),
    prior = c(
      brms::prior(normal(0, 1), class = "b"),
      brms::prior(lognormal(-1, 0.5), class = "sds"),
      brms::prior(lognormal(-1, 0.5), class = "sd"),
      brms::prior(lognormal(-1, 0.5), class = "sigma")
    ),
    chains = 4, cores = 4, iter = 4000, warmup = 2000, seed = 42,
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    file = file.path(here::here("03_models"), paste0("fit_disp_", model_key(metric_name, guild_name)))
  )
}

disp_fits <- purrr::map(
  setNames(
    purrr::cross2(metrics, guilds),
    purrr::map_chr(purrr::cross2(metrics, guilds), ~ model_key(.x[[1]], .x[[2]]))
  ),
  function(x) {
    cat(sprintf("  %s — %s\n", x[[1]], x[[2]]))
    fit_dispersion(x[[1]], x[[2]])
  }
)

# Step 8. Dispersion posterior probabilities
get_modal <- function(data, metric_name, guild_name, is_ogf_val = 0) {
  d <- data |> dplyr::filter(metric == metric_name, guild == guild_name,
                             is_ogf == is_ogf_val)
  list(
    crop = names(sort(table(d$crop_type), decreasing = TRUE))[1],
    topo = names(sort(table(d$topography), decreasing = TRUE))[1],
    soil = names(sort(table(d$soil_type), decreasing = TRUE))[1],
    cyear = mean(d$crop_year_model, na.rm = TRUE)
  )
}

disp_probs <- purrr::map_dfr(purrr::cross2(metrics, guilds), function(x) {
  m <- x[[1]]
  g <- x[[2]]
  fit <- disp_fits[[model_key(m, g)]]
  ref <- get_modal(disp.all, m, g)
  make_nd <- function(age_val) data.frame(
    age_s = age_val, crop_year_model = ref$cyear,
    crop_type = ref$crop, topography = ref$topo,
    soil_type = ref$soil, forest = NA
  )
  
  p_early <- brms::posterior_epred(fit, newdata = make_nd(0.1), re_formula = NA)
  p_late  <- brms::posterior_epred(fit, newdata = make_nd(0.9), re_formula = NA)
  
  tibble::tibble(
    guild = g, metric = m,
    P_dispersion_declines = round(mean(p_late < p_early), 3),
    disp_change_mean = round(mean(p_late - p_early), 3),
    disp_change_lci = round(quantile(p_late - p_early, 0.025), 3),
    disp_change_uci = round(quantile(p_late - p_early, 0.975), 3)
  )
}) |>
  dplyr::mutate(guild = factor(guild, levels = guilds),
                metric = factor(metric, levels = metrics))

# Step 9. Compute distance of each SF plot to OGF centroid
dist_ogf.all <- purrr::map_dfr(guilds, function(g) {
  purrr::map_dfr(metrics, function(m) {
    dm <- as.matrix(get_dist(g, m))
    ogf_idx <- which(data$categorie == "OGF")
    sf_idx <- which(data$categorie != "OGF")
    tibble::tibble(
      sample = data$sample[sf_idx],
      forest = data$forest[sf_idx],
      age = data$age[sf_idx],
      crop_type = data$crop_type[sf_idx],
      crop_year = data$crop_year[sf_idx],
      soil_type = data$soil_type[sf_idx],
      topography = data$topography[sf_idx],
      dist_to_ogf = rowMeans(dm[sf_idx, ogf_idx]),
      guild = g,
      metric = m
    )
  })
}) |>
  dplyr::mutate(
    age_s = age / max(age, na.rm = TRUE),
    crop_year_model = crop_year,
    crop_type = factor(crop_type, levels = c("cacao", "cacao_cafe", "cafe", "coton", "igname", "mais", "riz")),
    topography = factor(topography),
    soil_type = factor(soil_type),
    guild = factor(guild,  levels = guilds),
    metric = factor(metric, levels = metrics),
    dist_to_ogf = pmin(pmax(dist_to_ogf, 1e-4), 1 - 1e-4)
  )

# Step 10. Fit distance-to-OGF models
fit_dist_ogf <- function(metric_name, guild_name) {
  d <- dist_ogf.all |> dplyr::filter(metric == metric_name, guild == guild_name)
  brms::brm(
    dist_to_ogf ~ s(age_s) + crop_year_model + crop_type + topography + soil_type +
      (1 | forest),
    data  = d,
    family = brms::Beta(),
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(normal(0, 1), class = "Intercept"),
      prior(lognormal(-1, 0.5), class = "sds"),
      prior(lognormal(-1, 0.5), class = "sd"),
      prior(lognormal(0, 1), class = "phi")
    ),
    chains = 4, cores = 4, iter = 4000, warmup = 2000, seed = 42,
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    file = file.path(here::here("03_models"), paste0("fit_distogf_", model_key(metric_name, guild_name)))
  )
}

dist_ogf_fits <- purrr::map(
  setNames(
    purrr::cross2(metrics, guilds),
    purrr::map_chr(purrr::cross2(metrics, guilds), ~ model_key(.x[[1]], .x[[2]]))
  ),
  function(x) {
    cat(sprintf("  %s — %s\n", x[[1]], x[[2]]))
    fit_dist_ogf(x[[1]], x[[2]])
  }
)

# Step 11. Distance-to-OGF posterior probabilities
dist_ogf_probs <- purrr::map_dfr(purrr::cross2(metrics, guilds), function(x) {
  m <- x[[1]]; g <- x[[2]]
  fit <- dist_ogf_fits[[model_key(m, g)]]
  ref <- get_modal(
    dist_ogf.all |> dplyr::mutate(is_ogf = 0, crop_year_model = crop_year),
    m, g
  )
  
  make_nd <- function(age_val) data.frame(
    age_s = age_val, crop_year_model = ref$cyear,
    crop_type = ref$crop, topography = ref$topo,
    soil_type = ref$soil, forest = NA
  )
  
  p_early <- brms::posterior_epred(fit, newdata = make_nd(0.1), re_formula = NA)
  p_late <- brms::posterior_epred(fit, newdata = make_nd(0.9), re_formula = NA)
  
  tibble::tibble(
    guild = g, 
    metric = m,
    P_recovery = round(mean(p_late < p_early), 3),
    dist_change_mean = round(mean(p_late - p_early), 3),
    dist_change_lci = round(quantile(p_late - p_early, 0.025), 3),
    dist_change_uci = round(quantile(p_late - p_early, 0.975), 3),
    dist_early_mean = round(mean(p_early), 3),
    dist_late_mean = round(mean(p_late),  3)
  )
}) |>
  dplyr::mutate(guild = factor(guild, levels = guilds),
                metric = factor(metric, levels = metrics))

write.csv(turnover_probs, here::here("05_results/turnover_posteriors.csv"), row.names = FALSE)
write.csv(disp_probs, here::here("05_results/dispersion_posteriors.csv"), row.names = FALSE)
write.csv(dist_ogf_probs, here::here("05_results/dist_ogf_posteriors.csv"), row.names = FALSE)

get_preds <- function(fits, data_long, response_var, age_var = "age_s") {
  purrr::map_dfr(purrr::cross2(metrics, guilds), function(x) {
    m <- x[[1]]; g <- x[[2]]
    fit <- fits[[model_key(m, g)]]
    d <- data_long |> dplyr::filter(metric == m, guild == g)
    
    cyear_col <- if ("crop_year_model" %in% names(d)) "crop_year_model" else "crop_year"
    ref <- list(
      crop = names(sort(table(d$crop_type), decreasing = TRUE))[1],
      topo = names(sort(table(d$topography), decreasing = TRUE))[1],
      soil = names(sort(table(d$soil_type),  decreasing = TRUE))[1],
      cyear = mean(d[[cyear_col]], na.rm = TRUE)
    )
    
    newdat <- data.frame(
      age_s = seq(0.01, 0.99, length.out = 100),
      crop_year_model = ref$cyear,
      crop_type = ref$crop,
      topography = ref$topo,
      soil_type = ref$soil,
      forest = NA
    )
    
    tidybayes::add_epred_draws(fit, newdata = newdat, re_formula = NA) |>
      dplyr::mutate(guild = g, metric = m)
  }) |>
    dplyr::mutate(guild  = factor(guild,  levels = guilds),
                  metric = factor(metric, levels = metrics))
}

pred_dist <- get_preds(
  dist_ogf_fits, 
  dist_ogf.all |>
    dplyr::mutate(is_ogf = 0, crop_year_model = crop_year), "dist_to_ogf"
  )

p_dist_ogf <- ggplot2::ggplot(
  pred_dist,
  ggplot2::aes(x = age_s, y = .epred, colour = guild, fill = guild)
) +
  tidybayes::stat_lineribbon(.width = 0.95, alpha = 0.20, linewidth = 0.9) +
  ggplot2::geom_point(
    data    = dist_ogf.all,
    mapping = ggplot2::aes(x = age_s, y = dist_to_ogf, colour = guild),
    shape = 16, size = 1.8, alpha = 0.7, inherit.aes = FALSE
  ) +
  ggplot2::scale_colour_manual(values = guild_colours, labels = c("Bacteria", "AM fungi")) +
  ggplot2::scale_fill_manual(values = guild_fills, labels = c("Bacteria", "AM fungi")) +
  ggplot2::facet_wrap(guild ~ metric, scales = "free_y", ncol = 4) +
  ggplot2::scale_x_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = function(x) round(x * max(data$age, na.rm = TRUE))
  ) +
  ggplot2::labs(
    x = "Time since abandonment (years)", y = "Compositional distance to OGF",
    colour = "Guild", fill = "Guild",
    caption = "Shaded band: 95% CI | Predictions marginalised over confounders"
  ) +
  theme_Publication() +
  ggplot2::theme(legend.position = "bottom", strip.text = ggplot2::element_text(size = 9),
                 panel.spacing = ggplot2::unit(0.6, "lines"))

export_plot(p_dist_ogf, here::here("04_figures/Fig_dist_to_OGF.png"), width = 13, height = 8)

pred_disp <- get_preds(disp_fits, disp.all |> dplyr::filter(is_ogf == 0), "dispersion")

p_disp <- ggplot2::ggplot(
  pred_disp,
  ggplot2::aes(x = age_s, y = .epred, colour = guild, fill = guild)
) +
  tidybayes::stat_lineribbon(.width = 0.95, alpha = 0.20, linewidth = 0.9) +
  ggplot2::geom_point(
    data = disp.all |> dplyr::filter(is_ogf == 0),
    mapping = ggplot2::aes(x = age_s, y = dispersion, colour = guild),
    shape = 16, size = 1.8, alpha = 0.7, inherit.aes = FALSE
  ) +
  ggplot2::scale_colour_manual(values = guild_colours, labels = c("Bacteria", "AM fungi")) +
  ggplot2::scale_fill_manual(values = guild_fills,   labels = c("Bacteria", "AM fungi")) +
  ggplot2::facet_wrap(guild ~ metric, scales = "free_y", ncol = 4) +
  ggplot2::scale_x_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = function(x) round(x * max(data$age, na.rm = TRUE))
  ) +
  ggplot2::labs(
    x = "Time since abandonment (years)", y = "Beta-dispersion (distance to centroid)",
    colour = "Guild", fill = "Guild",
    caption = "Shaded band: 95% CI | Predictions marginalised over confounders"
  ) +
  theme_Publication() +
  ggplot2::theme(legend.position = "bottom", strip.text = ggplot2::element_text(size = 9),
                 panel.spacing = ggplot2::unit(0.6, "lines"))

export_plot(p_disp, here::here("04_figures/Fig_dispersion.png"), width = 13, height = 8)


# Step 12. Within-OGF vs within-SF dispersion comparison
compare_ogf_sf_dispersion <- function(dist_mat, data, guild_name, metric_name) {
  # Dispersion within OGF plots
  ogf_idx <- which(data$categorie == "OGF")
  sf_idx  <- which(data$categorie != "OGF")
  
  if (length(ogf_idx) < 2) return(NULL)
  
  dist_mat <- as.matrix(dist_mat)
  
  # Mean pairwise distance within OGF
  ogf_pairs <- dist_mat[ogf_idx, ogf_idx]
  mean_ogf_dist <- mean(ogf_pairs[lower.tri(ogf_pairs)])
  
  # Mean pairwise distance within SF
  sf_pairs <- dist_mat[sf_idx, sf_idx]
  mean_sf_dist <- mean(sf_pairs[lower.tri(sf_pairs)])
  
  # Betadisper within OGF (distance to OGF centroid)
  if (length(unique(data$forest[ogf_idx])) > 1) {
    bd_ogf <- vegan::betadisper(
      as.dist(dist_mat[ogf_idx, ogf_idx]),
      data$forest[ogf_idx]
    )
    mean_ogf_disp <- mean(bd_ogf$distances)
  } else {
    mean_ogf_disp <- NA
  }
  
  # Betadisper within SF
  bd_sf <- vegan::betadisper(
    as.dist(dist_mat[sf_idx, sf_idx]),
    data$forest[sf_idx]
  )
  mean_sf_disp <- mean(bd_sf$distances)
  
  tibble::tibble(
    guild = guild_name,
    metric = metric_name,
    mean_ogf_pw_dist = round(mean_ogf_dist, 3),
    mean_sf_pw_dist = round(mean_sf_dist, 3),
    mean_ogf_disp = round(mean_ogf_disp, 3),
    mean_sf_disp = round(mean_sf_disp, 3),
    n_ogf = length(ogf_idx),
    n_sf = length(sf_idx)
  )
}

ogf_sf_disp <- purrr::map_dfr(guilds, function(g) {
  purrr::map_dfr(metrics, function(m) {
    compare_ogf_sf_dispersion(get_dist(g, m), data, g, m)
  })
}) |>
  dplyr::mutate(
    guild = factor(guild, levels = guilds),
    metric = factor(metric, levels = metrics)
  )

write.csv(ogf_sf_disp, here::here("05_results/ogf_sf_dispersion_comparison.csv"), row.names = FALSE)

# Step 13. Sensitivity: turnover model restricted to between-site pairs only
fit_turnover_between_sites <- function(metric_name, guild_name) {
  d <- pairs.all |>
    dplyr::filter(
      metric == metric_name,
      guild == guild_name,
      forest_i != forest_j      # between-site pairs only
    )
  
  if (nrow(d) < 10) {
    warning(sprintf("Too few between-site pairs for %s %s", metric_name, guild_name))
    return(NULL)
  }
  
  brms::brm(
    dissimilarity ~ age_diff_s + same_crop + same_soil + same_topo +
      (1 | forest_i) + (1 | forest_j),
    data = d,
    family = Beta(),
    prior = c(
      brms::prior(normal(0, 1), class = "b"),
      brms::prior(normal(0, 1), class = "Intercept"),
      brms::prior(lognormal(0, 1), class = "phi"),
      brms::prior(lognormal(-1, 0.5), class = "sd")
    ),
    chains = 4, cores = 4, iter = 4000, warmup = 2000, seed = 42,
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    file = file.path(here::here("03_models"), paste0("fit_turnover_btwsites_", model_key(metric_name, guild_name)))
  )
}

turnover_btwsite_fits <- purrr::map(
  setNames(
    purrr::cross2(metrics, guilds),
    purrr::map_chr(purrr::cross2(metrics, guilds), ~ model_key(.x[[1]], .x[[2]]))
  ),
  function(x) fit_turnover_between_sites(x[[1]], x[[2]])
)

# Extract sensitivity probabilities
turnover_sensitivity <- purrr::map_dfr(
  purrr::cross2(metrics, guilds),
  function(x) {
    m <- x[[1]]; g <- x[[2]]
    fit <- turnover_btwsite_fits[[model_key(m, g)]]
    if (is.null(fit)) return(NULL)
    post <- brms::as_draws_df(fit)
    tibble::tibble(
      guild = g, metric = m,
      P_increase_between_sites = round(mean(post$b_age_diff_s > 0), 3),
      age_effect_mean = round(mean(post$b_age_diff_s), 3)
    )
  }
)

write.csv(turnover_sensitivity, here::here("05_results/turnover_sensitivity_betweensites.csv"), row.names = FALSE)