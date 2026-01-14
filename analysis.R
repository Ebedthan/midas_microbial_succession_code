library(ggplot2)
library(brms)
library(patchwork)

theme_Publication <- function(base_size = 14,
                              base_family = "helvetica") {
  
  library(ggthemes)
  library(ggplot2)
  library(grid)
  
  ggthemes::theme_foundation(
    base_size = base_size,
    base_family = base_family
  ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = rel(1.2),
        hjust = 0.5
      ),
      text = element_text(),
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),
      legend.background = element_rect(
        fill = "white",
        colour = NA
      ),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size = unit(0.2, "cm"),
      legend.margin = margin(t = 0, unit = "cm"),
      legend.title = element_text(face = "italic"),
      strip.background = element_rect(
        fill = "white",
        colour = NA
      ),
      strip.text = element_text(face = "bold"),
      plot.margin = unit(c(10, 5, 5, 5), "mm")
    )
}

export_plot <- function(myplot, file, width = 8, height = 6, units = "in", type = "tiff", bg = "white") {
  require(Cairo)
  require(cli)
  
  # Export plot using Cairo
  Cairo::Cairo(width = width,
               height = height,
               file = file,
               type = type,
               bg = bg,
               units = units,
               dpi = 300)
  plot(myplot)
  dev.off()
  
  cli::cli_alert_success(paste("Exporting plot to", file, "succeeded"))
  
  # Show plot
  myplot
}

# 0. ===========================================================================
## loading forests and metadata
data <- read.csv("data.csv")
asv.bact <- read.csv("asv_table.csv", row.names = 1, check.names = F)
asv.amf <- read.csv("amf_table.csv", row.names = 1, check.names = F)
asv.bact <- t(asv.bact)
asv.amf <- t(asv.amf)

## computing relative abundance
asv.bact.rel <- asv.bact / rowSums(asv.bact)
asv.amf.rel <- asv.amf / rowSums(asv.amf)

## computing Bray-Curtis distance 
bray.mat.bact <- as.matrix(vegan::vegdist(asv.bact.rel, method = "bray"))
bray.mat.amf <- as.matrix(vegan::vegdist(asv.amf.rel, method = "bray"))

# =========================================================
# Research question 1
# Do microbial alpha-diversity increase with successional
# age, and do bacteria and AM fungi differ in their 
# responses ?
# =========================================================

# Plot diversity vs. age group
data.div <- data |>
  dplyr::select(
    categorie, 
    dplyr::starts_with("amf"), 
    dplyr::starts_with("bacterial")
  ) |>
  tidyr::pivot_longer(
    cols = amf_richness:bacterial_simpson,
    names_to = "estimate",
    values_to = "value"
  )

facet.labels <- c(
  amf_richness        = "AMF richness",
  amf_shannon         = "AMF Shannon diversity",
  amf_simpson         = "AMF Simpson diversity",
  bacterial_richness  = "Bacterial richness",
  bacterial_shannon   = "Bacterial Shannon diversity",
  bacterial_simpson   = "Bacterial Simpson diversity"
)

p1 <- ggplot(data.div, aes(x = categorie, y = value)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.6) +
  facet_wrap(~estimate, scales = "free", labeller = labeller(estimate = facet.labels)) +
  labs(x = "Age category", y = NULL) +
  theme_Publication() 

export_plot(p1, "Fig2.png", width = 13, height = 8, type = "png")

# 
data.long <- data |>
  dplyr::select(
    site = forest,
    sample, categorie, age,
    bacterial_shannon, amf_shannon
  ) |>
  tidyr::pivot_longer(
    cols = c(bacterial_shannon, amf_shannon),
    names_to = "guild",
    values_to = "diversity"
  ) |>
  dplyr::mutate(
    guild = factor(
      guild,
      levels = c("bacterial_shannon", "amf_shannon"),
      labels = c("Bacteria", "AMF")
    )
  )

bayes.div <- brm(
  diversity ~ s(age, by = guild) + guild + (1 | site),
  data = data.long,
  family = lognormal(),
  prior = c(
    # regression coeeficients
    prior(normal(0, 1), class = "b"),
    # smooth term standard deviation
    prior(lognormal(-1, 0.5), class = "sds"),
    # random effect SDs
    prior(lognormal(-1, 0.5), class = "sd"),
    # residual SD
    prior(lognormal(-1, 0.5), class = "sigma")
  ),
  chains = 4,
  cores = 8,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

post.div <- as_draws_df(bayes.div)

# probability of increasing diversity with age
P.div.increase.bact <- mean(post.div$`bs_sage:guildBacteria_1` > 0)
P.div.increase.amf <- mean(post.div$`bs_sage:guildAMF_1` > 0)
P.div.increase.bact
P.div.increase.amf

# probability that AM fungi respond more strongly than
# bacteria
P.amf.more.successional <- mean(
  post.div$`bs_sage:guildAMF_1` > post.div$`bs_sage:guildBacteria_1`
)
P.amf.more.successional

# plot
newdat.div <- expand.grid(
  age = seq(min(data.long$age, na.rm = T),
            max(data.long$age, na.rm = T),
            length.out = 100),
  guild = levels(data.long$guild),
  site = NA
)

pred.div <- tidybayes::add_epred_draws(
  bayes.div,
  newdata = newdat.div,
  re_formula = NA
)

ggplot(pred.div, aes(x = age, y = .epred, color = guild)) +
  tidybayes::stat_lineribbon(
    aes(fill = guild),
    .width = 0.95,
    alpha = 0.3
  ) +
  labs(
    x = "Time since abandonment (years)",
    y = "Alpha diversity (Hill-based)",
    color = "Guild",
    fill = "Guild"
  ) +
  theme_Publication()

# ======================================================================
# Research question 2
# Do microbial communities reassemble toward old-growth forest 
# reference states during succession?
# ======================================================================
## presence-absence (q = 0)
asv.bact.pa <- asv.bact > 0
asv.amf.pa <- asv.amf > 0

## disimilarities
### Sorensen
dist.bact.q0 <- vegan::vegdist(asv.bact.pa, method = "bray")
dist.amf.q0 <- vegan::vegdist(asv.amf.pa, method = "bray")

### Horn
dist.bact.q1 <- vegan::vegdist(asv.bact.rel, method = "horn")
dist.amf.q1  <- vegan::vegdist(asv.amf.rel,  method = "horn")

### Morisita–Horn
dist.bact.q2 <- vegan::vegdist(asv.bact, method = "morisita")
dist.amf.q2  <- vegan::vegdist(asv.amf,  method = "morisita")

### Bray–Curtis
dist.bact.bray <- vegan::vegdist(asv.bact.rel, method = "bray")
dist.amf.bray  <- vegan::vegdist(asv.amf.rel,  method = "bray")

## NMDS
nmds.bact.q1 <- vegan::metaMDS(dist.bact.q1, k = 2, trymax = 100)
nmds.amf.q1  <- vegan::metaMDS(dist.amf.q1,  k = 2, trymax = 100)
datage <- data$age
datage[is.na(datage)] <- 100
p2 <- ggplot(
  nmds.bact.q1$points, 
  aes(
    x = MDS1, 
    y = MDS2, 
    colour = data$forest,
    size = datage
    )
  ) +
  geom_point() +
  ggsci::scale_color_aaas() +
  labs(x = "NMDS1", y = "NMDS2", colour = "Forests", size = "Successional age") +
  theme_Publication()

p3 <- ggplot(nmds.amf.q1$points, aes(x = MDS1, y = MDS2, colour = data$forest,
                                     size = datage)) +
  geom_point() +
  ggsci::scale_color_aaas() +
  labs(x = "NMDS1", y = "NMDS2", colour = "Forests", size = "Successional age") +
  theme_Publication()

p4 <- p2 + p3 + 
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

export_plot(p4, "Fig3.png", type = "png", width = 13, height = 8)

nmds.bact.q1$stress
nmds.amf.q1$stress
## tests of reassembly
### PERMANOVA
vegan::adonis2(
  dist.bact.q1 ~ categorie,
  data = data,
  permutations = 9999
)

vegan::adonis2(
  dist.amf.q1 ~ categorie,
  data = data,
  permutations = 9999
)

### PERMDIST
bd.bact.q1 <- vegan::betadisper(dist.bact.q1, data$categorie)
bd.amf.q1  <- vegan::betadisper(dist.amf.q1,  data$categorie)

anova(bd.bact.q1)
anova(bd.amf.q1)

vegan::permutest(bd.bact.q1, pairwise = TRUE, permutations = 9999)
vegan::permutest(bd.amf.q1,  pairwise = TRUE, permutations = 9999)

# ====================================================================
# Research question 3
# How rapidly do microbial communities reorganize 
# and stabilize following abandonment ?
# ====================================================================
compute.turnover <- function(data, bray.mat) {
  
  data.ord <- data |>
    dplyr::filter(!is.na(age)) |>
    dplyr::arrange(forest, age)
  
  turnover <- data.ord |>
    dplyr::group_by(forest) |>
    dplyr::mutate(
      next_sample = dplyr::lead(sample),
      next_age    = dplyr::lead(age)
    ) |>
    dplyr::filter(!is.na(next_sample)) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      turnover = bray.mat[sample, next_sample],
      age_mid  = (age + next_age) / 2
    ) |>
    dplyr::ungroup()
  
  turnover
}

turnover.bact <- compute.turnover(data, as.matrix(dist.bact.bray))
turnover.amf <- compute.turnover(data, as.matrix(dist.amf.bray))

## Turnover rates
turnover.long <- dplyr::bind_rows(
  turnover.bact |>
    dplyr::mutate(guild = "Bacteria"),
  turnover.amf |>
    dplyr::mutate(guild = "AMF")
) |>
  dplyr::mutate(
    guild = factor(guild),
    age_s = age_mid / max(age_mid, na.rm = TRUE)
  )

bayes.turnover <- brm(
  turnover ~ s(age_s, by = guild) + guild + (1 | forest),
  data = turnover.long,
  family = Beta(),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(lognormal(-1, 0.5), class = "sds"),
    prior(lognormal(0, 0.5), class = "phi")
  ),
  chains = 4,
  cores = 8,
  iter = 4000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)
summary(bayes.turnover)

post.bayes.turnover <- as_draws_df(bayes.turnover)

### probability of declining turnover with age
P.amf.turnover.decline <- mean(
  post.bayes.turnover$`bs_sage_s:guildAMF_1` < 0
)
P.amf.turnover.decline

P.bact.turnover.decline <- mean(
  post.bayes.turnover$`bs_sage_s:guildBacteria_1` < 0
  )
P.bact.turnover.decline

### cross-guild comparison
P.amf.slower.turnover <- mean(
  post.bayes.turnover$`bs_sage_s:guildAMF_1` >
    post.bayes.turnover$`bs_sage_s:guildBacteria_1`
)
P.amf.slower.turnover

### stability (dispersion)
disp.long <- data
disp.long$dispersion_bact <- vegan::betadisper(
  vegan::vegdist(asv.bact.rel, method = "bray"),
  data$forest
)$distances

disp.long$dispersion_amf <- vegan::betadisper(
  vegan::vegdist(asv.amf.rel, method = "bray"),
  data$forest
)$distances

disp.long <- disp.long |>
  dplyr::select(site = forest, age, categorie, dispersion_bact, dispersion_amf) |>
  tidyr::pivot_longer(
    cols = c(dispersion_bact, dispersion_amf),
    names_to = "guild",
    values_to = "dispersion"
  ) |>
  dplyr::mutate(
    guild = factor(
      guild,
      levels = c("dispersion_bact", "dispersion_amf"),
      labels = c("Bacteria", "AMF")
    ),
    age_s = age / max(age, na.rm = TRUE)
  )

bayes.disp <- brm(
  dispersion ~ age_s * guild + (1 | site),
  data = disp.long,
  family = lognormal(),
  prior = c(
    # fixed effects 
    prior(normal(0, 0.5), class = "b"),
    # random-effect SD (site-to-site variability)
    prior(lognormal(-1, 0.5), class = "sd"),
    # residual variation (within-site variability)
    prior(lognormal(-1, 0.5), class = "sigma")
  ),
  chains = 4,
  cores = 8,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)
summary(bayes.disp)

post.bayes.disp <- as_draws_df(bayes.disp)

### bacterial convergence
P.bact.disp.decline <- mean(post.bayes.disp$b_age_s < 0)
P.bact.disp.decline

### AMF convergence
P.amf.disp.decline <- mean(
  (post.bayes.disp$b_age_s + post.bayes.disp$`b_age_s:guildAMF`) < 0
)
P.amf.disp.decline

### cross-guild comparison
P.amf.more.convergence <- mean(
  post.bayes.disp$`b_age_s:guildAMF` < 0
)
P.amf.more.convergence

### Plot Beta dispersion vs time
newdisp <- disp.long |>
  dplyr::distinct(guild) |>
  tidyr::expand_grid(
    age_s = seq(
      min(disp.long$age_s, na.rm = TRUE),
      max(disp.long$age_s, na.rm = TRUE),
      length.out = 200
    ),
    site = NA
  )

pred.disp <- bayes.disp |>
  tidybayes::add_epred_draws(
    newdata = newdisp,
    re_formula = NA
  )

annot.disp <- data.frame(
  guild = c("AMF", "Bacteria"),
  label = c(
    "P(convergence) = 0.48",
    "P(convergence) = 0.37"
  ),
  x = 0.8,
  y = c(0.75, 0.7)
)

p7 <- ggplot() +
  geom_point(
    data = disp.long,
    aes(x = age_s, y = dispersion, colour = guild),
    alpha = 0.5, 
    size = 2 
  ) + geom_text(
    data = annot.disp,
    aes(x = x, y = y, label = label, colour = guild),
    hjust = 0,
    size = 4
  ) +
  tidybayes::stat_lineribbon(
    data = pred.disp,
    aes(x = age_s, y = .epred, fill = guild),
    .width = 0.95,
    alpha = 0.25
  ) +
  stat_summary(
    data = pred.disp,
    aes(x = age_s, y = .epred, colour = guild),
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  scale_colour_manual(
    values = c("Bacteria" = "#0072B2", "AMF" = "#D55E00")
  ) +
  scale_fill_manual(
    values = c("Bacteria" = "#0072B2", "AMF" = "#D55E00")
  ) +
  labs(
    x = "Normalized time since abandonment",
    y = "Beta-dispersion (distance to centroid)",
    colour = "Guild",
    fill = "Guild"
  ) +
  theme_Publication()

export_plot(p7, "Fig4.png", type = "png", width = 13, height = 8)

# ========================================================================
# Research question 4 
# Is microbial succession deterministic and stage-structured, 
# or diffuse and probabilistic ?
# ========================================================================

## Computing distance to reference
### computing distance to reference as centroid
# extracting OGF rows ids
forest.idx <- data$categorie == "OGF"
data.succ <- data
data.succ$dref_bact <- apply(
  bray.mat.bact[, forest.idx, drop = FALSE],
  1,
  mean
)
data.succ$dref_amf <- apply(
  bray.mat.amf[, forest.idx, drop = FALSE],
  1,
  mean
)

data.succ <- data.succ |>
  dplyr::filter(!is.na(age)) |>
  dplyr::select(sample, forest, categorie, age, dref_amf, dref_bact) |>
  tidyr::pivot_longer(
    cols = c(dref_bact, dref_amf),
    names_to = "guild",
    values_to = "dref"
  ) |>
  dplyr::mutate(
    guild = factor(
      guild,
      levels = c("dref_bact", "dref_amf"),
      labels = c("Bacteria", "AMF")
    )
  )

bayes.succ <- brm(
  dref ~ s(age, by = guild) + guild + (1 | forest),
  data = data.succ,
  family = Beta(),
  prior = c(
    # fixed effects
    prior(normal(0, 1), class = "b"),
    # smooth term
    prior(student_t(3, 0, 0.5), class = "sds"),
    # random effect sd
    prior(exponential(2), class = "sd"),
    # precision param
    prior(exponential(1), class = "phi")
  ),
  chains = 4,
  cores = 8,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)
summary(bayes.succ)

post.bayes.succ <- as_draws_df(bayes.succ)

# probability that AMF distance to forest decrease with age
P.amf.decline <- mean(post.bayes.succ$`bs_sage:guildAMF_1` < 0)
P.amf.decline

# probability that bacteria distance to forest decrease with age
P.bact.decline <- mean(post.bayes.succ$`bs_sage:guildBacteria_1` < 0)
P.bact.decline

# probability AMF age effect more negative than bacterial
P.amf.more.recovery <- mean(
  post.bayes.succ$`bs_sage:guildAMF_1` <
    post.bayes.succ$`bs_sage:guildBacteria_1`
)
P.amf.more.recovery

# Fig 5. distance to reference vs time plot
data.succ <- data.succ |>
  dplyr::mutate(similarity = 1 - dref)

newdat.succ <- data.succ |>
  dplyr::distinct(guild) |>
  tidyr::expand_grid(
    age = seq(
      min(data.succ$age, na.rm = TRUE),
      max(data.succ$age, na.rm = TRUE),
      length.out = 200
    ),
    forest = NA
  )

pred.dref <- bayes.succ |>
  tidybayes::add_epred_draws(
    newdata = newdat.succ,
    re_formula = NA
  ) |>
  dplyr::mutate(similarity = 1 - .epred)

annot <- data.frame(
  guild = c("AMF", "Bacteria"),
  label = c(
    "P(decline) = 0.80",
    "P(decline) = 0.67"
  ),
  x = max(data.succ$age) * 0.7,
  y = c(0.03, 0.05)
)

p5 <- ggplot() +
  geom_point(data = data.succ, aes(x = age, y = similarity, colour = guild), alpha = 0.5, size = 2) +
  geom_text(
    data = annot,
    aes(x = x, y = y, label = label, colour = guild),
    hjust = 0,
    size = 4
  ) +
  tidybayes::stat_lineribbon(data = pred.dref, aes(x = age, y = similarity, fill = guild), .width = 0.95, alpha = 0.25) +
  stat_summary(data = pred.dref, aes(x = age, y = similarity, colour = guild), fun = mean, geom = "line", linewidth = 1.2) +
  scale_colour_manual(
    values = c("Bacteria" = "#0072B2", "AMF" = "#D55E00"),
  ) +
  scale_fill_manual(
    values = c("Bacteria" = "#0072B2", "AMF" = "#D55E00"),
  ) +
  labs(
    x = "Time since abandonment (years)",
    y = "Similarity to reference forest (1 - Bray-Curtis)",
    colour = "Guild",
    fill = "Guild"
  ) +
  theme_Publication()

export_plot(p5, "Fig5.png", type = "png", width = 13, height = 8)

# 4. indicator taxa and successional microbial states ====================
group.stage <- factor(
  data$categorie,
  levels = c("[1-10]", "[11-20]", "[21-30]", "[30 <[", "OGF")
)

names(group.stage) <- data$sample

asv.bact.hel <- vegan::decostand(asv.bact, method = "hellinger")
asv.amf.hel <- vegan::decostand(asv.amf, method = "hellinger")

indval.bact <- indicspecies::multipatt(
  asv.bact.hel,
  group.stage,
  func = "IndVal.g",
  control = permute::how(nperm = 999)
)
indval.bact$sign$p.adj <- p.adjust(indval.bact$sign$p.value, method = "BH")

summary(indval.bact, indvalcomp = TRUE)
saveRDS(indval.bact, "indval.bact.rds")

indval.amf <- indicspecies::multipatt(
  asv.amf.hel,
  group.stage,
  func = "IndVal.g",
  control = permute::how(nperm = 999)
)
indval.amf$sign$p.adj <- p.adjust(indval.amf$sign$p.value, method = "BH")

summary(indval.amf, indvalcomp = TRUE)
saveRDS(indval.amf, "indval.amf.rds")

get.indicators <- function(indval.obj, alpha = 0.05) {
  
  res <- indval.obj$sign
  res$taxon <- rownames(res)
  res <- res |> 
    tidyr::pivot_longer(dplyr::starts_with("s."), names_to = "s", values_to = "estimate")
  
  data.frame(
    taxon = res$taxon,
    group = res$s,
    stat  = res$stat,
    p.adj  = res$p.adj
  ) |>
    dplyr::filter(p.adj <= alpha) |>
    dplyr::arrange(group, dplyr::desc(stat))
}

indicators.bact <- get.indicators(indval.bact, 0.05)
indicators.amf  <- get.indicators(indval.amf, 0.05)

forest.specialists.bact <- indicators.bact |>
  dplyr::filter(group == "Old forest")

forest.specialists.amf <- indicators.amf |>
  dplyr::filter(group == "Old forest")


# =======================================================================
# PLOTS
# =======================================================================
facet.labels.2 <- c(
  dref_amf    = "AM fungi",
  dref_bact   = "Bacteria"
)

p6 <- ggplot(data.traj, aes(age, dref, group = forest)) +
  geom_path(aes(colour = forest), alpha = 0.6) +
  geom_point(aes(colour = forest), size = 2) +
  facet_wrap(~ guild, ncol = 1, labeller = labeller(guild = facet.labels.2)) +
  labs(
    x = "Time since abandonment (years)",
    y = "Distance to reference forest",
    colour = "Forest"
  ) +
  ggsci::scale_color_lancet() +
  theme_Publication()

export_plot(p6, "Fig S1.png", type = "png", width = 13, height = 8)

p8 <- ggplot(disp.long, aes(age_s, dispersion, group = site)) +
  geom_path(aes(colour = site), alpha = 0.6) +
  geom_point(aes(colour = site), size = 2) +
  facet_wrap(~ guild, ncol = 1) +
  labs(
    x = "Normalized time since abandonment",
    y = "Beta-dispersion",
    colour = "Site"
  ) +
  ggsci::scale_color_lancet() +
  theme_Publication()

export_plot(p8, "FigS2.png", type = "png", width = 13, height = 8)

# comparative recovery rates plot
age.grid <- seq(
  min(data.traj$age, na.rm = TRUE),
  max(data.traj$age, na.rm = TRUE),
  length.out = 200
)

newdata <- expand.grid(
  age   = age.grid,
  guild = levels(data.traj$guild),
  forest = NA   # marginal over random effects
)

pred <- bayes.gam.guild |>
  tidybayes::add_epred_draws(newdata = newdata, re_formula = NA)

pred.diff <- pred |>
  dplyr::ungroup() |>
  dplyr::select(.draw, age, guild, .epred) |>
  tidyr::pivot_wider(
    names_from = guild,
    values_from = .epred
  ) |>
  dplyr::mutate(diff = dref_amf - dref_bact)

pred.diff.summary <- pred.diff |>
  dplyr::group_by(age) |>
  tidybayes::mean_qi(diff, .width = c(0.5, 0.8, 0.95))

p9 <- ggplot(pred.diff.summary, aes(x = age, y = diff)) +
  geom_ribbon(
    aes(ymin = .lower, ymax = .upper),
    fill = "grey70",
    alpha = 0.4
  ) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Time since abandonment (years)",
    y = expression(Delta~"Distance to forest (AMF – Bacteria)")
  ) +
  theme_Publication()

export_plot(p9, "Fig6.png", type = "png", width = 13, height = 8)


# Shannon diversity for plateau reaching
alpha.long <- data |>
  dplyr::filter(!is.na(age)) |>
  dplyr::select(
    age, bacterial_shannon,
    amf_shannon, successional_stage
  ) |>
  tidyr::pivot_longer(
    cols = c(bacterial_shannon, amf_shannon),
    names_to = "guild",
    values_to = "diversity"
  ) |>
  dplyr::mutate(
    guild = factor(
      guild,
      levels = c("bacterial_shannon", "amf_shannon"),
      labels = c("Bacteria", "AMF")
    )
  )
alpha.long.early <- alpha.long |>
  dplyr::filter(successional_stage %in% c("Early", "Mid"))

ggplot(alpha.long.early, aes(x = age, y = diversity)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 4),
    se = TRUE,
    linewidth = 1.2
  ) +
  facet_wrap(~ guild, scales = "free_y") +
  labs(
    x = "Time since abandonment (years)",
    y = "Alpha diversity (Hill-based)"
  ) +
  theme_Publication()

