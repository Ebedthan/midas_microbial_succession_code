# Microbial Secondary Succession in West African Post-Forest Landscapes

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R version](https://img.shields.io/badge/R-%3E%3D4.6-blue.svg)](https://www.r-project.org/)

This repository contains all data and R code associated with the manuscript:

> **Ebou A.E.T., Amani B.H.K., Toure G.P.T., Ehouman E., Zaoui S.V.,
> Toure A.D., Ndiaye S.M., Yapo S.C., Koffi A.B., Fossou R.K.,
> Aussenac R., Zézé A., Koua D.K., and Hérault B.**
> *Microbial secondary succession in West African post-forest landscapes.*

## Study overview

This study investigates the successional dynamics of soil bacterial and
arbuscular mycorrhizal fungal (AMF) communities along post-agricultural
chronosequences in Côte d'Ivoire, spanning six classified forests across
three phytogeographic zones (dry, semi-deciduous, and evergreen). Thirty
plots were sampled across early (1–10 yr), intermediate (11–20, 21–30 yr),
and late (>30 yr) secondary forests, with old-growth forest reference plots
at each site.

Community composition was characterised from 16S rRNA (bacteria) and 18S
rRNA (AM fungi) amplicon sequencing data. Alpha diversity, beta diversity,
compositional turnover, distance to old-growth reference states, and
indicator taxa were evaluated within a Bayesian hierarchical framework
using generalised additive models (GAMs) and mixed-effects beta regression
models.


## Repository structure

```
midas_microbial_succession_code/
│
├── 01_data/                  # Input data files
│   ├── data.csv              # Plot-level metadata (forest, age, crop type,
│   │                         #   soil type, topography, land-use history)
│   ├── asv_table.csv         # Bacterial ASV count table (ASVs × samples)
│   └── amf_table.csv         # AM fungal Virtual Taxa count table (VT × samples)
│
├── 02_R/                     # R analysis scripts
│   ├── richness_fit.R        # Bayesian hierarchical GAMs for alpha diversity
│   │                         #   (richness, Shannon, Simpson; negative binomial
│   │                         #   and lognormal families; guild × age smooths)
│   └── fit_beta.R            # Beta diversity analyses:
│                             #   - NMDS ordination (Horn dissimilarity)
│                             #   - PERMANOVA
│                             #   - Variance partitioning (vegan::varpart)
│                             #   - Pairwise turnover (Beta regression)
│                             #   - Beta-dispersion GAMs (betadisper + brms)
│                             #   - Distance-to-OGF GAMs (Beta family)
│                             #   - Sensitivity analyses (between-site pairs)
│                             #   - Within-OGF vs within-SF dispersion comparison
│
└── .gitignore
```

## Dependencies

All analyses were conducted in **R (≥ 4.3)**. The following packages are
required:

| Package | Version | Purpose |
|---|---|---|
| `brms` | ≥ 2.20 | Bayesian model fitting via Stan |
| `vegan` | ≥ 2.7 | Dissimilarity matrices, NMDS, PERMANOVA, betadisper, varpart |
| `tidybayes` | ≥ 3.0 | Posterior extraction and summarisation |
| `ggplot2` | ≥ 3.4 | Figures |
| `patchwork` | ≥ 1.2 | Figure composition |
| `purrr` | ≥ 1.0 | Functional iteration over model combinations |
| `dplyr` | ≥ 1.1 | Data manipulation |
| `tidyr` | ≥ 1.3 | Data reshaping |
| `entropart` | ≥ 1.6 | Hill-number diversity metrics |
| `indicspecies` | ≥ 1.8 | Indicator species analysis (IndVal) |
| `ggsci` | — | Colour palettes for figures |

Install all required packages with:

```r
install.packages(c(
  "brms", "vegan", "tidybayes", "ggplot2", "patchwork",
  "purrr", "dplyr", "tidyr", "entropart", "indicspecies", "ggsci"
))
```

Stan (the probabilistic programming language used by `brms`) must also be
installed. Follow the instructions at
[mc-stan.org](https://mc-stan.org/users/interfaces/rstan) for your
operating system.

## Reproducibility

All Bayesian models are cached to disk via the `file` argument of `brm()`.
Once fitted, models are reloaded automatically on subsequent runs without
refitting. Model convergence was verified for all models using:

- Potential scale reduction factor: $\hat{R} < 1.01$ for all parameters
- Bulk ESS > 400 and tail ESS > 400 for all parameters
- Visual inspection of trace plots

A `models/` directory will be created automatically in the working
directory to store fitted model objects (`.rds` files). These are excluded
from version control via `.gitignore` due to file size.

## Data

### `data.csv`

Plot-level metadata. Key columns:

| Column | Description |
|---|---|
| `sample` | Unique plot identifier |
| `forest` | Forest site (6 levels: Badénou, Foumbou, Haut-Sassandra, Téné, Niegré, Irobo) |
| `age` | Time since agricultural abandonment (years); `NA` for old-growth forest plots |
| `categorie` | Successional age category (`[1-10]`, `[11-20]`, `[21-30]`, `[30-<[`, `OGF`) |
| `crop_type` | Former crop type (e.g., `cacao`, `cafe`, `mais`) |
| `crop_year` | Years since last cultivation |
| `soil_type` | Dominant soil type |
| `topography` | Topographic position |
| `zone` | Phytogeographic zone (`dry`, `semi-deciduous`, `evergreen`) |

### `asv_table.csv`

Bacterial amplicon sequence variant (ASV) count table. Rows are ASVs,
columns are samples. Produced by the DADA2 pipeline from 16S rRNA V3–V4
amplicons (primers 344F / 802R), sequenced on Illumina MiSeq
(2 × 300 bp).

### `amf_table.csv`

AM fungal Virtual Taxa (VT) count table. Rows are VTs, columns are
samples. Produced by the gDAT pipeline from 18S rRNA SSU amplicons
(primers WANDA / AML2), with taxonomic assignment against the MaarjAM
database (97% identity, 95% alignment length).

## Usage

1. Clone the repository:

```bash
git clone https://github.com/Ebedthan/midas_microbial_succession_code.git
cd midas_microbial_succession_code
```

2. Open R and set the working directory to the repository root.

3. Run the alpha diversity analysis:

```r
source("02_R/richness_fit.R")
```

4. Run the beta diversity analysis:

```r
source("02_R/fit_beta.R")
```

Scripts are self-contained and will create `figures/` and `models/`
subdirectories automatically if they do not already exist.

## Citation

If you use this code or data, please cite the associated manuscript:

> Ebou A.E.T. *et al.* (2025). *Microbial succession in West African secondary forests: 
> rapid internal stabilisation without convergence toward old-growth
> reference states*

## Contact

- **Anicet E. T. Ebou** - ediman.ebou@inphb.ci  
  Equipe Bioinformatique et Biostatistiques,
  Laboratoire de Microbiologie, Biotechnologie et Bioinformatique,  
  Institut National Polytechnique Félix Houphouët-Boigny, Yamoussoukro,
  Côte d'Ivoire

- **Bruno Hérault** - bruno.herault@cirad.fr  
  CIRAD, UPR Forêts et Sociétés, Montpellier, France

## Licence

This repository is released under the [MIT Licence](LICENSE).
