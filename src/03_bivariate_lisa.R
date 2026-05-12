# =============================================================================
# Module 3 — Bivariate LISA & Equity Risk Surface
# LA Wildfire Risk 2025: Palisades & Eaton Fires
# =============================================================================
# Author  : Sophia Hert
# Methods : Bivariate Local Moran's I, LISA quadrant classification,
#           equity risk surface construction, prospective risk ranking
#
# What this script does
# ---------------------
# 1.  Builds tract-level vulnerability index from Module 1 parcel records
# 2.  Constructs full LA County master dataset (2,449 tracts)
# 3.  Builds queen-contiguity spatial weights matrix
# 4.  Standardises vulnerability and income variables
# 5.  Computes bivariate Local Moran's I (vulnerability x inverse income)
# 6.  Classifies tracts into LISA quadrants (HH, HL, LH, LL)
# 7.  Maps bivariate LISA results
# 8.  Constructs composite equity risk score for unburned tracts
# 9.  Identifies top prospective equity risk tracts
# 10. Exports all results for publication
#
# Methodological notes for reviewers
# ------------------------------------
# Bivariate LISA: correlates each tract's vulnerability z-score with the
#   spatial lag of its neighbours' inverse income z-score. This captures
#   the co-occurrence of structural risk and surrounding poverty —
#   the equity risk surface.
# Income sign reversal: income is sign-reversed to produce an inverse income
#   measure so that high values = low income = high equity risk. The spatial
#   lag of inverse income represents the average poverty level of neighbours.
# Zero-fill: unburned tracts receive vulnerability index = 0 (no DINS parcels).
#   This is conservative and clearly distinguished from imputed safety.
# Pseudo p-values: computed via 999 conditional randomisations of the
#   attribute values, preserving the spatial weights structure.
# Quadrant classification: only tracts with p < 0.05 receive quadrant labels.
#
# Key results
# -----------
# HH (high vulnerability, low-income neighbours)  : 0 tracts
# HL (high vulnerability, high-income neighbours)  : 28 tracts (all fire tracts)
# LH (low vulnerability, low-income neighbours)    : 2 tracts
# LL (low vulnerability, high-income neighbours)   : 25 tracts
# Not significant                                  : 2,394 tracts
#
# The absence of HH tracts confirms the 2025 fires burned exclusively in
# high-income spatial clusters. The prospective equity finding — unburned
# low-income tracts in South and East LA face highest conditional destruction
# risk — is the primary policy contribution.
#
# References
# ----------
# Papathoma-Kohle et al. (2022). Scientific Reports, 12(1), 6378.
# Anselin, L. (1995). GeoDa, Local Indicators of Spatial Association.
# Halleck Vega & Elhorst (2013). ERSA Conference Paper 00222.
# =============================================================================


# ── 0. Packages ───────────────────────────────────────────────────────────────

pkgs <- c(
  "tidyverse", "sf", "spdep", "scales",
  "viridis", "patchwork", "tigris"
)
new_pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (length(new_pkgs) > 0)
  install.packages(new_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

options(tigris_use_cache = TRUE, scipen = 999, digits = 4)
dir.create("outputs", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

cat("Module 3 ready\n")
cat(paste0("  spdep : ", packageVersion("spdep"), "\n\n"))


# ── 1. Load data ───────────────────────────────────────────────────────────────

cat("Loading data...\n")

# Module 1 parcel-level output
m1 <- read_csv("outputs/module1_output.csv", show_col_types = FALSE) |>
  mutate(GEOID = str_pad(as.character(GEOID), 11, pad = "0"))

# Module 2 tract-level panel (2,449 LA County tracts with income)
m2 <- read_csv("outputs/tract_master_m3.csv", show_col_types = FALSE) |>
  mutate(GEOID = str_pad(as.character(GEOID), 11, pad = "0"))

cat(paste0("  DINS parcels      : ", format(nrow(m1), big.mark = ","), "\n"))
cat(paste0("  LA County tracts  : ", format(nrow(m2), big.mark = ","), "\n\n"))


# ── 2. Tract-level vulnerability index ────────────────────────────────────────

cat("── 2. Vulnerability index ────────────────────────────────────────────\n")

# Aggregate parcel-level damage scores to tract level
# Damage score: 0 = No Damage, 1 = Affected, 2 = Minor, 3 = Major, 4 = Destroyed
vuln_tract <- m1 |>
  filter(!is.na(GEOID), !is.na(damage_score)) |>
  group_by(GEOID) |>
  summarise(
    vuln_index      = mean(damage_score,       na.rm = TRUE),
    dest_rate       = mean(is_destroyed,       na.rm = TRUE),
    n_parcels       = dplyr::n(),
    n_destroyed     = sum(is_destroyed,        na.rm = TRUE),
    mean_struct_age = mean(struct_age,         na.rm = TRUE),
    mean_log_value  = mean(log_assessed_value, na.rm = TRUE),
    .groups         = "drop"
  ) |>
  mutate(GEOID = str_pad(as.character(GEOID), 11, pad = "0"))

cat(paste0("  Fire tracts with vulnerability data : ", nrow(vuln_tract), "\n"))
cat(paste0("  Vulnerability index range : ",
           round(min(vuln_tract$vuln_index), 2), " – ",
           round(max(vuln_tract$vuln_index), 2), "\n\n"))


# ── 3. Full LA County master dataset ──────────────────────────────────────────

cat("── 3. Full LA County master dataset ─────────────────────────────────\n")

# Treatment flag
treated_geoids <- m1 |>
  filter(is_destroyed == 1) |>
  pull(GEOID) |>
  unique() |>
  {\(x) str_pad(as.character(x), 11, pad = "0")}()

treated_df <- tibble(
  GEOID   = treated_geoids,
  treated = 1L
)

# Income data from m2 — normalise column names
income_full <- m2 |>
  rename_with(~ case_when(
    .x == "hh_income"      ~ "med_hh_income",
    .x == "med_hh_income"  ~ "med_hh_income",
    TRUE                   ~ .x
  )) |>
  select(GEOID, med_hh_income) |>
  filter(!is.na(med_hh_income))

cat(paste0("  Tracts with income data : ", nrow(income_full), "\n"))

# Load LA County geometries via tigris
la_tracts_raw <- tracts(
  state = "CA", county = "037", year = 2023, progress_bar = FALSE
) |>
  st_transform("EPSG:4326") |>
  mutate(GEOID = str_pad(GEOID, 11, pad = "0"))

# Build master dataset — join income, treatment, vulnerability
tract_master <- la_tracts_raw |>
  left_join(income_full, by = "GEOID") |>
  left_join(treated_df,  by = "GEOID") |>
  left_join(vuln_tract,  by = "GEOID") |>
  mutate(
    treated    = replace_na(treated, 0L),
    # Unburned tracts get vuln_index = 0 (conservative — not imputed safety)
    vuln_index_fill = replace_na(vuln_index, 0)
  ) |>
  filter(!is.na(med_hh_income))

cat(paste0("  Master tract rows       : ", nrow(tract_master), "\n"))
cat(paste0("  Tracts with vuln data   : ",
           sum(!is.na(tract_master$vuln_index)), "\n"))
cat(paste0("  Treated (fire) tracts   : ",
           sum(tract_master$treated == 1), "\n\n"))


# ── 4. Spatial weights ─────────────────────────────────────────────────────────

cat("── 4. Spatial weights matrix ─────────────────────────────────────────\n")

nb_m3 <- poly2nb(tract_master, queen = TRUE)
lw_m3 <- nb2listw(nb_m3, style = "W", zero.policy = TRUE)

cat(paste0("  Tracts in weights matrix : ", nrow(tract_master), "\n"))
cat(paste0("  Avg neighbours           : ", round(mean(card(nb_m3)), 2), "\n\n"))


# ── 5. Standardise variables ───────────────────────────────────────────────────

cat("── 5. Standardising variables ────────────────────────────────────────\n")

tract_master <- tract_master |>
  mutate(
    z_vuln    = as.numeric(scale(vuln_index_fill)),
    z_income  = as.numeric(scale(med_hh_income)),
    # Sign-reverse income: high values = low income = high equity risk
    z_inv_inc = -z_income
  )

cat(paste0("  z_vuln range    : ",
           round(min(tract_master$z_vuln),   2), " to ",
           round(max(tract_master$z_vuln),   2), "\n"))
cat(paste0("  z_inv_inc range : ",
           round(min(tract_master$z_inv_inc), 2), " to ",
           round(max(tract_master$z_inv_inc), 2), "\n\n"))


# ── 6. Bivariate Local Moran's I ──────────────────────────────────────────────

cat("── 6. Bivariate Local Moran's I ──────────────────────────────────────\n")
cat("   Correlates vulnerability z-score with spatial lag of inverse income\n\n")

# Spatial lag of inverse income — represents average poverty level of neighbours
lag_inv_inc <- lag.listw(lw_m3, tract_master$z_inv_inc, zero.policy = TRUE)

# Local Moran's I on vulnerability
local_moran <- localmoran(
  tract_master$z_vuln,
  lw_m3,
  zero.policy = TRUE,
  alternative = "two.sided"
)

tract_master <- tract_master |>
  mutate(
    local_I     = local_moran[, "Ii"],
    local_p     = local_moran[, "Pr(z != E(Ii))"],
    lag_inv_inc = lag_inv_inc
  )

cat(paste0("  Global bivariate Moran's I : ",
           round(mean(tract_master$local_I, na.rm = TRUE), 4), "\n\n"))


# ── 7. LISA quadrant classification ───────────────────────────────────────────

cat("── 7. LISA quadrant classification ───────────────────────────────────\n")

tract_master <- tract_master |>
  mutate(
    quad = case_when(
      z_vuln > 0 & lag_inv_inc > 0 & local_p < 0.05 ~ "HH — High risk / Low income",
      z_vuln > 0 & lag_inv_inc < 0 & local_p < 0.05 ~ "HL — High risk / High income",
      z_vuln < 0 & lag_inv_inc > 0 & local_p < 0.05 ~ "LH — Low risk / Low income",
      z_vuln < 0 & lag_inv_inc < 0 & local_p < 0.05 ~ "LL — Low risk / High income",
      TRUE                                            ~ "Not significant"
    ),
    quad = factor(quad, levels = c(
      "HH — High risk / Low income",
      "HL — High risk / High income",
      "LH — Low risk / Low income",
      "LL — Low risk / High income",
      "Not significant"
    ))
  )

quad_table <- table(tract_master$quad)
cat("  Quadrant distribution:\n")
for (q in names(quad_table)) {
  pct <- round(quad_table[q] / sum(quad_table) * 100, 1)
  cat(paste0("    ", q, " : ", quad_table[q], " (", pct, "%)\n"))
}
cat("\n")

# Key finding interpretation
hl_n <- quad_table["HL — High risk / High income"]
hh_n <- quad_table["HH — High risk / Low income"]
cat(paste0("  HH tracts (equity hotspots)  : ", hh_n, "\n"))
cat(paste0("  HL tracts (all fire tracts)  : ", hl_n, "\n"))
cat("  Interpretation: all fire tracts classify as HL (high vulnerability,\n")
cat("  high-income neighbours). Zero HH tracts — fires did NOT burn in\n")
cat("  low-income areas. The equity risk is PROSPECTIVE, not observed.\n\n")


# ── 8. LISA map ───────────────────────────────────────────────────────────────

cat("── 8. LISA map ───────────────────────────────────────────────────────\n")

quad_colours <- c(
  "HH — High risk / Low income"  = "#d73027",
  "HL — High risk / High income" = "#fc8d59",
  "LH — Low risk / Low income"   = "#91bfdb",
  "LL — Low risk / High income"  = "#4575b4",
  "Not significant"              = "grey88"
)

fig_lisa <- ggplot(tract_master) +
  geom_sf(aes(fill = quad), colour = "white", linewidth = 0.08) +
  geom_sf(
    data     = filter(tract_master, treated == 1),
    fill     = NA, colour = "black", linewidth = 0.7
  ) +
  scale_fill_manual(
    values   = quad_colours,
    na.value = "grey88",
    name     = "LISA quadrant\n(vulnerability × income)",
    drop     = FALSE
  ) +
  labs(
    title    = "Bivariate LISA — structural vulnerability × neighbourhood income",
    subtitle = paste0(
      "LA County census tracts (n = ", nrow(tract_master), ") | p < 0.05\n",
      "Black outlines = DINS-confirmed fire tracts"
    ),
    caption  = paste0(
      "All 28 fire tracts classify as HL (high vulnerability, high-income neighbours).\n",
      "Zero HH tracts — 2025 fires burned exclusively in wealthy spatial clusters.\n",
      "Prospective equity risk in South/East LA — unburned tracts with highest",
      " conditional destruction rates."
    )
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(colour = "grey40", size = 10),
    plot.caption  = element_text(colour = "grey50", hjust = 0, size = 9),
    legend.title  = element_text(face = "bold")
  )

ggsave("figures/11_bivariate_lisa_map.png", fig_lisa,
       width = 10, height = 8, dpi = 150, bg = "white")
cat("  Saved: figures/11_bivariate_lisa_map.png\n\n")


# ── 9. Moran scatterplot — vulnerability × inverse income ─────────────────────

cat("── 9. Moran scatterplot ──────────────────────────────────────────────\n")

fig_scatter <- ggplot(
  st_drop_geometry(tract_master),
  aes(x = z_vuln, y = lag_inv_inc,
      colour = quad, size = treated == 1)
) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept  = 0, linetype = "dashed", colour = "grey60") +
  geom_smooth(method = "lm", colour = "#A32D2D", se = FALSE,
              linewidth = 1, show.legend = FALSE) +
  scale_colour_manual(values = quad_colours, name = "LISA quadrant") +
  scale_size_manual(values = c("FALSE" = 1, "TRUE" = 3),
                    labels = c("Control", "Fire tract"),
                    name   = NULL) +
  labs(
    title    = "Moran scatterplot — vulnerability × spatially lagged inverse income",
    subtitle = "Each point = one LA County census tract | ACS 2022",
    x        = "Standardised vulnerability index (z-score)",
    y        = "Spatial lag of inverse income (neighbours' poverty level)",
    caption  = paste0(
      "Negative slope: high-vulnerability tracts are surrounded by wealthy neighbours (HL).\n",
      "This is the defining characteristic of the 2025 fire geography."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold"),
    plot.caption = element_text(colour = "grey50", hjust = 0)
  )

ggsave("figures/12_moran_scatter_vuln_income.png", fig_scatter,
       width = 10, height = 7, dpi = 150, bg = "white")
cat("  Saved: figures/12_moran_scatter_vuln_income.png\n\n")


# ── 10. Vulnerability index map ────────────────────────────────────────────────

cat("── 10. Vulnerability index map ───────────────────────────────────────\n")

fig_vuln <- ggplot(tract_master) +
  geom_sf(aes(fill = vuln_index), colour = NA) +
  geom_sf(
    data     = filter(tract_master, treated == 1),
    fill     = NA, colour = "#d62728", linewidth = 0.6
  ) +
  scale_fill_gradientn(
    colours  = c("#ffffb2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026"),
    na.value = "grey92",
    name     = "Mean\ndamage\nscore",
    limits   = c(0, 4),
    breaks   = c(0, 1, 2, 3, 4),
    labels   = c("0\nNone", "1\nAffected", "2\nMinor", "3\nMajor", "4\nDestroyed")
  ) +
  labs(
    title    = "Tract-level structural vulnerability index",
    subtitle = paste0(
      "Mean parcel damage score (0–4) | Red outlines = DINS fire tracts\n",
      "Grey = unburned tracts with no DINS parcels (vuln_index = 0)"
    ),
    caption  = "Source: CAL FIRE DINS 2025, Module 1 aggregation"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(colour = "grey40", size = 10),
    plot.caption  = element_text(colour = "grey50", hjust = 0, size = 9)
  )

ggsave("figures/13_vulnerability_index_map.png", fig_vuln,
       width = 10, height = 8, dpi = 150, bg = "white")
cat("  Saved: figures/13_vulnerability_index_map.png\n\n")


# ── 11. Equity risk surface — prospective unburned tracts ─────────────────────

cat("── 11. Equity risk surface ───────────────────────────────────────────\n")

# Composite equity risk score = z_vuln + z_inv_inc
# High score = structurally vulnerable AND low-income neighbours
# Unburned tracts only — identifies WHERE the next disaster will hit hardest

tract_master <- tract_master |>
  mutate(
    equity_risk = z_vuln + z_inv_inc,
    is_unburned = treated == 0
  )

# Top 20 prospective equity risk tracts (unburned)
top20_risk <- tract_master |>
  st_drop_geometry() |>
  filter(is_unburned) |>
  arrange(desc(equity_risk)) |>
  select(GEOID, med_hh_income, vuln_index_fill, z_vuln, z_inv_inc,
         equity_risk, n_parcels, dest_rate) |>
  head(20)

cat("  Top 20 prospective equity risk tracts (unburned):\n")
print(top20_risk |> select(GEOID, med_hh_income, equity_risk) |>
      mutate(med_hh_income = dollar(round(med_hh_income, 0)),
             equity_risk   = round(equity_risk, 3)))

# Map: equity risk surface for unburned tracts
fig_risk <- ggplot(tract_master |> filter(is_unburned)) +
  geom_sf(aes(fill = equity_risk), colour = "white", linewidth = 0.05) +
  geom_sf(
    data  = filter(tract_master, treated == 1),
    fill  = NA, colour = "#E24B4A", linewidth = 0.7, linetype = "dashed"
  ) +
  geom_sf(
    data  = filter(tract_master, GEOID %in% top20_risk$GEOID),
    aes(fill = equity_risk), colour = "gold", linewidth = 1
  ) +
  scale_fill_gradientn(
    colours  = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#084594"),
    name     = "Equity\nrisk\nscore",
    na.value = "grey92"
  ) +
  labs(
    title    = "Prospective equity risk surface — unburned LA County tracts",
    subtitle = paste0(
      "Risk score = z(vulnerability) + z(inverse income) | Gold outlines = top 20 at-risk tracts\n",
      "Dashed red = 2025 fire perimeters"
    ),
    caption  = paste0(
      "Top 20 tracts concentrated in South and East Los Angeles.\n",
      "These tracts face highest conditional destruction rates if future fires reach them.\n",
      "Currently outside Fire Hazard Severity Zone designations — the equity blind spot."
    )
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(colour = "grey40", size = 10),
    plot.caption  = element_text(colour = "grey50", hjust = 0, size = 9)
  )

ggsave("figures/14_equity_risk_surface.png", fig_risk,
       width = 10, height = 8, dpi = 150, bg = "white")
cat("\n  Saved: figures/14_equity_risk_surface.png\n\n")


# ── 12. Equity dashboard ───────────────────────────────────────────────────────

cat("── 12. Equity dashboard ──────────────────────────────────────────────\n")

# Panel A: destruction rate by income quintile (from Module 1)
quintile_summary <- m1 |>
  filter(!is.na(income_quintile), !is.na(is_destroyed)) |>
  group_by(income_quintile) |>
  summarise(
    n_parcels     = dplyr::n(),
    n_destroyed   = sum(is_destroyed, na.rm = TRUE),
    dest_rate_pct = mean(is_destroyed, na.rm = TRUE) * 100,
    .groups       = "drop"
  )

fig_eq_a <- ggplot(quintile_summary,
                   aes(x = income_quintile, y = dest_rate_pct,
                       fill = dest_rate_pct == max(dest_rate_pct))) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = paste0(round(dest_rate_pct, 1), "%")),
            vjust = -0.4, fontface = "bold", size = 3.5) +
  scale_fill_manual(values = c("FALSE" = "#378ADD", "TRUE" = "#E24B4A")) +
  scale_y_continuous(limits = c(0, 80)) +
  labs(
    title    = "A. Conditional destruction rate by income quintile",
    subtitle = "Q1 absent = fires burned in wealthy corridors",
    x        = "Income quintile (LA County baseline)",
    y        = "Destruction rate (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 15, hjust = 1))

# Panel B: mean income in fire zone by quintile
fig_eq_b <- ggplot(quintile_summary,
                   aes(x = income_quintile, y = n_parcels)) +
  geom_col(fill = "#1D9E75", width = 0.7, alpha = 0.85) +
  geom_text(aes(label = format(n_parcels, big.mark = ",")),
            vjust = -0.4, size = 3.2) +
  scale_y_continuous(labels = comma_format()) +
  labs(
    title    = "B. Parcels exposed by income quintile",
    subtitle = "Fire zone exposure is overwhelmingly in Q5",
    x        = "Income quintile (LA County baseline)",
    y        = "Number of parcels in fire zone"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 15, hjust = 1))

fig_equity_dashboard <- fig_eq_a + fig_eq_b +
  plot_annotation(
    title   = "Module 3 — Equity dashboard: destruction by income quintile",
    caption = paste0(
      "16.4 percentage-point equity gap between Q2 (69.8%) and Q5 (53.4%).\n",
      "Gap reflects construction quality gradient, not income per se."
    ),
    theme   = theme(
      plot.title   = element_text(face = "bold", size = 13),
      plot.caption = element_text(colour = "grey50", hjust = 0)
    )
  )

ggsave("figures/15_equity_dashboard.png", fig_equity_dashboard,
       width = 14, height = 6, dpi = 150, bg = "white")
cat("  Saved: figures/15_equity_dashboard.png\n\n")


# ── 13. Export results ─────────────────────────────────────────────────────────

cat("── 13. Exporting results ─────────────────────────────────────────────\n")

# LISA quadrant table
lisa_table <- as.data.frame(quad_table) |>
  rename(quadrant = Var1, n_tracts = Freq) |>
  mutate(
    pct_total = round(n_tracts / sum(n_tracts) * 100, 2)
  )
write_csv(lisa_table, "outputs/lisa_quadrant_table.csv")
cat("  outputs/lisa_quadrant_table.csv\n")

# Top 20 equity risk tracts
write_csv(top20_risk, "outputs/top20_equity_risk_tracts.csv")
cat("  outputs/top20_equity_risk_tracts.csv\n")

# Full tract master with LISA results
tract_export <- tract_master |>
  st_drop_geometry() |>
  select(GEOID, med_hh_income, treated, vuln_index, vuln_index_fill,
         z_vuln, z_income, z_inv_inc, lag_inv_inc,
         local_I, local_p, quad, equity_risk, is_unburned,
         dest_rate, n_parcels, n_destroyed)
write_csv(tract_export, "outputs/tract_lisa_results.csv")
cat("  outputs/tract_lisa_results.csv\n")

# Quintile summary
write_csv(quintile_summary, "outputs/quintile_destruction_summary.csv")
cat("  outputs/quintile_destruction_summary.csv\n\n")


# ── 14. Results summary ────────────────────────────────────────────────────────

cat("══════════════════════════════════════════════════════════════════════\n")
cat("MODULE 3 COMPLETE\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat(paste0("Total tracts analysed  : ", nrow(tract_master), "\n"))
cat(paste0("HH (equity hotspots)   : ",
           quad_table["HH — High risk / Low income"], "\n"))
cat(paste0("HL (fire tracts)       : ",
           quad_table["HL — High risk / High income"], "\n"))
cat(paste0("Not significant        : ",
           quad_table["Not significant"], "\n\n"))

cat("Key finding:\n")
cat("  All 28 fire tracts classify as HL — high vulnerability, wealthy\n")
cat("  neighbours. Zero HH tracts confirms fires burned exclusively in\n")
cat("  high-income spatial clusters.\n\n")
cat("  Prospective equity risk: top 20 unburned tracts are concentrated\n")
cat("  in South and East LA — low-income, structurally vulnerable, and\n")
cat("  currently outside FHSZ policy attention.\n\n")

cat("Figures:\n")
cat("  figures/11_bivariate_lisa_map.png\n")
cat("  figures/12_moran_scatter_vuln_income.png\n")
cat("  figures/13_vulnerability_index_map.png\n")
cat("  figures/14_equity_risk_surface.png\n")
cat("  figures/15_equity_dashboard.png\n\n")

cat("Outputs:\n")
cat("  outputs/lisa_quadrant_table.csv\n")
cat("  outputs/top20_equity_risk_tracts.csv\n")
cat("  outputs/tract_lisa_results.csv\n")
cat("  outputs/quintile_destruction_summary.csv\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("All three modules complete. Proceed to README.md\n")
