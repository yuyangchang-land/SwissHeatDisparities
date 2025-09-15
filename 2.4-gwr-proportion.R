# ================================
# ================================
# Bar chart of exposure categories from GWR t-values
# (with socio + physical variables)
# ================================
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# ---- Inputs ----
base_dir <- "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/revision_r1/results_gwr_multi"
run      <- "with_socio"  # change to "no_socio" for robustness

# Socio variables (reversed inside group)
vars_socio <- rev(c("income_ln","age65_79","age80_","percen_eu","percen_noneu","socioassis","ave_living"))

# Physical variables (reversed inside group)
vars_phys  <- rev(c("dem_mean","evi_mean","bldg_mean","albd_mean","imp_mean"))

# Full order: socio first, physical after
vars_keep <- c(vars_phys, vars_socio)

# Which DVs to plot (includes CHEI)
dvs <- c("CHEI_pca", "LST_z", "HWD_z", "HWP_z")

# Facet labels
dv_labels <- c(
  CHEI_pca = "CHEI",
  LST_z    = "Land surface temperature",
  HWD_z    = "Heat warning days",
  HWP_z    = "Heatwave probability"
)

# Variable labels (match keys above)
var_labels <- c(
  # socio
  income_ln    = "Income (log)",
  age65_79     = "Age 65–79 (%)",
  age80_       = "Age 80+ (%)",
  percen_eu    = "EU foreigners (%)",
  percen_noneu = "Non–EU foreigners (%)",
  socioassis   = "Social assistance (%)",
  ave_living   = "Avg. living space",
  # physical
  dem_mean     = "Elevation",
  evi_mean     = "EVI",
  bldg_mean    = "Building height",
  albd_mean    = "Albedo",
  imp_mean     = "Imperviousness"
)

# ---- Helper to read one DV's local results (coeffs + t-values) ----
read_local <- function(dv) {
  path <- file.path(base_dir, paste0(run, "_", dv), paste0(dv, "_gwr_results.csv"))
  if (!file.exists(path)) stop("Not found: ", path)
  read_csv(path, show_col_types = FALSE)
}

# ---- Build long summary of categories by DV & variable ----
df_summary <- lapply(dvs, function(dv) {
  df <- read_local(dv)
  
  # Only keep variables that actually exist in this results file
  present <- intersect(vars_keep, names(df))
  tcols   <- paste0(present, "_t")
  needed  <- intersect(c(present, tcols), names(df))
  if (!length(needed)) return(NULL)
  
  dd <- df[, needed]
  
  # Long format: one row per municipality x variable
  long <- bind_rows(lapply(present, function(v) {
    tcol <- paste0(v, "_t")
    tibble(
      DV        = dv,
      Variable  = v,
      t_value   = if (tcol %in% names(dd)) dd[[tcol]] else NA_real_
    )
  }))
  
  # Classify exposure category from t-values
  long <- long %>%
    mutate(Exposure_Category = case_when(
      t_value >  1.96 ~ "Overexposed",
      t_value < -1.96 ~ "Underexposed",
      TRUE            ~ "No significant disparity"
    ))
  
  # Proportions per variable (ignore NAs)
  long %>%
    filter(!is.na(t_value)) %>%
    count(DV, Variable, Exposure_Category, name = "n") %>%
    group_by(DV, Variable) %>%
    mutate(Total = sum(n),
           Proportion = n / Total) %>%
    ungroup()
}) %>% bind_rows()

# ---- Relabel + enforce order ----
# Drop label keys that aren't present to avoid NA labels
vars_present_all <- intersect(vars_keep, unique(df_summary$Variable))

df_summary <- df_summary %>%
  mutate(
    Dataset  = factor(DV, levels = dvs, labels = dv_labels[dvs]),
    Variable = factor(Variable,
                      levels = vars_present_all,                 # keep intended order
                      labels = var_labels[vars_present_all]),
    Exposure_Category = factor(Exposure_Category,
                               levels = c("Underexposed","No significant disparity","Overexposed"))
  )

# ---- Plot ----
p <- ggplot(df_summary, aes(x = Variable, y = Proportion, fill = Exposure_Category)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Dataset, nrow = 1) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Underexposed" = "#536976",
      "No significant disparity" = "gray95",
      "Overexposed" = "#f5af19"
    )
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = NULL,
    y = "Proportion of municipalities (%)",
    fill = NULL,
    title = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

print(p)

# Optional save
# ggsave(file.path(base_dir, paste0("bar_exposure_categories_", run, "_with_physical.png")),
#        p, width = 13, height = 6, dpi = 300)