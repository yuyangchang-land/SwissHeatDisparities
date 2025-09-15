library(dplyr)
library(ggplot2)
library(forcats)
library(tibble)
library(patchwork)
# ---- 0) Path ----
data_path <- "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/revision_r1/muni_all_vars_final.csv"

# ---- 1) Load ----
dat <- read_csv(data_path, show_col_types = FALSE)
colnames(dat)
# Helper: z-standardize
std <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

# ---- 2) Replace "X" with NA in character cols, then coerce needed columns to numeric ----
raw_needed <- c(
  # outcomes
  "lst_rbn","avgt25_","htwv_p_",
  # physical drivers
  "dem_men","imprv_m","evi_men","bhght_m","albd_mn",
  # socio-demographics
  "income","age65_","age80_","forgn19","for_eu","socisss","av_lvng",
  # weights
  "num_rsd"
)

dat <- dat %>%
  mutate(across(where(is.character), ~na_if(., "X"))) %>%           # "X" -> NA safely
  mutate(across(all_of(raw_needed), ~suppressWarnings(as.numeric(.))))

# ---- 3) Derive socio fields + weights ----
dat <- dat %>%
  mutate(
    # EU share as proportion
    for_eu_prop = ifelse(for_eu > 1, for_eu/100, for_eu),
    # age split
    age65_79    = age65_ - age80_,
    # foreigners split
    percen_eu    = forgn19 * for_eu_prop,
    percen_noneu = forgn19 * (1 - for_eu_prop),
    # weights = log(pop + 1)
    w = log(pmax(num_rsd, 0) + 1),
    # log-income
    income_ln = log(pmax(income, 0) + 1)
  )

# ---- 4) Build composite heat index (CHEI_pca) from LST/HWD/HWP ----
dat <- dat %>%
  mutate(
    LST_z = std(lst_rbn),
    HWD_z = std(avgt25_),
    HWP_z = std(htwv_p_)
  )

X3 <- dat %>% select(LST_z, HWD_z, HWP_z)

# keep only rows with valid data
ok <- complete.cases(X3) & apply(as.matrix(X3), 1, function(r) all(is.finite(r)))
stopifnot(sum(ok) >= 3L)  # enough rows for PCA

# PCA
pca_fit <- prcomp(as.matrix(X3[ok, ]), center = TRUE, scale. = FALSE)

# 1) Equal weights
w_equal <- rep(1/3, 3)
dat$CHEI_equal_raw <- as.numeric(as.matrix(X3) %*% w_equal)

# scale to [0,1]
rng1 <- range(dat$CHEI_equal_raw, na.rm = TRUE)
dat$CHEI_equal <- (dat$CHEI_equal_raw - rng1[1]) / (rng1[2] - rng1[1])

# 2) PCA weights (first PC loadings normalized to sum = 1)
w_pca <- pca_fit$rotation[,1] / sum(pca_fit$rotation[,1])

dat$CHEI_pca_raw <- as.numeric(as.matrix(X3) %*% w_pca)

# scale to [0,1]
rng2 <- range(dat$CHEI_pca_raw, na.rm = TRUE)
dat$CHEI_pca <- (dat$CHEI_pca_raw - rng2[1]) / (rng2[2] - rng2[1])

# Check weights
cat("Equal weights (LST, HWD, HWP):", round(w_equal, 3), "\n")
cat("PCA weights (LST, HWD, HWP):", round(w_pca, 3), "\n")

# ---- 5) Final drop of incomplete rows for modeling ----
physical_vars <- c("dem_men","imprv_m","evi_men","bhght_m","albd_mn")
socio_vars    <- c("income_ln","age65_79","age80_","percen_eu","percen_noneu","socisss","av_lvng")

need_for_models <- c("CHEI_pca", physical_vars, socio_vars, "w", "num_rsd")

dat <- dat %>%
  mutate(across(all_of(need_for_models), ~ifelse(is.infinite(.), NA, .))) %>%
  tidyr::drop_na(all_of(need_for_models))

# ---- 6) Quick sanity prints ----
cat("Rows remaining after cleaning:", nrow(dat), "\n")
cat("Equal weights (LST, HWD, HWP):", round(w_pca, 3), "\n")



library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(Hmisc)

# -----------------------------
# 1) Settings
# -----------------------------
socio_vars <- c("income_ln","age65_79","age80_",
                "percen_eu","percen_noneu","socisss","av_lvng")

var_labels <- c(
  income_ln    = "Income (log)",
  age65_79     = "Age 65–79",
  age80_       = "Age 80+",
  percen_eu    = "EU foreigners",
  percen_noneu = "Non-EU foreigners",
  socisss   = "Social assistance",
  av_lvng   = "Living space per capita"
)

# Heat indicators to plot (left-to-right order)
indicator_vars <- c("CHEI_equal", "lst_rbn", "avgt25_", "htwv_p_")
indicator_labels <- c(
  CHEI_equal  = "CHEI",
  lst_urban = "Land surface temperature",
  heat_d    = "Heat warning days",
  heat_p    = "Heatwave probability"
)

# -----------------------------
# 2) Make quartiles and long data
# -----------------------------
dat_q <- dat %>%
  mutate(across(all_of(socio_vars), ~ntile(., 4), .names = "{.col}_q"))

# Long socio (quartiles) + long indicators
plot_df <- dat_q %>%
  pivot_longer(ends_with("_q"), names_to = "socio_var", values_to = "quartile") %>%
  mutate(socio_var = sub("_q$", "", socio_var)) %>%
  pivot_longer(all_of(indicator_vars), names_to = "indicator", values_to = "heat_value")

# Summary per socio_var × quartile × indicator
summary_df <- plot_df %>%
  group_by(socio_var, quartile, indicator) %>%
  summarise(
    mean_heat = mean(heat_value, na.rm = TRUE),
    se_heat   = sd(heat_value,   na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# -----------------------------
# 3) Build y-axis order and shading
# -----------------------------
# Exact y-axis order: for each socio var, Q1..Q4
y_levels <- unlist(lapply(socio_vars, function(v) paste(var_labels[v], paste0("Q", 1:4))),
                   use.names = FALSE)

summary_df <- summary_df %>%
  mutate(
    socio_lab = var_labels[socio_var],
    indicator = factor(indicator, levels = indicator_vars, labels = indicator_labels),
    y_axis    = paste(socio_lab, paste0("Q", quartile)),
    y_axis    = factor(y_axis, levels = rev(y_levels))   # reverse so Q1 at top
  )

# Shaded blocks (one 4-row band per socio var), using the reversed levels used for plotting
y_levels_rev <- rev(y_levels)
shade_df <- tibble(socio_var = socio_vars) %>%
  mutate(
    socio_lab  = var_labels[socio_var],
    block_rows = lapply(socio_lab, function(lbl) which(y_levels_rev %in% paste(lbl, paste0("Q", 1:4)))),
    ymin       = sapply(block_rows, min) - 0.5,
    ymax       = sapply(block_rows, max) + 0.5,
    fill_col   = rep(c("gray95","white"), length.out = length(socio_vars))
  )

# -----------------------------
# 4) Plot all indicators together
# -----------------------------
p_all <- ggplot(summary_df, aes(x = mean_heat, y = y_axis)) +
  # background shading
  geom_rect(data = shade_df,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = fill_col),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_fill_identity() +
  # points + SE bars
  geom_point(size = 1.5, color = "brown") +
  geom_errorbarh(
    aes(xmin = mean_heat - se_heat, xmax = mean_heat + se_heat),
    height = 0,  # removes vertical short lines
    color = "brown"
  )+
  facet_wrap(~ indicator, ncol = length(indicator_vars), scales = "free_x") +
  scale_color_brewer(palette = "Dark2", name = "Quartile") + 
  labs(
    x = NULL,  # each facet has different scale; keep x label off or set a generic label
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(hjust = 1),
    strip.text = element_text(face = "bold")
  )

# Show it
print(p_all)



#########weight population

# =============================
# Libraries
# =============================
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(Hmisc)   # wtd.mean, wtd.var

# =============================
# 0) Settings (edit as needed)
# =============================
socio_vars <- c("income_ln","age65_79","age80_",
                "percen_eu","percen_noneu","socisss","av_lvng")

var_labels <- c(
  income_ln    = "Income (log)",
  age65_79     = "Age 65–79",
  age80_       = "Age 80+",
  percen_eu    = "EU foreigners",
  percen_noneu = "Non-EU foreigners",
  socisss      = "Social assistance",
  av_lvng      = "Living space per capita"
)

# Heat indicators (in facet order) and pretty labels
indicator_vars <- c("CHEI_equal", "lst_rbn", "avgt25_", "htwv_p_")
indicator_labels <- c(
  CHEI_equal = "CHEI",
  lst_rbn    = "Land surface temperature",
  avgt25_    = "Heat warning days",
  htwv_p_    = "Heatwave probability"
)

# =============================
# 1) Make log(pop) weights and an ID
# =============================
# Assumes dat already exists in your workspace and has 'num_rsd' (resident population)
dat <- dat %>%
  mutate(
    muni_id  = row_number(),
    num_rsd  = ifelse(is.na(num_rsd) | num_rsd <= 0, NA_real_, num_rsd),
    w_logpop = log(pmax(num_rsd, 1))   # guard against zeros
  )

# =============================
# 2) Quartiles and long data
# =============================
dat_q <- dat %>%
  mutate(across(all_of(socio_vars), ~ntile(., 4), .names = "{.col}_q"))

plot_df <- dat_q %>%
  # keep ID, weights, socio vars, their quartiles, and indicators
  select(muni_id, w_logpop, all_of(socio_vars), ends_with("_q"), all_of(indicator_vars)) %>%
  # socio quartiles to long
  pivot_longer(ends_with("_q"), names_to = "socio_var", values_to = "quartile") %>%
  mutate(socio_var = sub("_q$", "", socio_var)) %>%
  # indicators to long
  pivot_longer(all_of(indicator_vars), names_to = "indicator", values_to = "heat_value")

# =============================
# 3) Weighted summaries
# =============================
weighted_se_mean <- function(x, w){
  ok <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w)
  if (!any(ok)) return(NA_real_)
  var_w <- Hmisc::wtd.var(x[ok], weights = w[ok], na.rm = TRUE)
  sw    <- sum(w[ok])
  sw2   <- sum(w[ok]^2)
  Neff  <- ifelse(sw2 > 0, (sw^2) / sw2, NA_real_)
  sqrt(var_w / Neff)
}

#summary_df <- plot_df %>%
#  group_by(socio_var, quartile, indicator) %>%
#  summarise(
#    mean_heat = Hmisc::wtd.mean(heat_value, weights = w_logpop, na.rm = TRUE),
#    se_heat   = weighted_se_mean(heat_value, w_logpop),
#    .groups   = "drop"
#  )

summary_df <- plot_df %>%
  group_by(socio_var, quartile, indicator) %>%
  summarise(
    mean_heat = Hmisc::wtd.mean(heat_value, weights = dat$num_rsd[muni_id], na.rm = TRUE),
    se_heat   = weighted_se_mean(heat_value, dat$num_rsd[muni_id]),
    .groups   = "drop"
  )
# =============================
# 4) Y-axis order and shading bands
# =============================
y_levels <- unlist(lapply(socio_vars, function(v) paste(var_labels[v], paste0("Q", 1:4))),
                   use.names = FALSE)

summary_df <- summary_df %>%
  mutate(
    socio_lab = var_labels[socio_var],
    indicator = factor(indicator, levels = indicator_vars,
                       labels = indicator_labels[indicator_vars]),
    y_axis    = paste(socio_lab, paste0("Q", quartile)),
    y_axis    = factor(y_axis, levels = rev(y_levels))   # reverse so Q1 at top
  )

y_levels_rev <- rev(y_levels)
shade_df <- tibble(socio_var = socio_vars) %>%
  mutate(
    socio_lab  = var_labels[socio_var],
    block_rows = lapply(socio_lab, function(lbl) which(y_levels_rev %in% paste(lbl, paste0("Q", 1:4)))),
    ymin       = sapply(block_rows, min) - 0.5,
    ymax       = sapply(block_rows, max) + 0.5,
    fill_col   = rep(c("gray95","white"), length.out = length(socio_vars))
  )

# =============================
# 5) Plot
# =============================
p_all <- ggplot(summary_df, aes(x = mean_heat, y = y_axis)) +
  geom_rect(data = shade_df,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = fill_col),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_fill_identity() +
  geom_point(size = 1.5, color = "brown") +
  geom_errorbarh(
    aes(xmin = mean_heat - se_heat, xmax = mean_heat + se_heat),
    height = 0,
    color = "brown"
  ) +
  facet_wrap(~ indicator, ncol = length(indicator_vars), scales = "free_x") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y       = element_text(hjust = 1),
    strip.text        = element_text(face = "bold")
  )

print(p_all)

# Save high-res
ggsave("quartile_heat_exposure_weighted_logpop.png", p_all, width = 12, height = 7, dpi = 300)

# =============================
# 6) Caption helper (paste to manuscript)
# =============================
cat(
  "Figure caption suggestion:\n",
  "Quartile differences in average heat exposure across socio-demographic groups.\n",
  "Municipalities are divided into quartiles (Q1–Q4) for each socio-demographic variable,\n",
  "where Q1 is the lowest and Q4 the highest. Points show means weighted by the log of\n",
  "municipal population; horizontal bars show ±1 standard error of the weighted mean.\n",
  "Facets report CHEI, land surface temperature, heat-warning days, and heatwave probability.\n"
)

# 1) Confirm weights vary
summary(dat$w_logpop)           # should not be constant
length(unique(dat$w_logpop))    # should be > 1

# 2) Compare unweighted vs weighted means for one facet
tmp <- subset(plot_df, indicator == "CHEI_equal" & socio_var == "income_ln" & quartile == 4)
mean(tmp$heat_value, na.rm = TRUE)                           # unweighted
Hmisc::wtd.mean(tmp$heat_value, weights = tmp$w_logpop)      # weighted by log(pop)
Hmisc::wtd.mean(tmp$heat_value, weights = dat$num_rsd[tmp$muni_id])  # raw pop (stronger effect)

# 3) Sanity: effective sample size used in SE (should be < n when weights unequal)
weighted_se_mean(tmp$heat_value, tmp$w_logpop)


##########roubtness check for CHEI pca and equal
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(Hmisc)   # wtd.mean, wtd.var

# -----------------------------
# 0) Ensure weights are ready
# -----------------------------
dat <- dat %>%
  mutate(
    muni_id = row_number(),
    num_rsd = ifelse(is.na(num_rsd) | num_rsd <= 0, NA_real_, num_rsd)  # guard against zeros
  )

# -----------------------------
# 1) Heat indicators: only CHEI_equal and CHEI_pca
# -----------------------------
indicator_vars <- c("CHEI_equal", "CHEI_pca")
indicator_labels <- c(
  CHEI_equal = "CHEI (Equal weights)",
  CHEI_pca   = "CHEI (PCA weights)"
)

# Keep your socio vars + labels from above
socio_vars <- c("income_ln","age65_79","age80_",
                "percen_eu","percen_noneu","socisss","av_lvng")

var_labels <- c(
  income_ln    = "Income (log)",
  age65_79     = "Age 65–79",
  age80_       = "Age 80+",
  percen_eu    = "EU foreigners",
  percen_noneu = "Non-EU foreigners",
  socisss      = "Social assistance",
  av_lvng      = "Living space per capita"
)

# -----------------------------
# 2) Make quartiles and long data
# -----------------------------
dat_q <- dat %>%
  mutate(across(all_of(socio_vars), ~ntile(., 4), .names = "{.col}_q"))

plot_df <- dat_q %>%
  select(muni_id, num_rsd, ends_with("_q"), all_of(indicator_vars)) %>%
  pivot_longer(ends_with("_q"), names_to = "socio_var", values_to = "quartile") %>%
  mutate(socio_var = sub("_q$", "", socio_var)) %>%
  pivot_longer(all_of(indicator_vars), names_to = "indicator", values_to = "heat_value")

# -----------------------------
# 3) Weighted summaries (raw pop weights)
# -----------------------------
weighted_se_mean <- function(x, w){
  ok <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w)
  if (!any(ok)) return(NA_real_)
  var_w <- Hmisc::wtd.var(x[ok], weights = w[ok], na.rm = TRUE)
  sw    <- sum(w[ok])
  sw2   <- sum(w[ok]^2)
  Neff  <- ifelse(sw2 > 0, (sw^2) / sw2, NA_real_)
  sqrt(var_w / Neff)
}

summary_df <- plot_df %>%
  group_by(socio_var, quartile, indicator) %>%
  summarise(
    mean_heat = Hmisc::wtd.mean(heat_value, weights = num_rsd, na.rm = TRUE),
    se_heat   = weighted_se_mean(heat_value, num_rsd),
    .groups   = "drop"
  )

# -----------------------------
# 4) Build y-axis order and shading
# -----------------------------
y_levels <- unlist(lapply(socio_vars, function(v) paste(var_labels[v], paste0("Q", 1:4))),
                   use.names = FALSE)

summary_df <- summary_df %>%
  mutate(
    socio_lab = var_labels[socio_var],
    indicator = factor(indicator, levels = indicator_vars, labels = indicator_labels),
    y_axis    = paste(socio_lab, paste0("Q", quartile)),
    y_axis    = factor(y_axis, levels = rev(y_levels))   # Q1 at top
  )

y_levels_rev <- rev(y_levels)
shade_df <- tibble(socio_var = socio_vars) %>%
  mutate(
    socio_lab  = var_labels[socio_var],
    block_rows = lapply(socio_lab, function(lbl) which(y_levels_rev %in% paste(lbl, paste0("Q", 1:4)))),
    ymin       = sapply(block_rows, min) - 0.5,
    ymax       = sapply(block_rows, max) + 0.5,
    fill_col   = rep(c("gray95","white"), length.out = length(socio_vars))
  )

# -----------------------------
# 5) Plot CHEI_equal vs CHEI_pca (weighted)
# -----------------------------
p_chei <- ggplot(summary_df, aes(x = mean_heat, y = y_axis)) +
  geom_rect(data = shade_df,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = fill_col),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_fill_identity() +
  geom_point(size = 1.5, color = "brown") +
  geom_errorbarh(
    aes(xmin = mean_heat - se_heat, xmax = mean_heat + se_heat),
    height = 0,
    color = "brown"
  ) +
  facet_wrap(~ indicator, ncol = 2, scales = "free_x") +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(hjust = 1),
    strip.text = element_text(face = "bold")
  )

print(p_chei)

# Save high-res
ggsave("CHEI_equal_vs_pca_weighted_rawpop.png", p_chei,
       width = 11, height = 7, dpi = 300)