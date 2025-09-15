
# Load required libraries
library(sf)
library(ggplot2)
library(dplyr)
library(classInt)
library(tidyr)
library(cowplot)

# === Load Municipality Shapefile ===
shp_path <- "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/revision_r1/muni_with_zonal_means.shp"
df_sf <- st_read(shp_path)
colnames(df_sf)
df_sf <- st_read(shp_path) %>%
  select(id, lst_rbn, age65_, age80_, forgn19, for_eu,
         htwv_p_, avgt25_, socisss, av_lvng, geometry) %>%   # keep raw fields needed
  filter(lst_rbn != "X", age65_ != "X", htwv_p_ != "X", avgt25_ != "X") %>%
  mutate(
    # convert to numeric
    Urban_LST_ = as.numeric(lst_rbn),
    percen_65  = as.numeric(age65_),
    age80_     = as.numeric(age80_),
    forgn19    = as.numeric(forgn19),
    for_eu     = as.numeric(for_eu),
    socisss    = as.numeric(socisss),
    avgt25_    = as.numeric(avgt25_),
    htwv_p_    = as.numeric(htwv_p_),
    av_lvng    = as.numeric(av_lvng),
    
    # derive 65–79 group
    age65_79     = percen_65 - age80_,
    
    # foreigners split
    for_eu_prop  = ifelse(for_eu > 1, for_eu/100, for_eu),
    percen_eu    = forgn19 * for_eu_prop,
    percen_noneu = forgn19 * (1 - for_eu_prop),
    
    # rename heat indicators
    heatwave_p = htwv_p_,
    heatwave_d = avgt25_
  )
# === Load Swiss National Boundary ===
boundary_path <- "/Users/yuchang/PhD/PhD_research/Swiss_green_mortality/data_map/Boundary_admin/swissBOUNDARIES3D_1_5_TLM_KANTONSGEBIET.shp"
swiss_boundary <- st_read(boundary_path)
swiss_boundary <- st_zm(swiss_boundary, drop = TRUE, what = "ZM")
# === Define Color Schemes Based on Reference Images ===
color_palettes <- list(
  "Urban_LST_" = c("#e8e8e8", "#dfb0d6", "#be64ac",  
                   "#ace4e4", "#a5add3", "#8c62aa",
                   "#5ac8c8", "#5698b9", "#3b4994"),
  
  "heatwave_d" = c("#e8e8e8", "#E4ACAC", "#C85A5A",  
                   "#B0D5DF", "#AD9EA5", "#985356",  
                   "#64A0C8", "#626366", "#3B3B3B"), 
  
  "heatwave_p" = c("#e8e8e8", "#F0E0A0", "#E69F00", 
                   "#B5D3E7", "#B0B0B0", "#B15928",
                   "#4191C6", "#2C5173", "#000000")  # White → Light Blue → Medium Blue (Bottom Row)
)

# === Function to Create Bivariate Map ===
create_bivariate_map <- function(df, var_x, var_y, x_label, y_label, color_palette, title) {
  
  # Classify Data into Three Categories
  df <- df %>%
    mutate(
      x_class = cut(!!sym(var_x), breaks = quantile(!!sym(var_x), probs = seq(0, 1, length.out = 4), na.rm = TRUE), include.lowest = TRUE),
      y_class = cut(!!sym(var_y), breaks = quantile(!!sym(var_y), probs = seq(0, 1, length.out = 4), na.rm = TRUE), include.lowest = TRUE)
    )
  
  # Create Bivariate Color Scale
  bivariate_colors <- expand.grid(
    x_class = levels(df$x_class),
    y_class = levels(df$y_class)
  ) %>%
    mutate(fill_color = color_palette)
  
  # Merge color scale into data
  df <- df %>%
    left_join(bivariate_colors, by = c("x_class", "y_class"))
  
  # Plot the Map with Swiss Boundary
  map_plot <- ggplot() +
    geom_sf(data = df, aes(fill = fill_color), color = NA) +  
    geom_sf(data = swiss_boundary, fill = NA, color = "darkgray", size = 0.5) +  
    scale_fill_identity() +
    labs(title = title) +  # Add title for subfigure
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),  
      plot.title = element_text(size = 12, face = "bold"),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  # Create Tiny Square Legend
  legend_plot <- ggplot(bivariate_colors) +
    geom_tile(aes(x = as.numeric(x_class), y = as.numeric(y_class), fill = fill_color), width = 0.5, height = 0.5) +  
    scale_fill_identity() +
    
    # Move Y-axis label further left with smaller text
    annotate("text", x = -0.3, y = 2, label = y_label, size = 3, fontface = "bold", angle = 90) +
    
    # Move X-axis label further down with smaller text
    annotate("text", x = 2, y = -0.3, label = x_label, size = 3, fontface = "bold") +
    
    # Shorter arrows
    annotate("segment", x = 0.8, xend = 3, y = 0.2, yend = 0.2, arrow = arrow(length = unit(0.2, "cm")), size = 0.8) +
    annotate("segment", x = 0.2, xend = 0.2, y = 0.8, yend = 3, arrow = arrow(length = unit(0.2, "cm")), size = 0.8) +
    
    theme_void() +  # Removes grid, axis ticks, and background
    coord_fixed()
  
  # Combine the Map and Tiny Legend
  final_plot <- cowplot::plot_grid(
    map_plot,
    legend_plot, 
    ncol = 2, rel_widths = c(2, 0.3)# Reduce legend width
  )
  
  return(final_plot)
}

# === Generate Three Bivariate Maps without Titles ===
map_LST <- create_bivariate_map(df_sf, "Urban_LST_", "percen_noneu", "LST", "Non-EU foreigner", color_palettes[["Urban_LST_"]], title = NULL)
map_HeatwaveD <- create_bivariate_map(df_sf, "heatwave_d", "percen_noneu", "HWD", "Non-EU foreigner", color_palettes[["heatwave_d"]], title = NULL)
map_HeatwaveP <- create_bivariate_map(df_sf, "heatwave_p", "percen_noneu", "HWP", "Non-EU foreigner", color_palettes[["heatwave_p"]], title = NULL)

# === Combine Three Maps into One Big Figure (Closer Legends, No Labels) ===
final_combined_plot <- cowplot::plot_grid(
  map_LST, map_HeatwaveD, map_HeatwaveP, 
  ncol = 1, align = "hv", rel_heights = c(1, 1, 1),  # Ensure equal heights
  hjust = -0.5  # Moves the legends closer to the maps
)

# Display the final combined figure
print(final_combined_plot)





