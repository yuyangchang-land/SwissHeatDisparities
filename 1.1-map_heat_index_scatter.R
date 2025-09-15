
library(dplyr)
library(scatterplot3d)

file_path <- "/Users/yuchang/Desktop/processed_living_area_with_code.csv"
df <- read.csv(file_path)

canton_summary <- df %>%
  group_by(canton) %>%
  summarise(
    mean_heatwave_prob = mean(heatwave_probability_mean, na.rm = TRUE),
    mean_heatwave_days = mean(heatwave_days_mean, na.rm = TRUE),
    mean_urban_LST = mean(Urban_LST_Mean, na.rm = TRUE)
  )

scatter_3d <- scatterplot3d(canton_summary$mean_heatwave_prob, 
                            canton_summary$mean_heatwave_days, 
                            canton_summary$mean_urban_LST,
                            pch = 19, color = "blue",
                            xlab = "Heatwave Probability",
                            ylab = "Heatwave Days",
                            zlab = "Urban LST Mean",
                            main = "3D Scatter Plot of Heat Stress Indicators by Canton")

text(scatter_3d$xyz.convert(canton_summary$mean_heatwave_prob, 
                            canton_summary$mean_heatwave_days, 
                            canton_summary$mean_urban_LST), 
     labels = canton_summary$canton, cex = 0.8, pos = 3, col = "red")

print(canton_summary)


# Load necessary libraries ######map x-y scatter fig with lst as color category
library(ggplot2)
library(dplyr)
library(ggrepel)  # For better label positioning

# Set file path (update this to your actual file location)
file_path <- "/Users/yuchang/Desktop/processed_living_area_with_code.csv"

# Load the dataset
df <- read.csv(file_path)

# Remove rows where the canton column is empty or NA
df_cleaned <- df %>%
  filter(!is.na(canton) & canton != "")

# Summarize data by canton
canton_summary <- df_cleaned %>%
  group_by(code, canton) %>%
  summarise(
    mean_heatwave_prob = mean(heatwave_probability_mean, na.rm = TRUE),
    mean_heatwave_days = mean(heatwave_days_mean, na.rm = TRUE),
    mean_urban_LST = mean(Urban_LST_Mean, na.rm = TRUE)
  )

# Customizable Titles (Modify These)
x_axis_title <- "Heatwave probability"  
y_axis_title <- "Heat warning days"  
colorbar_title <- "LST(°C)"  

# Scatter plot with borders and adjustable titles
ggplot(canton_summary, aes(x = mean_heatwave_prob, 
                           y = mean_heatwave_days, 
                           color = mean_urban_LST)) +
  geom_point(size = 3, alpha = 0.8) +  # Points
  scale_color_gradient(low = "blue", high = "red", name = colorbar_title) +  # Customizable colormap title
  labs(title = NULL,
       x = x_axis_title,  # Customizable x-axis title
       y = y_axis_title) +  # Customizable y-axis title
  geom_text_repel(aes(label = canton), size = 3, box.padding = 0.5, max.overlaps = 20) +  # Adjusted labels
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 12),  # Increase axis title size
    legend.title = element_text(size = 12),  # Adjust size of colorbar title
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add border
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.ticks = element_line(size = 0.5),  # Add tick marks
    axis.ticks.length = unit(5, "pt")  # Make tick marks slightly longer
  )

# Save plot as an image
ggsave("/Users/yuchang/Desktop/heatwave_scatter_plot_with_borders.png", width = 8, height = 6, dpi = 300)

# Print message
cat("Scatter plot saved as 'heatwave_scatter_plot_with_borders.png'\n")
