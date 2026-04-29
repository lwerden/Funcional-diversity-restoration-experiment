#mappa costa ricaaa

# Load necessary libraries
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(ggspatial)
library(colorBlindness)
cvdPlot(c("#E69F00", "#56B4E9", "#009E73", "#D55E00"))

#FINAL VERSION!!!
# Get the map of Costa Rica
costa_rica <- ne_countries(scale = "medium", country = "Costa Rica", returnclass = "sf")

# Create a data frame for your sites
sites <- data.frame(
  name = c("FAB", "LLFS", "FCEA", "GAM"),  # Site names
  lon = c(-83.48893223605518, -82.9257529032521, -82.963622851001, -83.20150803146555),  # Longitudes
  lat = c(9.020756202074873, 8.740409777327136, 8.799968738174542, 8.701355351990747 )   # Latitudes
)

# Convert the sites data to an sf object
sites_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326)

# Plot the base map and add site points
Costa_rica_plot<- ggplot() +
  geom_sf(data = costa_rica, fill = "#F4F1E8",
          color = "grey20"
  ) +  # Base map
  geom_sf(data = sites_sf, aes(fill = name), color = "grey20", size = 3, shape = 21) +# Overlay sites
  scale_fill_manual(name = "Sites:",
                      values = c(
                        "FAB" = "olivedrab4", 
                        "FCEA" ="darkorange", 
                        "LLFS" = "firebrick4",
                        "GAM" = "royalblue1"
                       )) +
  annotation_scale(location = "bl", width_hint = 0.3) +  # Add scale bar
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering) +  # Compass
  theme_classic() +
  ggtitle("Map of Costa Rica with Site Locations") +
  theme_classic() +
  theme(legend.position = c(0.13, 0.28),
        legend.background = element_rect(fill = alpha("white", 0.8)),
        legend.title = element_text(face = "bold"))
  #theme(
  #  panel.background = element_rect(fill = "#EAF2F8", color = NA)
  #)
cvdPlot(plot = last_plot(), layout = c("olivedrab4", "goldenrod2", "firebrick4", "royalblue1"))
ggsave("~/Desktop/costarica4.jpeg", plot=Costa_rica_plot, width = 5, height = 5, dpi = 800)
ggsave("~/Desktop/costaricamap.svg", plot=Costa_rica_plot)