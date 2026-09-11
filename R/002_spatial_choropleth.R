# ============================================================================
# 02_spatial_choropleth_moran.R
# National choropleth maps (Figures 1-3) 

# NOTE: this script needs a GADM Nigeria state-boundary shapefile at
# data/gadm41_NGA_1.shp (the .shp plus its .dbf/.shx/.prj siblings). That file
# was not among the three uploads (annual_data.xlsx, Cases_rainfal_data.xlsx,
# climate_features.csv), so this stage cannot run yet - the guard below fails
# fast with a clear message instead of a cryptic sf::st_read() error.
# ============================================================================

source("000_setup.R")

if (!file.exists(SHAPEFILE)) {
  stop(
    "Shapefile not found at '", SHAPEFILE, "'. Download the Nigeria admin-1 ",
    "boundaries (GADM v4.1, NGA, level 1) and place gadm41_NGA_1.shp plus its ",
    ".dbf/.shx/.prj companions in the data/ folder before running this script."
  )
}

prepped           <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
annual_long        <- prepped$annual_long
annual_cumulative  <- prepped$annual_cumulative

gdf <- sf::st_read(SHAPEFILE, quiet = TRUE)

add_category <- function(df, case_col = "Cases") {
  df %>% mutate(Category = cut(.data[[case_col]], breaks = CASE_BREAKS,
                               labels = CASE_LABELS, right = TRUE))
}

# ---------------------------------------------------------------------------
# Figure 3: annual national choropleths, 2018-2025, focal states outlined
# ---------------------------------------------------------------------------
plot_annual_panels <- function(gdf, annual_long) {
  years <- sort(unique(annual_long$Year))
  
  panels <- lapply(years, function(yr) {
    yr_data <- annual_long |>
      dplyr::filter(Year == yr)
    merged <- gdf %>%
      left_join(yr_data, by = c("NAME_1" = "State")) %>%
      mutate(Cases = replace_na(Cases, 0)) %>%
      add_category()
    
    ggplot(merged) +
      geom_sf(aes(fill = Category), color = "#BDC3C7", linewidth = 0.15) +
      geom_sf(data = merged |>
                dplyr::filter(NAME_1 %in% FOCAL_STATES),
              fill = NA, color = "black", linewidth = 0.6) +
      scale_fill_manual(values = setNames(CASE_COLORS, CASE_LABELS),
                        drop = FALSE, na.value = "white") +
      labs(title = as.character(yr)) +
      theme_void() +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12,
                                      color = "#2C3E50"))
  })
  
  # shared legend built from a dummy plot (cowplot is loaded centrally in
  # 00_setup.R now, so cowplot::get_legend() resolves)
  legend_plot <- ggplot(gdf %>% mutate(Category = factor(CASE_LABELS[1], levels = CASE_LABELS))) +
    geom_sf(aes(fill = Category)) +
    scale_fill_manual(values = setNames(CASE_COLORS, CASE_LABELS), name = NULL, drop = FALSE) +
    theme_void() + theme(legend.position = "right")
  legend <- cowplot::get_legend(legend_plot)
  
  combined <- patchwork::wrap_plots(panels, ncol = 4)
  ggsave(file.path(FIG_DIR, "national_annual_choropleth.png"), combined,
         width = 19, height = 11, dpi = 220, bg = "white")
  cat("Saved", file.path(FIG_DIR, "national_annual_choropleth.png"), "\n")
}

# ---------------------------------------------------------------------------
# Figure 2: national cumulative choropleth, focal states labeled
# ---------------------------------------------------------------------------
plot_cumulative <- function(gdf, annual_cumulative) {
  merged <- gdf %>%
    left_join(annual_cumulative, by = c("NAME_1" = "State")) %>%
    mutate(Total = replace_na(Total, 0)) %>%
    add_category(case_col = "Total")
  
  centroids <- merged |>
    dplyr::filter(NAME_1 %in% FOCAL_STATES) %>%
    sf::st_centroid() %>%
    mutate(lon = sf::st_coordinates(.)[, 1], lat = sf::st_coordinates(.)[, 2])
  
  p <- ggplot(merged) +
    geom_sf(aes(fill = Category), color = "#BDC3C7", linewidth = 0.15) +
    geom_sf(data = merged %>% filter(NAME_1 %in% FOCAL_STATES),
            fill = NA, color = "black", linewidth = 0.6) +
    geom_text(data = centroids, aes(x = lon, y = lat, label = NAME_1),
              fontface = "bold", size = 2.8) +
    scale_fill_manual(values = setNames(CASE_COLORS, CASE_LABELS), name = NULL, drop = FALSE) +
    theme_void() +
    theme(legend.position = "right")
  
  ggsave(file.path(FIG_DIR, "national_cumulative_choropleth.png"), p,
         width = 8, height = 9.5, dpi = 200, bg = "white")
  cat("Saved", file.path(FIG_DIR, "national_cumulative_choropleth.png"), "\n")
}

# ---------------------------------------------------------------------------
# Figure 1: study area / ecological clusters (Southern Forest vs Northern
# Savannah), rest of country in neutral gray
# ---------------------------------------------------------------------------
plot_cluster_highlight <- function(gdf) {
  gdf2 <- gdf %>%
    mutate(Cluster = case_when(
      NAME_1 %in% SOUTHERN_CLUSTER ~ "Southern Forest Cluster",
      NAME_1 %in% NORTHERN_CLUSTER ~ "Northern Savannah Cluster",
      TRUE ~ "Other states (not in study)"
    ))
  
  colors <- c("Southern Forest Cluster" = "#2E7D32",
              "Northern Savannah Cluster" = "#C0392B",
              "Other states (not in study)" = "#EAECEE")
  
  centroids <- gdf2 |>
    dplyr::filter(NAME_1 %in% FOCAL_STATES) %>%
    sf::st_centroid() %>%
    mutate(lon = sf::st_coordinates(.)[, 1], lat = sf::st_coordinates(.)[, 2])
  
  p <- ggplot(gdf2) +
    geom_sf(aes(fill = Cluster), color = "#7F8C8D", linewidth = 0.2) +
    geom_text(data = centroids, aes(x = lon, y = lat, label = NAME_1),
              color = "white", fontface = "bold", size = 3) +
    scale_fill_manual(values = colors, name = NULL) +
    theme_void() +
    theme(legend.position = c(0.15, 0.1))
  
  ggsave(file.path(FIG_DIR, "study_area_ecological_clusters.png"), p,
         width = 8, height = 9, dpi = 200, bg = "white")
  cat("Saved", file.path(FIG_DIR, "study_area_ecological_clusters.png"), "\n")
}

# ---------------------------------------------------------------------------
# Run the figure generation
# ---------------------------------------------------------------------------
plot_cluster_highlight(gdf)
plot_cumulative(gdf, annual_cumulative)
plot_annual_panels(gdf, annual_long)

cat("Spatial figures \n")
