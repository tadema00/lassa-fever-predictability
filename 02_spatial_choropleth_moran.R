# ============================================================================
# 02_spatial_choropleth_moran.R
#
# National choropleth maps (Figures 1-3) + Global Moran's I spatial
# autocorrelation test (Section 2.9 / Section 3.1).
#
# Focal states are HIGHLIGHTED (outlined), never masked to zero, and the
# six-state cohort uses Plateau (not Benue).
#
# NOTE: this script needs a GADM Nigeria state-boundary shapefile at
# data/gadm41_NGA_1.shp (the .shp plus its .dbf/.shx/.prj siblings). That file
# was not among the three uploads (annual_data.xlsx, Cases_rainfal_data.xlsx,
# climate_features.csv), so this stage cannot run yet - the guard below fails
# fast with a clear message instead of a cryptic sf::st_read() error.
# ============================================================================

source("00_setup.R")

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
# Global Moran's I under Queen contiguity, per year, formally testing the
# Section 3.1 claim that LF's national footprint is not a single spatially
# contiguous cluster.
# ---------------------------------------------------------------------------
compute_morans_i_by_year <- function(gdf, annual_long) {
  nb <- spdep::poly2nb(gdf, queen = TRUE)
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  
  years <- sort(unique(annual_long$Year))
  results <- map_dfr(years, function(yr) {
    yr_data <- annual_long |>
      dplyr::filter(Year == yr)
    merged <- gdf %>%
      left_join(yr_data, by = c("NAME_1" = "State")) %>%
      mutate(Cases = replace_na(Cases, 0))
    test <- spdep::moran.test(merged$Cases, lw, zero.policy = TRUE)
    tibble(
      Year = yr,
      MoranI = unname(test$estimate["Moran I statistic"]),
      p_value = test$p.value
    )
  })
  results
}

morans_results <- compute_morans_i_by_year(gdf, annual_long)
print(morans_results)
write_csv(morans_results, file.path(TAB_DIR, "morans_i_by_year.csv"))

morans_plot <- ggplot(morans_results, aes(x = Year, y = MoranI)) +
  geom_col(aes(fill = p_value < 0.05)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#95A5A6"),
                    labels = c(`TRUE` = "p < 0.05", `FALSE` = "n.s."), name = NULL) +
  labs(title = "Global Moran's I, national Lassa fever case counts (Queen contiguity)",
       y = "Moran's I", x = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(FIG_DIR, "morans_i_by_year.png"), morans_plot, width = 7, height = 4.5, dpi = 200)

# ---------------------------------------------------------------------------
# Run the figure generation
# ---------------------------------------------------------------------------
plot_cluster_highlight(gdf)
plot_cumulative(gdf, annual_cumulative)
plot_annual_panels(gdf, annual_long)

cat("Spatial figures + Moran's I test complete.\n")

