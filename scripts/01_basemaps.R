library(tidyverse)
library(viridis)
library(marmap)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(patchwork)
library(ggnewscale)
library(ggrepel)

# Figura 1: Mapa de muestreo ---------------------------------------------------
# Coordenadas ---------------------------------------------------------------- #
# Generales
gen_xlim <- c(-18.1, 28.8)
gen_ylim <- c(-3.8, 46.8)

# España + Baleares
esp_xlim <- c(-10.5, 5.5)
esp_ylim <- c(35.0, 44.8)

# Canarias: archipiélago completo
can_xlim <- c(-18.8, -13.6)
can_ylim <- c(27.1, 29.9)

# Turquía / mar Egeo: centrado en Karaburun
tur_xlim <- c(23.2, 29.8)
tur_ylim <- c(37.1, 40.2)

# Ghana: centrado en Tema
gha_xlim <- c(-2.5, 2.5)
gha_ylim <- c(4.2, 7.2)

# Paleta de azules------------------------------------------------------------ #
blues <- c("#63727a",
           "#6a7a86",
           "#728291",
           "#7e90a2",
           "#94a6bc",
           "#a9bcd5",
           "#b9cae0",
           "#c6d3e3",
           "#d3dce6")

# Puntos de muestreo --------------------------------------------------------- #
coords <- read_csv("data/processed/coordinates.csv")

coords_gen <- coords %>%
  filter(lon >= gen_xlim[1], lon <= gen_xlim[2],
         lat >= gen_ylim[1], lat <= gen_ylim[2])

coords_esp <- coords %>%
  filter(lon >= esp_xlim[1], lon <= esp_xlim[2],
         lat >= esp_ylim[1], lat <= esp_ylim[2])

coords_can <- coords %>%
  filter(lon >= can_xlim[1], lon <= can_xlim[2],
         lat >= can_ylim[1], lat <= can_ylim[2])

coords_tur <- coords %>%
  filter(lon >= tur_xlim[1], lon <= tur_xlim[2],
         lat >= tur_ylim[1], lat <= tur_ylim[2])

coords_gha <- coords %>%
  filter(lon >= gha_xlim[1], lon <= gha_xlim[2],
         lat >= gha_ylim[1], lat <= gha_ylim[2])

# Escala de localidades común (evitar colores repetidos entre paneles) ------- #
all_localities <- sort(unique(coords$locality))

scale_locality <- scale_color_viridis_d(option = "plasma",
                                        limits = all_localities,
                                        drop   = FALSE,
                                        name   = "Locality")

# Batimetría ----------------------------------------------------------------- #
bathy_gen <- getNOAA.bathy(lon1 = gen_xlim[1],
                           lon2 = gen_xlim[2],
                           lat1 = gen_ylim[1],
                           lat2 = gen_ylim[2],
                           res = 2,
                           keep = TRUE,
                           path = "data/external/")

bathy_gen_xyz <- as.xyz(bathy_gen)

bathy_esp <- getNOAA.bathy(lon1 = esp_xlim[1],
                           lon2 = esp_xlim[2],
                           lat1 = esp_ylim[1],
                           lat2 = esp_ylim[2],
                           res = 1,
                           keep = TRUE,
                           path = "data/external/")

bathy_esp_xyz <- as.xyz(bathy_esp)

bathy_can <- getNOAA.bathy(lon1 = can_xlim[1],
                           lon2 = can_xlim[2],
                           lat1 = can_ylim[1],
                           lat2 = can_ylim[2],
                           res = 1,
                           keep = TRUE,
                           path = "data/external/")

bathy_can_xyz <- as.xyz(bathy_can)

bathy_tur <- getNOAA.bathy(lon1 = tur_xlim[1],
                           lon2 = tur_xlim[2],
                           lat1 = tur_ylim[1],
                           lat2 = tur_ylim[2],
                           res = 2,
                           keep = TRUE,
                           path = "data/external/")

bathy_tur_xyz <- as.xyz(bathy_tur)

bathy_gha <- getNOAA.bathy(lon1 = gha_xlim[1],
                           lon2 = gha_xlim[2],
                           lat1 = gha_ylim[1],
                           lat2 = gha_ylim[2],
                           res = 2,
                           keep = TRUE,
                           path = "data/external/")

bathy_gha_xyz <- as.xyz(bathy_gha)

# Keep only sea for fill
bathy_gen_xyz <- bathy_gen_xyz %>%
  filter(V3 < 0)

bathy_esp_xyz <- bathy_esp_xyz %>%
  filter(V3 < 0)

bathy_can_xyz <- bathy_can_xyz %>%
  filter(V3 < 0)

bathy_tur_xyz <- bathy_tur_xyz %>%
  filter(V3 < 0)

bathy_gha_xyz <- bathy_gha_xyz %>%
  filter(V3 < 0)

# Profundidad global (misma escala de azules en todos los paneles) ----------- #
depth_min <- min(c(bathy_gen_xyz$V3,
                   bathy_esp_xyz$V3,
                   bathy_can_xyz$V3,
                   bathy_tur_xyz$V3,
                   bathy_gha_xyz$V3),
                 na.rm = TRUE)

depth_breaks <- c(0,
                  -500,
                  -1000,
                  -2000,
                  -3000,
                  -4000,
                  -5000)

scale_depth_common <- scale_fill_gradientn(colours = blues,
                                           limits = c(depth_min, 0),
                                           breaks = depth_breaks,
                                           labels = as.character(depth_breaks),
                                           name = "Depth (m)",
                                           guide = guide_colorbar(title.position = "top",
                                                                  title.hjust = 0.5,
                                                                  direction = "vertical",
                                                                  reverse = FALSE,
                                                                  barwidth = unit(0.4, "cm"),
                                                                  barheight = unit(8, "cm")))

# Tierra -----------------------------------------------------------------------
world <- ne_countries(scale = "medium",
                      returnclass = "sf")

# COMMON THEME -----------------------------------------------------------------
theme_map <- theme_classic() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.ticks.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        panel.grid.major = element_line(color = "grey90",
                                        linewidth = 0.25),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1, 1, 1, 1))

# CAPA COMÚN DE PUNTOS + ETIQUETAS (con ggrepel) ------------------------------

point_labels_gen <- list(geom_point(data = coords_gen,
                                    aes(x = lon,
                                        y = lat,
                                        color = locality),
                                    size = 2.5,
                                    alpha = 0.95))
                         # geom_text_repel(data = coords_gen,
                         #                 aes(x = lon,
                         #                     y = lat,
                         #                     label = locality),
                         #                 size = 4,
                         #                 color = "black",
                         #                 fontface = "plain",
                         #                 max.overlaps = Inf,
                         #                 box.padding = 0.35,
                         #                 point.padding = 0.25,
                         #                 segment.color = "grey50",
                         #                 segment.size = 0.5,
                         #                 min.segment.length = 0))

point_labels_esp <- list(geom_point(data = coords_esp,
                                    aes(x = lon,
                                        y = lat,
                                        color = locality),
                                    size = 3,
                                    alpha = 0.95),
                         # geom_label_repel(data = coords_esp,
                         #                  aes(x = lon,
                         #                      y = lat,
                         #                      label = locality),
                         #                  size = 2.6,
                         #                  color = "black",
                         #                  fill = scales::alpha("white", 0.75),
                         #                  label.size = 0.15,
                         #                  show.legend = FALSE,
                         #                  max.overlaps = Inf,
                         #                  box.padding = 0.25,
                         #                  point.padding = 0.2,
                         #                  segment.color = "grey35",
                         #                  segment.size = 0.25,
                         #                  min.segment.length = 0))
                         geom_text_repel(data = coords_esp,
                                         aes(x = lon,
                                             y = lat,
                                             label = locality),
                                         size = 4,
                                         color = "black",
                                         fontface = "plain",
                                         max.overlaps = Inf,
                                         box.padding = 0.35,
                                         point.padding = 0.25,
                                         segment.color = "grey50",
                                         segment.size = 0.5,
                                         min.segment.length = 0))

point_labels_can <- list(geom_point(data = coords_can,
                                    aes(x = lon,
                                        y = lat,
                                        color = locality),
                                    size = 4,
                                    alpha = 0.95),
                         geom_text_repel(data = coords_can,
                                         aes(x = lon,
                                             y = lat,
                                             label = locality),
                                         size = 4,
                                         color = "black",
                                         fontface = "plain",
                                         max.overlaps = Inf,
                                         box.padding = 0.35,
                                         point.padding = 0.25,
                                         segment.color = "grey50",
                                         segment.size = 0.5,
                                         min.segment.length = 0))

point_labels_tur <- list(geom_point(data = coords_tur,
                                    aes(x = lon,
                                        y = lat,
                                        color = locality),
                                    size = 4,
                                    alpha = 0.95),
                         geom_text_repel(data = coords_tur,
                                         aes(x = lon,
                                             y = lat,
                                             label = locality),
                                         size = 4,
                                         color = "black",
                                         fontface = "plain",
                                         max.overlaps = Inf,
                                         box.padding = 0.35,
                                         point.padding = 0.25,
                                         segment.color = "grey50",
                                         segment.size = 0.5,
                                         min.segment.length = 0))

point_labels_gha <- list(geom_point(data = coords_gha,
                                    aes(x = lon,
                                        y = lat,
                                        color = locality),
                                    size = 4,
                                    alpha = 0.95),
                         geom_text_repel(data = coords_gha,
                                         aes(x = lon,
                                             y = lat,
                                             label = locality),
                                         size = 4,
                                         color = "black",
                                         fontface = "plain",
                                         max.overlaps = Inf,
                                         box.padding = 0.35,
                                         point.padding = 0.25,
                                         segment.color = "grey50",
                                         segment.size = 0.5,
                                         min.segment.length = 0))

# GENERAL ---------------------------------------------------------------
p_gen <- ggplot() +
  geom_tile(data = bathy_gen_xyz,
            aes(x = V1,
                y = V2,
                fill = V3)) +
  geom_contour(data = bathy_gen_xyz,
               aes(x = V1,
                   y = V2,
                   z = V3),
               breaks = c(-200,
                          -1000,
                          -3000),
               color = "#8f99a3",
               linewidth = 0.2) +
  geom_sf(data = world,
          fill = "lightgrey",
          color = "grey50",
          linewidth = 0.3) +
  point_labels_gen +
  annotation_north_arrow(location = "tl",
                         which_north = "true",
                         height = grid::unit(0.8, "cm"),
                         width  = grid::unit(0.8, "cm"),
                         style = north_arrow_orienteering(text_size = 7,
                                                          line_width = 0.6)) +
  annotation_scale(location = "bl",
                   width_hint = 0.25,
                   text_cex = 0.7) +
  coord_sf(xlim = gen_xlim,
           ylim = gen_ylim,
           expand = FALSE) +
  scale_fill_gradientn(colours = blues,
                       limits  = c(depth_min, 0),
                       breaks  = depth_breaks,
                       labels  = as.character(depth_breaks),
                       guide   = "none") +
  scale_locality +
  guides(color = "none") +
  labs(x = "Longitude", 
       y = "Latitude") +
  theme_map

p_gen

ggsave(plot = p_gen,
       filename = "panel_A_general.png",
       device = "png",
       path = "results/figures/maps/",
       units = "mm",
       width = 148,
       height = 80,
       dpi = 300)

# SPAIN + BALEARICS ------------------------------------------------------------
p_esp <- ggplot() +
  geom_tile(data = bathy_esp_xyz,
            aes(x = V1,
                y = V2,
                fill = V3)) +
  geom_contour(data = bathy_esp_xyz,
               aes(x = V1,
                   y = V2,
                   z = V3),
               breaks = c(-200,
                          -1000),
               color = "#8f99a3",
               linewidth = 0.2) +
  geom_sf(data = world,
          fill = "lightgrey",
          color = "grey50",
          linewidth = 0.3) +
  point_labels_esp  +
  ggspatial::annotation_scale(location = "bl",
                              width_hint = 0.25,
                              text_cex = 0.7) +
  coord_sf(xlim = esp_xlim,
           ylim = esp_ylim,
           expand = FALSE) +
  scale_depth_common +
  scale_locality +
  guides(color = "none") +
  labs(x = "Longitude",
       y = "Latitude") +
  theme_map

p_esp

ggsave(plot = p_esp,
       filename = "panel_B_esp_baleares.png",
       device = "png",
       path = "results/figures/maps/",
       units = "mm",
       width = 148,
       height = 95,
       dpi = 300)

# CANARY ISLANDS ---------------------------------------------------------------
p_can <- ggplot() +
  geom_tile(data = bathy_can_xyz,
            aes(x = V1,
                y = V2,
                fill = V3)) +
  geom_contour(data = bathy_can_xyz,
               aes(x = V1,
                   y = V2,
                   z = V3),
               breaks = c(-200,
                          -1000),
               color = "#8f99a3",
               linewidth = 0.2) +
  geom_sf(data = world,
          fill = "lightgrey",
          color = "grey50",
          linewidth = 0.3) +
  point_labels_can  +
  annotation_scale(location = "bl",
                   width_hint = 0.25,
                   text_cex = 0.7) +
  coord_sf(xlim = can_xlim,
           ylim = can_ylim,
           expand = FALSE) +
  scale_fill_gradientn(colours = blues,
                       limits  = c(depth_min, 0),
                       breaks  = depth_breaks,
                       labels  = as.character(depth_breaks),
                       guide   = "none") +
  scale_locality +
  guides(color = "none") +
  labs(x = "Longitude",
       y = "Latitude") +
  theme_map

p_can

ggsave(plot = p_can,
       filename = "panel_C_canarias.png",
       device = "png",
       path = "results/figures/maps/",
       units = "mm",
       width = 148,
       height = 80,
       dpi = 300)

# TURQUIA ----------------------------------------------------------------
p_tur <- ggplot() +
  geom_tile(data = bathy_tur_xyz,
            aes(x = V1,
                y = V2,
                fill = V3)) +
  geom_contour(data = bathy_tur_xyz,
               aes(x = V1,
                   y = V2,
                   z = V3),
               breaks = c(-200,
                          -1000,
                          -3000),
               color = "#8f99a3",
               linewidth = 0.2) +
  geom_sf(data = world,
          fill = "lightgrey",
          color = "grey50",
          linewidth = 0.3) +
  point_labels_tur  +
  annotation_scale(location = "bl",
                   width_hint = 0.25,
                   text_cex = 0.7) +
  coord_sf(xlim = tur_xlim,
           ylim = tur_ylim,
           expand = FALSE) +
  scale_fill_gradientn(colours = blues,
                       limits  = c(depth_min, 0),
                       breaks  = depth_breaks,
                       labels  = as.character(depth_breaks),
                       guide   = "none")  +
  scale_locality +
  guides(color = "none") +
  labs(x = "Longitude",
       y = "Latitude") +
  theme_map

p_tur

ggsave(plot = p_tur,
       filename = "panel_D_turquia.png",
       device = "png",
       path = "results/figures/maps/",
       units = "mm",
       width = 148,
       height = 80,
       dpi = 300)

# GULF OF GUINEA ---------------------------------------------------------------
p_gha <- ggplot() +
  geom_tile(data = bathy_gha_xyz,
            aes(x = V1,
                y = V2,
                fill = V3)) +
  geom_contour(data = bathy_gha_xyz,
               aes(x = V1,
                   y = V2,
                   z = V3),
               breaks = c(-200,
                          -1000,
                          -3000),
               color = "#8f99a3",
               linewidth = 0.2) +
  geom_sf(data = world,
          fill = "lightgrey",
          color = "grey50",
          linewidth = 0.3) +
  point_labels_gha  +
  annotation_scale(location = "bl",
                   width_hint = 0.25,
                   text_cex = 0.7) +
  coord_sf(xlim = gha_xlim,
           ylim = gha_ylim,
           expand = FALSE) +
  scale_fill_gradientn(colours = blues,
                       limits  = c(depth_min, 0),
                       breaks  = depth_breaks,
                       labels  = as.character(depth_breaks),
                       guide   = "none") +
  scale_locality +
  guides(color = "none") +
  labs(x = "Longitude", 
       y = "Latitude") +
  theme_map

p_gha

ggsave(plot = p_gha,
       filename = "panel_E_golfo_guinea.png",
       device = "png",
       path = "results/figures/maps/",
       units = "mm",
       width = 148,
       height = 80,
       dpi = 300)

design <- "
AABBBB
AABBBB
CCDDEE
"

final_map <- p_gen + p_esp + p_can + p_tur + p_gha +
  plot_layout(design  = design,
              widths  = rep(1, 6),
              heights = c(1, 1, 1.25),
              guides  = "collect") +
  plot_annotation(tag_levels = "A")

final_map

ggsave(plot = final_map,
       filename = "Figura1_mapa_muestreo_A4.png",
       device = "png",
       path = "results/figures/maps/",
       units = "mm",
       width = 297,
       height = 210,
       dpi = 300)

rm(list = ls())
gc()
