# ============================================================================ #
# Genome-wide SNP discovery and a reduced diagnostic panel for geographic
# assignment of common octopus (Octopus vulgaris) in the Bay of Biscay and
# adjacent fishing regions
# 
# Data Analysis
# Author: Marina Parrondo Lombardía (parrondomarina@proton.me)
# ============================================================================ #

suppressPackageStartupMessages({
library(tidyverse)
library(ggokabeito)
library(vcfR)
library(snpAIMeR)
library(yaml)
library(adegenet)
library(hierfstat)
library(patchwork)
library(viridis)
})

# Load the VCF output from mPCRselect
# Cargo el vcf que sale de mPCRselect
vcf_fst_pi <- read.vcfR("data/processed/octopus.fst_pi.vcf.gz")
vcf_fst_pi

# Convert to Genind/ Convierte a genind
genind <- vcfR2genind(vcf_fst_pi)

# Assigns populations in the same order as the VCF
# Asigna poblaciones en el mismo orden que el VCF

# Importar y recuperar los nombres de columna originales
popmap <- readr::read_csv("data/processed/popmap.csv",
                          col_types = readr::cols(.default = readr::col_character())) |>
  dplyr::rename(id = Sample,
                pop = Population)

# Emparejar con los individuos del genind
idx <- match(indNames(genind), popmap$id)

stopifnot(!anyDuplicated(popmap$id),
          !anyNA(idx),
          !anyNA(popmap$pop[idx]),
          all(nzchar(popmap$pop[idx])))

# Asignar poblaciones y verificar
pop(genind) <- factor(popmap$pop[idx])
table(pop(genind), useNA = "ifany")

# Preselect top SNPs by Fst for snpAIMeR ---------------------------------------
# Preseleccionar los SNP top según el valor Fst para snpAIMeR ---------------- #

# Convert genind to hierfstat format
# Convertir genind a formato hierfstat 
hf_data <- hierfstat::genind2hierfstat(genind,
                                       pop = as.integer(droplevels(pop(genind))))

# Calculate locus-specific statistics (includes locus-specific Fst)
# Calcular estadísticos por locus (incluye Fst por locus)
bs <- hierfstat::basic.stats(hf_data, diploid = TRUE)

# Extract Fst by locus and sort
# Extraer Fst por locus y ordenar
fst_df <- data.frame(locus = locNames(genind),
                     Fst = bs$perloc$Fst) |>
  dplyr::filter(is.finite(Fst)) |>
  dplyr::arrange(dplyr::desc(Fst))

stopifnot(nrow(fst_df) >= 15L)

# Top 15 SNPs by Fst
# Top 15 SNPs por Fst
print(head(fst_df, 15))

# Subsett genind to top 15
# Subsetear genind a top 15
top15 <- fst_df$locus[seq_len(15L)]
genind_top15 <- genind[loc = top15]

stopifnot(adegenet::nLoc(genind_top15) == 15L,
          setequal(adegenet::locNames(genind_top15), top15))

# PCA of the top 15 SNPs
# PCA de los top 15 SNPs
x <- tab(genind_top15,
         NA.method = "mean")

pca_res <- dudi.pca(x,
                    scannf = FALSE,
                    nf = 3)

pca_df <- data.frame(PC1 = pca_res$li[,1],
                     PC2 = pca_res$li[,2],
                     pop = pop(genind_top15))

ggplot(pca_df,
       aes(PC1,
           PC2,
           color = pop)) +
  geom_point(size = 3,
             alpha = 0.7) +
  theme_classic() +
  labs(title = "PCA — Top 15 SNPs por Fst")

# Export to STRUCTURE / Exporta a STRUCTURE
genind2structure <- function(obj, file = "", pops = FALSE){
  if(!"genind" %in% class(obj)){
    warning("Function was designed for genind objects.")
  }
  pl <- max(obj@ploidy)
  S <- adegenet::nInd(obj)
  tab <- data.frame(ind = rep(adegenet::indNames(obj), each = pl))
  if(pops){
    popnums <- 1:adegenet::nPop(obj)
    names(popnums) <- as.character(unique(adegenet::pop(obj)))
    popcol <- rep(popnums[as.character(adegenet::pop(obj))], each = pl)
    tab <- cbind(tab, data.frame(pop = popcol))
  }
  loci <- adegenet::locNames(obj)
  tab <- cbind(tab, matrix(-9, nrow = dim(tab)[1], ncol = adegenet::nLoc(obj),
                           dimnames = list(NULL,loci)))
  for(L in loci){
    thesegen <- obj@tab[,grep(paste("^", L, "\\.", sep = ""),
                              dimnames(obj@tab)[[2]]), drop = FALSE]
    al <- 1:dim(thesegen)[2]
    for(s in 1:S){
      if(all(!is.na(thesegen[s,]))){
        tabrows <- (1:dim(tab)[1])[tab[[1]] == adegenet::indNames(obj)[s]]
        tabrows <- tabrows[1:sum(thesegen[s,])]
        tab[tabrows,L] <- rep(al, times=thesegen[s,])
      }
    }
  }
  write.table(tab,
              file = file,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
}

# Exportar a STRUCTURE
genind2structure(genind_top15,
                 file = "data/processed/octopus_top15.str",
                 pops = TRUE)

# Guardar la correspondencia entre códigos y poblaciones
# Utiliza el mismo orden que la función genind2structure()
population_order <- unique(as.character(pop(genind_top15)))

population_key <- data.frame(code = seq_along(population_order),
                             population = population_order)

readr::write_csv(population_key,
                 "data/processed/octopus_top15_population_key.csv")

# Eliminar el encabezado para el archivo de entrada de snpAIMeR
lines <- readLines("data/processed/octopus_top15.str")

writeLines(lines[-1],
           "data/processed/octopus_top15_fixed.str")

rm(lines)

# Quick check / Verificación rápida
str_top15 <- read.table("data/processed/octopus_top15_fixed.str",
                        header = FALSE,
                        sep = "\t")
stopifnot(nrow(str_top15) == 2L * adegenet::nInd(genind_top15),
          ncol(str_top15) == adegenet::nLoc(genind_top15) + 2L,
          !anyNA(str_top15))

rm(str_top15)

# Configuration / Configuración
config <- list(min_range = 1L,
               max_range = 15L,
               assignment_rate_threshold = 0.9,
               cross_validation_replicates = 1000L,
               working_directory = normalizePath("results/snpAIMeR/"),
               structure_file = normalizePath("data/processed/octopus_top15_fixed.str"),
               number_of_individuals = 452L,
               number_of_loci = 15L,
               one_data_row_per_individual = FALSE,
               column_sample_IDs = 1L,
               column_population_assignments = 2L,
               column_other_info = NULL,
               row_markernames = 0L, # sin fila de nombres de marcadores
               no_genotype_character = -9L,
               optional_population_info = NULL,
               genotype_character_separator = NULL)

write_yaml(config,
           "data/processed/snpAIMeR_config.yaml")

# Run / Ejecutamos
# snpAIMeR("non-interactive",
#          "data/processed/snpAIMeR_config.yaml",
#          verbose = TRUE)

all_comb <- read.csv("results/snpAIMeR/All_combinations_assign_rate.csv")
above <- read.csv("results/snpAIMeR/Above_threshold_assign_rate.csv")
panel <- read.csv("results/snpAIMeR/Panel_size_assign_rate.csv")

# Combinations evaluated / Combinaciones evaluadas
nrow(all_comb)
# Combinations above the 0.9 threshold
# Combinaciones sobre umbral 0.9
nrow(above)
# Panel size / Tamaño de panel
print(panel)

head(all_comb)
colnames(all_comb)

# Top 10 best global combinations
# Top 10 mejores combinaciones globales
head(all_comb[order(-all_comb$avg_success_rate), ], 10)

# Add a column with the number of loci
# Añadir columna con número de loci
all_comb$panel_size <- sapply(strsplit(all_comb$marker, ", "), length)

# Best combination for each panel size
# Mejor combinación por cada tamaño de panel
best_per_size <- do.call(rbind, lapply(split(all_comb, all_comb$panel_size), 
                                       function(x) x[which.max(x$avg_success_rate), ]))
print(best_per_size)

# Table S3. Best-performing combination at each panel size ---------------------
# Añadir identificadores originales a los códigos L01–L15

output_dir <- "results/tables/supplementary_tables"

dir.create(output_dir,
           recursive = TRUE,
           showWarnings = FALSE)

# Recuperar el orden de los loci del archivo STRUCTURE original
# Usar el archivo CON encabezado, no octopus_top15_fixed.str.
structure_file <- "data/processed/octopus_top15.str"

header <- readLines(structure_file,
                    n = 1)

stopifnot(length(header) == 1L)

header_fields <- strsplit(header, "\t", fixed = TRUE)[[1]]

# Las dos primeras columnas son ind y pop.
stopifnot(length(header_fields) > 2L,
          identical(header_fields[1:2], c("ind", "pop")))

original_loci <- header_fields[-(1:2)]

stopifnot(length(original_loci) == 15L,
          !anyNA(original_loci),
          all(nzchar(original_loci)),
          !anyDuplicated(original_loci))

# Crear la correspondencia L01–L15 → identificador original 
marker_key <- data.frame(marker_code = sprintf("L%02d",
                                               seq_along(original_loci)),
                         locus = original_loci,
                         stringsAsFactors = FALSE)

# Añadir Fst y posición en el ranking de diferenciación.
stopifnot(all(c("locus", "Fst") %in% names(fst_df)),
          !anyDuplicated(fst_df$locus),
          all(original_loci %in% fst_df$locus))

fst_ranking <- fst_df |>
  dplyr::filter(is.finite(Fst)) |>
  dplyr::arrange(dplyr::desc(Fst)) |>
  dplyr::mutate(Fst_rank = dplyr::row_number()) |>
  dplyr::select(locus, Fst, Fst_rank)

marker_key <- marker_key |>
  dplyr::left_join(fst_ranking,
                   by = "locus")

stopifnot(!anyNA(marker_key$Fst),
          !anyNA(marker_key$Fst_rank),
          setequal(marker_key$Fst_rank, seq_len(15L)))

# Mapa con los nombres originales: conserva ":".
locus_map <- setNames(marker_key$locus,
                      marker_key$marker_code)

print(marker_key, row.names = FALSE)

# Traducir los códigos de cada combinación de best_per_size ---------------

best <- as.data.frame(best_per_size)

stopifnot(all(c("marker", "avg_success_rate", "panel_size") %in% names(best)),
  !anyNA(best$marker),
  all(is.finite(best$avg_success_rate)))

marker_lists <- lapply(strsplit(as.character(best$marker), ",", fixed = TRUE),
                       trimws)

stopifnot(all(vapply(marker_lists,
                     function(z) {
                       all(nzchar(z)) &&
                         !anyDuplicated(z) &&
                         all(z %in% names(locus_map))},
                     logical(1))))

best$panel_size <- as.integer(as.character(best$panel_size))

stopifnot(!anyNA(best$panel_size),
          all(lengths(marker_lists) == best$panel_size),
          !anyDuplicated(best$panel_size),
          setequal(best$panel_size, seq_len(15L)))

# Añadir los loci originales en el mismo orden que los códigos.
best$original_loci <- vapply(marker_lists,
                             function(codes) {
                               paste(unname(locus_map[codes]),
                                     collapse = "; ")},
                             character(1))

# Normalizar la presentación de los códigos.
best$marker_codes <- vapply(marker_lists,
                            paste,
                            character(1),
                            collapse = ", ")

# Tabla S3
table_s3 <- best |>
  dplyr::arrange(panel_size) |>
  dplyr::transmute(`Panel size (SNPs)` = panel_size,
                   `Marker codes` = marker_codes,
                   `Original locus identifiers` = original_loci,
                   `Mean correct assignment rate` = avg_success_rate)

print(tibble::as_tibble(table_s3),
      n = Inf,
      width = Inf)

# Exportar
readr::write_csv(table_s3,
  file.path(output_dir,
    "Table_S3_best_combinations_with_locus_names.csv"))

readr::write_csv(marker_key,
                 file.path(output_dir,
                           "Marker_code_locus_Fst_correspondence.csv"))


# Seleccionar la mejor combinación de diez SNPs
best_10 <- all_comb |>
  dplyr::filter(panel_size == 10L,
                is.finite(avg_success_rate)) |>
  dplyr::slice_max(avg_success_rate,
                   n = 1,
                   with_ties = FALSE)

stopifnot(nrow(best_10) == 1L)

# Recuperar los códigos de marcadores de la combinación seleccionada
panel_10 <- trimws(strsplit(best_10$marker[[1]], ",", fixed = TRUE)[[1]])

# Leer el encabezado original y excluir las columnas ind y pop
header <- readLines("data/processed/octopus_top15.str",
                    n = 1)

original_loci <- strsplit(header,
                          "\t",
                          fixed = TRUE)[[1]][-(1:2)]

# Asociar los códigos de snpAIMeR con los nombres originales
locus_map <- setNames(original_loci,
                      sprintf("L%02d", seq_along(original_loci)))

stopifnot(length(original_loci) == adegenet::nLoc(genind_top15),
          !anyDuplicated(original_loci),
          setequal(original_loci, adegenet::locNames(genind_top15)),
          length(panel_10) == 10L,
          !anyDuplicated(panel_10),
          all(panel_10 %in% names(locus_map)))

# Recuperar los loci originales y crear el subconjunto
panel_loci <- unname(locus_map[panel_10])

genind_panel10 <- genind_top15[, loc = panel_loci]

stopifnot(adegenet::nLoc(genind_panel10) == 10L,
          adegenet::nInd(genind_panel10) == adegenet::nInd(genind_top15),
          setequal(adegenet::locNames(genind_panel10), panel_loci))

# Resumen del panel seleccionado
panel10_summary <- data.frame(marker = panel_10,
                              locus = panel_loci)

print(panel10_summary)
print(best_10$avg_success_rate)
genind_panel10

# Orden explícito de las poblaciones
pop_order <- c("northern_iberian_atlantic",
               "southern_iberian_atlantic",
               "western_mediterranean",
               "macaronesia")

# Etiquetas vinculadas a los nombres originales
pop_labels <- c(northern_iberian_atlantic = "Northern Iberian Atlantic",
                southern_iberian_atlantic = "Southern Iberian Atlantic",
                western_mediterranean = "Western Mediterranean",
                macaronesia = "Macaronesia")

stopifnot(!anyNA(pop(genind_top15)),
          setequal(unique(as.character(pop(genind_top15))),
                   pop_order))

# Seleccionar la mejor combinación para un tamaño de panel
get_best_panel <- function(panel_size, all_comb, locus_map) {
  
  best <- all_comb |>
    dplyr::filter(.data$panel_size == .env$panel_size,
                  is.finite(avg_success_rate)) |>
    dplyr::slice_max(avg_success_rate,
                     n = 1,
                     with_ties = FALSE)
  
  stopifnot(nrow(best) == 1L)
  
  markers <- trimws(strsplit(best$marker[[1]], ",", fixed = TRUE)[[1]])
  
  stopifnot(length(markers) == panel_size,
            !anyDuplicated(markers),
            all(markers %in% names(locus_map)))
  
  loci <- unname(locus_map[markers])
  
  stopifnot(!anyNA(loci),
            !anyDuplicated(loci))
  
  list(size = panel_size,
       markers = markers,
       loci = loci,
       assignment_rate = best$avg_success_rate[[1]])
}

# PCA y gráfico de un panel seleccionado
run_pca_panel <- function(genind_obj, selected_panel) {
  
  stopifnot(all(selected_panel$loci %in% adegenet::locNames(genind_obj)))
  
  genind_panel <- genind_obj[, loc = selected_panel$loci]
  
  stopifnot(adegenet::nLoc(genind_panel) == selected_panel$size,
            adegenet::nInd(genind_panel) == adegenet::nInd(genind_obj))
  
  x <- adegenet::tab(genind_panel,
                     NA.method = "mean")
  
  pca <- ade4::dudi.pca(x,
                        center = TRUE,
                        scale = TRUE,
                        scannf = FALSE,
                        nf = 2)
  
  pca_df <- data.frame(sample_id = adegenet::indNames(genind_panel),
                       PC1 = pca$li[, 1],
                       PC2 = pca$li[, 2],
                       pop = factor(as.character(adegenet::pop(genind_panel)),
                                    levels = pop_order))
  
  var_exp <- 100 * pca$eig / sum(pca$eig)
  
  p <- ggplot(pca_df,
              aes(PC1, PC2, color = pop)) +
    geom_point(alpha = 0.7,
               size = 4) +
    ggokabeito::scale_color_okabe_ito(order = c(5, 1, 3, 7),
                                      limits = pop_order,
                                      breaks = pop_order,
                                      labels = unname(pop_labels[pop_order]),
                                      drop = FALSE) +
    labs(title = paste(selected_panel$size, "SNPs"),
         subtitle = sprintf("Assignment rate: %.2f%%",
                            100 * selected_panel$assignment_rate),
         x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2]),
         color = "Population") +
    theme_classic()
  
  ggsave(filename = paste0("results/snpAIMeR/PCA_panel",
                           selected_panel$size,
                           ".png"),
         plot = p,
         width = 9,
         height = 6,
         dpi = 300)
  
  list(genind = genind_panel,
       pca = pca,
       plot = p,
       variance = var_exp,
       selection = selected_panel)
}

# Tamaños de panel que se compararán
panel_sizes <- c(3L, 5L, 10L)

# Recuperar las mejores combinaciones actuales
selected_panels <- lapply(panel_sizes,
                          get_best_panel,
                          all_comb = all_comb,
                          locus_map = locus_map)

names(selected_panels) <- as.character(panel_sizes)

# Calcular las PCA
pca_results <- lapply(selected_panels,
  function(selected_panel) {
    run_pca_panel(genind_obj = genind_top15,
                  selected_panel = selected_panel)
  }
)

# Tabla de correspondencias de los paneles comparados
panels_summary <- dplyr::bind_rows(
  lapply(selected_panels, function(z) {
    data.frame(panel_size = z$size,
               marker = z$markers,
               locus = z$loci,
               assignment_rate = z$assignment_rate)
  })
)

print(panels_summary)

readr::write_csv(panels_summary,
                 "results/snpAIMeR/PCA_panels_marker_summary.csv")

# Figura conjunta
fig_pca <- (patchwork::wrap_plots(
  lapply(pca_results, function(z) z$plot),
  nrow = 1) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A")) &
  theme(legend.position = "bottom",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        plot.tag = element_text(size = 14))

fig_pca

ggsave("results/snpAIMeR/PCA_panels_comparison.png",
       plot = fig_pca,
       width = 20,
       height = 7,
       units = "in",
       dpi = 600)

# ---------------------------------------------------------------------------- #
# Evaluación complementaria de paneles fijos: 3, 5 y 10 SNPs
# No repite la selección de marcadores dentro de entrenamiento.
# ---------------------------------------------------------------------------- #

out_dir <- "results/snpAIMeR/complementary_CV"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cv_replicates <- 1000L
pairwise_replicates <- 500L
train_fraction <- 0.75
base_seed <- 2026L

unit_order <- c("northern_iberian_atlantic",
                "southern_iberian_atlantic",
                "western_mediterranean",
                "macaronesia")

unit_labels <- c(northern_iberian_atlantic = "Northern Iberian Atlantic",
                 southern_iberian_atlantic = "Southern Iberian Atlantic",
                 western_mediterranean = "Western Mediterranean",
                 macaronesia = "Macaronesia")

stopifnot(exists("genind_top15"),
          exists("all_comb"),
          exists("locus_map"),
          adegenet::nInd(genind_top15) == 452L,
          !anyNA(adegenet::pop(genind_top15)),
          setequal(unique(as.character(adegenet::pop(genind_top15))),
                   unit_order),
          all(genind_top15@ploidy == 2L),
          !anyDuplicated(adegenet::indNames(genind_top15)))

sample_sizes <- data.frame(unit = unit_order,
                           n = as.integer(table(factor(
                             as.character(adegenet::pop(genind_top15)),
                             levels = unit_order))))

stopifnot(all(sample_sizes$n >= 4L),
          sum(sample_sizes$n) == 452L)

readr::write_csv(sample_sizes,
                 file.path(out_dir, "sample_sizes.csv"))

print(sample_sizes)

# Combinaciones: comprobaciones y dispersión por tamaño
stopifnot(all(c("marker", "avg_success_rate") %in% names(all_comb)),
          !anyNA(all_comb$marker),
          all(nzchar(all_comb$marker)),
          all(is.finite(all_comb$avg_success_rate)),
          all(all_comb$avg_success_rate >= 0),
          all(all_comb$avg_success_rate <= 1))

marker_sets <- strsplit(as.character(all_comb$marker),
                        ",",
                        fixed = TRUE)

marker_sets <- lapply(marker_sets, trimws)

stopifnot(all(vapply(marker_sets,
                     function(z) {
                       all(nzchar(z)) &&
                         !anyDuplicated(z) &&
                         all(z %in% names(locus_map))},
                     logical(1))))

all_comb$panel_size <- lengths(marker_sets)

combination_keys <- vapply(marker_sets,
                           function(z) paste(sort(z), collapse = "|"),
                           character(1))

stopifnot(!anyDuplicated(combination_keys))

expected_loci <- adegenet::nLoc(genind_top15)

stopifnot(expected_loci == 15L,
          nrow(all_comb) == 2^expected_loci - 1)

combination_counts <- all_comb |>
  dplyr::count(panel_size, name = "observed") |>
  dplyr::mutate(expected = choose(expected_loci, panel_size))

stopifnot(setequal(combination_counts$panel_size, seq_len(expected_loci)),
          all(combination_counts$observed == combination_counts$expected))

panel_dispersion <- all_comb |>
  dplyr::group_by(panel_size) |>
  dplyr::summarise(n_comb = dplyr::n(),
                   mean = mean(avg_success_rate),
                   sd = if (dplyr::n() > 1L) sd(avg_success_rate) else NA_real_,
                   q025 = quantile(avg_success_rate, 0.025, names = FALSE),
                   q975 = quantile(avg_success_rate, 0.975, names = FALSE),
                   min = min(avg_success_rate),
                   max = max(avg_success_rate),
                   .groups = "drop")

readr::write_csv(panel_dispersion,
                 file.path(out_dir, "panel_size_dispersion.csv"))

print(panel_dispersion)

# Comprobar la tabla de medias de snpAIMeR
if (exists("panel") &&
    is.data.frame(panel) &&
    all(c("panel_size", "avg_success_rate") %in% names(panel))) {
  
  stopifnot(!anyDuplicated(panel$panel_size))
  
  mean_check <- panel |>
    dplyr::select(panel_size,
                  snpAIMeR_mean = avg_success_rate) |>
    dplyr::full_join(panel_dispersion |>
                       dplyr::select(panel_size,
                                     recalculated_mean = mean),
                     by = "panel_size") |>
    dplyr::arrange(panel_size) |>
    dplyr::mutate(difference = snpAIMeR_mean - recalculated_mean,
                  matches = !is.na(difference) & abs(difference) < 1e-6)
  
  print(mean_check)
  
  readr::write_csv(mean_check,
                   file.path(out_dir, "panel_size_mean_check.csv"))
}

best_by_size <- all_comb |>
  dplyr::group_by(panel_size) |>
  dplyr::slice_max(avg_success_rate,
                   n = 1,
                   with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::arrange(panel_size)

readr::write_csv(best_by_size,
                 file.path(out_dir, "best_combination_per_size.csv"))

p_dispersion <- ggplot(panel_dispersion,
                       aes(panel_size, mean)) +
  geom_ribbon(aes(ymin = q025,
                  ymax = q975),
              alpha = 0.18) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_point(data = best_by_size,
             aes(x = panel_size,
                 y = avg_success_rate),
             inherit.aes = FALSE,
             shape = 17,
             size = 2.4) +
  geom_hline(yintercept = 0.90,
             linetype = "dashed",
             linewidth = 0.4) +
  scale_x_continuous(breaks = seq_len(expected_loci)) +
  labs(x = "Panel size (number of loci)",
       y = "Correct assignment rate",
       caption = paste("Circles: mean across marker combinations.",
                       "Band: empirical 2.5th–97.5th percentiles across combinations.",
                       "Triangles: best combination. Dashed line: 0.90 threshold.")) +
  theme_classic()

ggsave(file.path(out_dir, "panel_size_with_dispersion.png"),
       plot = p_dispersion,
       width = 8,
       height = 5,
       dpi = 300)

# Recuperar los paneles actuales sin listas manuales de loci
panel_sizes <- c(3L, 5L, 10L)

selected_panels <- setNames(
  lapply(panel_sizes, function(s) {
    
    best <- best_by_size |>
      dplyr::filter(panel_size == s)
    
    stopifnot(nrow(best) == 1L)
    
    markers <- trimws(strsplit(best$marker[[1]], ",",
                               fixed = TRUE)[[1]])
    
    stopifnot(length(markers) == s,
              !anyDuplicated(markers),
              all(markers %in% names(locus_map)))
    
    loci <- unname(locus_map[markers])
    
    stopifnot(!anyNA(loci),
              !anyDuplicated(loci),
              all(loci %in% adegenet::locNames(genind_top15)))
    
    list(size = s,
         markers = markers,
         loci = loci,
         snpAIMeR_rate = best$avg_success_rate[[1]])
  }),
  as.character(panel_sizes))

selected_panel_table <- dplyr::bind_rows(
  lapply(selected_panels, function(z) {
    data.frame(panel_size = z$size,
               marker = z$markers,
               locus = z$loci,
               snpAIMeR_rate = z$snpAIMeR_rate)
  }))

readr::write_csv(
  selected_panel_table,
  file.path(out_dir, "selected_panels.csv")
)

print(selected_panel_table)


# Una única función de validación DAPC
cv_dapc_fixed <- function(
    gi,
    loci,
    units,
    n_rep = 1000L,
    train_frac = 0.75,
    balance = FALSE,
    n_bal = NULL,
    prediction_priors = c("proportional", "uniform"),
    max_pcs = 20L,
    seed = 1L) {
  
  stopifnot(n_rep >= 1L,
            train_frac > 0,
            train_frac < 1,
            max_pcs >= 1L,
            !anyDuplicated(loci),
            all(loci %in% adegenet::locNames(gi)),
            !anyNA(adegenet::pop(gi)),
            all(prediction_priors %in% c("proportional", "uniform")))
  
  prediction_priors <- unique(prediction_priors)
  
  keep <- as.character(adegenet::pop(gi)) %in% units
  gi <- gi[keep, loc = loci]
  
  stopifnot(adegenet::nLoc(gi) == length(loci),
            !anyDuplicated(adegenet::indNames(gi)))
  
  groups <- factor(as.character(adegenet::pop(gi)),
                   levels = units)
  
  sizes <- table(groups)
  
  stopifnot(length(units) >= 2L,
            all(sizes >= 4L))
  
  if (balance) {
    if (is.null(n_bal)) {
      n_bal <- min(sizes)
    }
    
    stopifnot(length(n_bal) == 1L,
              is.finite(n_bal),
              n_bal == as.integer(n_bal),
              n_bal >= 4L,
              n_bal <= min(sizes))
    
    n_bal <- as.integer(n_bal)
  }
  
  # Sin imputación previa: las medias se estiman dentro de cada entrenamiento.
  X <- adegenet::tab(gi,
                     freq = TRUE,
                     NA.method = "asis")
  
  storage.mode(X) <- "double"
  
  stopifnot(nrow(X) == adegenet::nInd(gi),
            !any(is.infinite(X)))
  
  ids <- adegenet::indNames(gi)
  set.seed(seed)
  
  predictions <- list()
  diagnostic_rows <- list()
  
  for (r in seq_len(n_rep)) {
    
    warning_messages <- character()
    fit_info <- list(n_train = NA_integer_,
                     n_test = NA_integer_,
                     n_variables = NA_integer_,
                     n_pca_requested = NA_integer_,
                     n_pca_retained = NA_integer_,
                     n_da_retained = NA_integer_)
    
    result <- tryCatch(
      withCallingHandlers({
        
        idx_use <- if (balance) {
          unlist(lapply(units, function(u) {
              sample(which(groups == u), size = n_bal)
            }),
            use.names = FALSE
          )
        } else {
          seq_along(groups)
        }
        
        idx_train <- unlist(
          lapply(units, function(u) {
            
            candidates <- idx_use[groups[idx_use] == u]
            n_train <- floor(train_frac * length(candidates))
            
            if (
              n_train < 2L ||
              n_train >= length(candidates)
            ) {
              stop("Partición entrenamiento/test no válida.")
            }
            
            sample(candidates, size = n_train)
          }),
          use.names = FALSE
        )
        
        idx_test <- setdiff(idx_use, idx_train)
        
        fit_info$n_train <- length(idx_train)
        fit_info$n_test <- length(idx_test)
        
        x_train <- X[idx_train, , drop = FALSE]
        x_test <- X[idx_test, , drop = FALSE]
        
        # Descartar variables completamente missing en entrenamiento.
        usable <- colSums(!is.na(x_train)) > 0L
        
        x_train <- x_train[, usable, drop = FALSE]
        x_test <- x_test[, usable, drop = FALSE]
        
        if (!ncol(x_train)) {
          stop("No hay variables observadas en entrenamiento.")
        }
        
        training_means <- colMeans(x_train, na.rm = TRUE)
        
        for (j in seq_len(ncol(x_train))) {
          x_train[is.na(x_train[, j]), j] <- training_means[j]
          x_test[is.na(x_test[, j]), j] <- training_means[j]
        }
        
        # Descartar variables invariantes usando solo entrenamiento.
        training_variance <- apply(x_train, 2, var)
        
        variable <- is.finite(training_variance) &
          training_variance > 0
        
        x_train <- x_train[, variable, drop = FALSE]
        x_test <- x_test[, variable, drop = FALSE]
        
        fit_info$n_variables <- ncol(x_train)
        
        if (!ncol(x_train)) {
          stop("No quedan variables polimórficas en entrenamiento.")
        }
        
        stopifnot(all(is.finite(x_train)),
                  all(is.finite(x_test)),
                  identical(colnames(x_train),
                            colnames(x_test)))
        
        n_pca_requested <- min(max_pcs,
                               ncol(X) - 1L,
                               ncol(x_train),
                               nrow(x_train) - 1L)
        
        if (n_pca_requested < 1L) {
          stop("No hay componentes principales disponibles.")
        }
        
        fit_info$n_pca_requested <- n_pca_requested
        
        training_groups <- factor(groups[idx_train],
                                  levels = units)
        
        fit <- adegenet::dapc(x_train,
                              grp = training_groups,
                              n.pca = n_pca_requested,
                              n.da = length(units) - 1L,
                              center = TRUE,
                              scale = FALSE,
                              pca.info = TRUE)
        
        fit_info$n_pca_retained <- fit$n.pca
        fit_info$n_da_retained <- fit$n.da
        
        fit_groups <- levels(fit$grp)
        
        stopifnot(setequal(fit_groups, units))
        
        prediction_results <- list()
        prediction_log <- list()
        
        for (prior_name in prediction_priors) {
          
          prior_prob <- if (prior_name == "uniform") {
            rep(1 / length(fit_groups), length(fit_groups))
          } else {
            as.numeric(
              table(training_groups)[fit_groups]
            ) / length(training_groups)
          }
          
          pr <- tryCatch({
            
            prediction <- stats::predict(
              fit,
              newdata = x_test,
              prior = prior_prob
            )
            
            assigned <- as.character(prediction$assign)
            posterior <- prediction$posterior
            
            stopifnot(length(assigned) == length(idx_test),
                      !anyNA(assigned),
                      all(assigned %in% units),
                      nrow(posterior) == length(idx_test),
                      all(is.finite(posterior)),
                      all(abs(rowSums(posterior) - 1) < 1e-6))
            
            data.frame(rep = r,
                       prior = prior_name,
                       ind = ids[idx_test],
                       true = as.character(groups[idx_test]),
                       assigned = assigned,
                       post_max = apply(posterior, 1, max),
                       correct = assigned == as.character(groups[idx_test]),
                       stringsAsFactors = FALSE)
            
          }, error = function(e) e)
          
          if (inherits(pr, "error")) {
            prediction_log[[prior_name]] <- data.frame(
              prior = prior_name,
              status = "prediction_error",
              message = conditionMessage(pr)
            )
          } else {
            prediction_results[[prior_name]] <- pr
            
            prediction_log[[prior_name]] <- data.frame(
              prior = prior_name,
              status = "ok",
              message = ""
            )
          }
        }
        
        list(predictions = dplyr::bind_rows(prediction_results),
             log = dplyr::bind_rows(prediction_log))
        
      }, warning = function(w) {
        warning_messages <<- c(
          warning_messages,
          conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }),
      error = function(e) e
    )
    
    if (inherits(result, "error")) {
      
      diagnostic_rows[[r]] <- data.frame(rep = r,
                                         prior = prediction_priors,
                                         status = "fit_or_preprocessing_error",
                                         message = conditionMessage(result),
                                         warnings = paste(unique(warning_messages),
                                                          collapse = " | "),
                                         as.data.frame(fit_info),
                                         stringsAsFactors = FALSE)
      
    } else {
      
      predictions[[r]] <- result$predictions
      
      diagnostic_rows[[r]] <- cbind(
        data.frame(rep = r),
        result$log,
        data.frame(warnings = paste(unique(warning_messages),
                                    collapse = " | ")),
        as.data.frame(fit_info))
    }
  }
  
  list(predictions = dplyr::bind_rows(predictions),
       diagnostics = dplyr::bind_rows(diagnostic_rows),
       settings = list(units = units,
                       loci = loci,
                       n_rep = n_rep,
                       train_frac = train_frac,
                       balance = balance,
                       n_bal = if (balance) n_bal else NA_integer_,
                       prediction_priors = prediction_priors,
                       max_pcs = max_pcs,
                       seed = seed,
                       unit_sizes = sizes))
}


# Ejecutar todos los análisis con la misma función
cv_results <- list()

add_cv_run <- function(
    key,
    loci,
    units,
    panel_size,
    comparison,
    comparison_id,
    design,
    n_rep,
    seed) {
  
  message("Ejecutando: ", key)
  
  ans <- cv_dapc_fixed(gi = genind_top15,
                       loci = loci,
                       units = units,
                       n_rep = n_rep,
                       train_frac = train_fraction,
                       balance = design == "balanced",
                       prediction_priors = c("proportional", "uniform"),
                       seed = seed)
  
  metadata <- list(panel_size = panel_size,
                   comparison = comparison,
                   comparison_id = comparison_id,
                   design = design)
  
  ans$diagnostics <- ans$diagnostics |>
    dplyr::mutate(panel_size = metadata$panel_size,
                  comparison = metadata$comparison,
                  comparison_id = metadata$comparison_id,
                  design = metadata$design)
  
  if (nrow(ans$predictions)) {
    ans$predictions <- ans$predictions |>
      dplyr::mutate(panel_size = metadata$panel_size,
                    comparison = metadata$comparison,
                    comparison_id = metadata$comparison_id,
                    design = metadata$design)
  }
  
  ans$settings$metadata <- metadata
  
  # Checkpoint de cada ejecución.
  saveRDS(ans,
          file.path(out_dir, paste0(key, ".rds")))
  ans
}

# Cuatro unidades: 1,000 réplicas por diseño y panel.
# La misma semilla proporciona las mismas particiones para los tres paneles.
for (s in panel_sizes) {
  for (design in c("original", "balanced")) {
    
    key <- paste0("panel", s, "_four_units_", design)
    
    cv_results[[key]] <- add_cv_run(key = key,
                                    loci = selected_panels[[as.character(s)]]$loci,
                                    units = unit_order,
                                    panel_size = s,
                                    comparison = "four_units",
                                    comparison_id = "all_four_units",
                                    design = design,
                                    n_rep = cv_replicates,
                                    seed = base_seed)
  }
}

# Seis pares: 500 réplicas por diseño y panel.
unit_pairs <- combn(unit_order, 2, simplify = FALSE)

for (j in seq_along(unit_pairs)) {
  
  units_pair <- unit_pairs[[j]]
  pair_id <- paste(units_pair, collapse = "__")
  
  for (s in panel_sizes) {
    for (design in c("original", "balanced")) {
      
      key <- paste0("panel", s,
                    "_pair", j,
                    "_", design)
      
      cv_results[[key]] <- add_cv_run(key = key,
                                      loci = selected_panels[[as.character(s)]]$loci,
                                      units = units_pair,
                                      panel_size = s,
                                      comparison = "pairwise",
                                      comparison_id = pair_id,
                                      design,
                                      n_rep = pairwise_replicates,
                                      seed = base_seed + j)
    }
  }
}

predictions <- dplyr::bind_rows(lapply(cv_results, function(z) z$predictions))

diagnostics <- dplyr::bind_rows(lapply(cv_results, function(z) z$diagnostics))

saveRDS(predictions,
        file.path(out_dir, "individual_predictions.rds"))

write_csv(diagnostics,
          file.path(out_dir, "replicate_diagnostics.csv"))

replicate_status <- diagnostics |>
  dplyr::group_by(panel_size,
                  comparison,
                  comparison_id,
                  design,
                  prior) |>
  dplyr::summarise(requested_replicates = dplyr::n(),
                   successful_replicates = sum(status == "ok"),
                   failed_replicates = sum(status != "ok"),
                   replicates_with_warnings = sum(nzchar(warnings)),
                   .groups = "drop")

write_csv(replicate_status,
          file.path(out_dir, "replicate_status.csv"))

print(replicate_status)

if (any(replicate_status$failed_replicates > 0L)) {
  warning("Hay réplicas fallidas. Revisar Replicate_diagnostics.csv ",
          "antes de interpretar las tasas.")
}

if (!nrow(predictions)) {
  stop("No se obtuvo ninguna predicción válida.")
}

# Comparaciones entre paneles/diseños/priors:
# usar solamente réplicas válidas en todas las configuraciones
# del mismo ámbito y comparación.
required_configurations <- diagnostics |>
  dplyr::distinct(comparison,
                  comparison_id,
                  panel_size,
                  design,
                  prior) |>
  dplyr::count(comparison,
               comparison_id,
               name = "n_configurations")

common_replicates <- diagnostics |>
  dplyr::group_by(comparison,
                  comparison_id,
                  rep) |>
  dplyr::summarise(all_ok = all(status == "ok"),
                   n_observed = dplyr::n(),
                   .groups = "drop") |>
  dplyr::left_join(required_configurations,
                   by = c("comparison", "comparison_id")) |>
  dplyr::filter(all_ok,
                n_observed == n_configurations) |>
  dplyr::select(comparison,
                comparison_id,
                rep)

common_counts <- common_replicates |>
  dplyr::count(comparison,
               comparison_id,
               name = "common_successful_replicates")

write_csv(common_counts,
          file.path(out_dir, "common_successful_replicates.csv"))

print(common_counts)

expected_comparisons <- diagnostics |>
  dplyr::distinct(comparison, comparison_id)

missing_comparisons <- expected_comparisons |>
  dplyr::anti_join(common_counts,
                   by = c("comparison", "comparison_id"))

if (nrow(missing_comparisons)) {
  print(missing_comparisons)
  stop("Alguna comparación no tiene réplicas válidas comunes. ",
       "Revisar los diagnósticos.")
}

predictions_used <- predictions |>
  dplyr::semi_join(common_replicates,
                   by = c("comparison", "comparison_id", "rep"))

saveRDS(predictions_used,
        file.path(out_dir, "individual_predictions_common_replicates.rds"))

# Dimensiones realmente utilizadas en los ajustes válidos.
dimension_summary <- diagnostics |>
  dplyr::semi_join(common_replicates,
                   by = c("comparison", "comparison_id", "rep")) |>
  dplyr::filter(status == "ok",
                prior == "proportional") |>
  dplyr::group_by(panel_size,
                  comparison,
                  comparison_id,
                  design) |>
  dplyr::summarise(n_replicates = dplyr::n(),
                   variables_min = min(n_variables),
                   variables_max = max(n_variables),
                   pca_requested_min = min(n_pca_requested),
                   pca_requested_max = max(n_pca_requested),
                   pca_retained_min = min(n_pca_retained),
                   pca_retained_max = max(n_pca_retained),
                   da_retained_min = min(n_da_retained),
                   da_retained_max = max(n_da_retained),.groups = "drop")

write_csv(dimension_summary,
          file.path(out_dir, "retained_dimensions.csv"))


# Tasas por réplica, por unidad y globales
configuration_columns <- c("panel_size",
                           "comparison",
                           "comparison_id",
                           "design",
                           "prior")

per_unit_replicate <- predictions_used |>
  dplyr::group_by(dplyr::across(dplyr::all_of(c(configuration_columns, "rep", "true")))) |>
  dplyr::summarise(n_tests = dplyr::n(),
                   n_correct = sum(correct),
                   rate = mean(correct),
                   .groups = "drop")

per_unit_summary <- per_unit_replicate |>
  dplyr::group_by(dplyr::across(dplyr::all_of(c(configuration_columns, "true")))) |>
  dplyr::summarise(n_replicates = dplyr::n(),
                   total_tests = sum(n_tests),
                   total_correct = sum(n_correct),
                   mean_rate = mean(rate),
                   pooled_rate = sum(n_correct) / sum(n_tests),
                   sd_between_splits = sd(rate),
                   q025 = quantile(rate, 0.025,
                                   names = FALSE),
                   q975 = quantile(rate, 0.975,
                                   names = FALSE),
                   .groups = "drop")

pooled_replicate <- predictions_used |>
  dplyr::group_by(dplyr::across(dplyr::all_of(c(configuration_columns, "rep")))) |>
  dplyr::summarise(n_tests = dplyr::n(),
                   n_correct = sum(correct),
                   pooled_rate = mean(correct),
                   .groups = "drop")

macro_replicate <- per_unit_replicate |>
  dplyr::group_by(dplyr::across(dplyr::all_of(c(configuration_columns, "rep")))) |>
  dplyr::summarise(macro_rate = mean(rate),
                   .groups = "drop")

overall_replicate <- pooled_replicate |>
  dplyr::left_join(macro_replicate,
                   by = c(configuration_columns, "rep"))

overall_summary <- overall_replicate |>
  tidyr::pivot_longer(cols = c(pooled_rate,
                               macro_rate),
                      names_to = "metric",
                      values_to = "rate") |>
  dplyr::group_by(dplyr::across(dplyr::all_of(c(configuration_columns, "metric")))) |>
  dplyr::summarise(n_replicates = dplyr::n(),
                   mean_rate = mean(rate),
                   sd_between_splits = sd(rate),
                   q025 = quantile(rate, 0.025, names = FALSE),
                   q975 = quantile(rate, 0.975, names = FALSE),
                   .groups = "drop")

write_csv(per_unit_replicate,
          file.path(out_dir, "per_unit_rates_by_replicate.csv"))

write_csv(per_unit_summary,
          file.path(out_dir, "per_unit_assignment_summary.csv"))

write_csv(overall_replicate,
          file.path(out_dir, "overall_rates_by_replicate.csv"))

write_csv(overall_summary,
          file.path(out_dir, "overall_assignment_summary.csv"))

# Matrices de confusión acumuladas, incluidas celdas con cero predicciones

make_confusion <- function(z) {
  
  units <- unit_order[unit_order %in% unique(z$true)]
  
  counts <- as.data.frame(table(True = factor(z$true,
                                              levels = units),
                                Assigned = factor(z$assigned,
                                                  levels = units)),
                          stringsAsFactors = FALSE)
  
  counts |>
    dplyr::mutate(True = as.character(True),
                  Assigned = as.character(Assigned)) |>
    dplyr::group_by(True) |>
    dplyr::mutate(row_total = sum(Freq),
                  proportion = Freq / row_total) |>
    dplyr::ungroup()
}

confusion_long <- predictions_used |>
  dplyr::group_by(dplyr::across(dplyr::all_of(configuration_columns))) |>
  dplyr::group_modify(~ make_confusion(.x)) |>
  dplyr::ungroup()

readr::write_csv(confusion_long,
                 file.path(out_dir,
                           "confusion_matrices_long.csv"))

# Iberian north-south: directional false positives -----------------------------
# Norte-sur ibérico: falsos positivos direccionales

north <- "northern_iberian_atlantic"
south <- "southern_iberian_atlantic"
ns_id <- paste(c(north, south), collapse = "__")

ns_false_positive_replicate <- predictions_used |>
  dplyr::filter(comparison == "pairwise",
                comparison_id == ns_id,
                true == south) |>
  dplyr::group_by(panel_size,
                  design,
                  prior,
                  rep) |>
  dplyr::summarise(n_southern_tests = dplyr::n(),
                   n_assigned_north = sum(assigned == north),
                   false_positive_rate = mean(assigned == north),
                   .groups = "drop")

ns_false_positive_summary <- ns_false_positive_replicate |>
  dplyr::group_by(panel_size,
                  design,
                  prior) |>
  dplyr::summarise(n_replicates = dplyr::n(),
                   total_southern_tests = sum(n_southern_tests),
                   total_assigned_north = sum(n_assigned_north),
                   mean_false_positive_rate = mean(false_positive_rate),
                   pooled_false_positive_rate = sum(n_assigned_north) / sum(n_southern_tests),
                   sd_between_splits = sd(false_positive_rate),
                   q025 = quantile(false_positive_rate,
                                   0.025,
                                   names = FALSE),
                   q975 = quantile(false_positive_rate,
                                   0.975,
                                   names = FALSE),
                   .groups = "drop")

write_csv(ns_false_positive_replicate,
          file.path(out_dir,
                    "north_south_false_positives_by_replicate.csv"))

write_csv(ns_false_positive_summary,
          file.path(out_dir,
                    "north_south_false_positive_summary.csv"))

# Figure 2 ---------------------------------------------------------------------
# Figuras de las cuatro unidades

scenario_labels <- c(original_proportional = "Original sizes / proportional prior",
                     original_uniform = "Original sizes / uniform prior",
                     balanced_proportional = "Balanced sizes / proportional prior")

per_unit_plot_data <- per_unit_summary |>
  dplyr::filter(comparison == "four_units") |>
  dplyr::mutate(scenario = paste(design, prior, sep = "_")) |>
  dplyr::filter(scenario %in% names(scenario_labels)) |>
  dplyr::mutate(scenario = factor(scenario,
                                  levels = names(scenario_labels),
                                  labels = unname(scenario_labels)),
                true = factor(true, levels = rev(unit_order)),
                panel_label = factor(paste(panel_size, "SNPs"),
                                     levels = paste(panel_sizes, "SNPs")))

# Labels for the three scenarios
scenario_labels <- c(original_proportional = "Original sizes / proportional prior",
                     original_uniform = "Original sizes / uniform prior",
                     balanced_proportional = "Balanced sizes / proportional prior")

# Okabe-Ito colours for geographic units
unit_colours <- setNames(ggokabeito::palette_okabe_ito(order = c(5, 1, 3, 7)),
                         names(unit_labels)[match(c("Northern Iberian Atlantic",
                                                    "Southern Iberian Atlantic",
                                                    "Western Mediterranean",
                                                    "Macaronesia"),
                                                  unname(unit_labels))])

# Shapes for scenarios:
# circle = 16, square = 15, triangle = 17
scenario_shapes <- c("Original sizes / proportional prior" = 16,
                     "Original sizes / uniform prior" = 15,
                     "Balanced sizes / proportional prior" = 17)

p_per_unit <- ggplot(per_unit_plot_data,
                     aes(x = true,
                         y = mean_rate,
                         colour = true,
                         shape = scenario)) +
  geom_linerange(aes(ymin = q025,
                     ymax = q975),
                 position = position_dodge(width = 0.6),
                 linewidth = 0.6) +
  geom_point(position = position_dodge(width = 0.6),
             size = 3) +
  geom_hline(yintercept = 0.90,
             linetype = "dashed",
             linewidth = 0.7) +
  facet_wrap(~ panel_label) +
  scale_x_discrete(labels = function(z) unname(unit_labels[z])) +
  scale_shape_manual(values = scenario_shapes,
                     breaks = names(scenario_shapes),
                     name = "Scenario",
                     guide = guide_legend(order = 1)) +
  scale_colour_manual(values = unit_colours,
                      breaks = names(unit_colours),
                      labels = function(z) unname(unit_labels[z]),
                      name = "Geographic Unit",
                      guide = guide_legend(order = 2)) +
  coord_flip(ylim = c(0, 1)) +
  labs(x = NULL,
       y = "Correct assignment rate",
       caption = paste("Points: mean across repeated partitions.",
                       "Bars: empirical 2.5th–97.5th percentiles.")) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.box = "vertical")

ggsave(file.path(out_dir, "per_unit_assignment_by_scenario.png"),
       plot = p_per_unit,
       width = 15,
       height = 6,
       dpi = 300)

# Confusión con tamaños originales y prior proporcional
confusion_plot_data <- confusion_long |>
  dplyr::filter(comparison == "four_units",
                design == "original",
                prior == "proportional") |>
  dplyr::mutate(True = factor(True, levels = rev(unit_order)),
                Assigned = factor(Assigned, levels = unit_order),
                panel_label = factor(paste(panel_size, "SNPs"),
                                     levels = paste(panel_sizes, "SNPs")))

p_confusion <- ggplot(confusion_plot_data,
                      aes(Assigned, True, fill = proportion)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", proportion),
                colour = proportion > 0.5),
            size = 4) +
  scale_colour_manual(values = c("FALSE" = "black",
                                 "TRUE" = "white"),
                      guide = "none") +
  scale_fill_viridis_c(limits = c(0, 1),
                       option = "rocket",
                       direction = -1) +
  scale_x_discrete(labels = function(z) unname(unit_labels[z])) +
  scale_y_discrete(labels = function(z) unname(unit_labels[z])) +
  facet_wrap(~ panel_label) +
  labs(x = "Assigned origin",
       y = "True origin",
       fill = "Proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))

ggsave(file.path(out_dir, "confusion_matrices_original_proportional.png"),
       plot = p_confusion,
       width = 16,
       height = 6,
       dpi = 300)

# Combine Figure 2 panels
p_figure2 <- (p_per_unit / p_confusion) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(out_dir, "Figure_2.png"),
       plot = p_figure2,
       width = 16,
       height = 12,
       dpi = 300)

# Supplemental Figure 2 --------------------------------------------------------
# Figura por pares: tasa macro, paneles y diseños

pair_labels <- setNames(vapply(unit_pairs, function(z) {
    paste(unname(unit_labels[z]), collapse = " vs ")
  }, character(1)),
  vapply(unit_pairs, paste, character(1), collapse = "__"))

pairwise_plot_data <- overall_summary |>
  dplyr::filter(comparison == "pairwise",
                metric == "macro_rate",
                prior == "proportional") |>
  dplyr::mutate(comparison_id = factor(comparison_id,
                                       levels = names(pair_labels)),
                panel_label = factor(paste(panel_size, "SNPs"),
                                     levels = paste(panel_sizes, "SNPs")),
                design = factor(design,
                                levels = c("original",
                                           "balanced"),
                                labels = c("Original sizes",
                                           "Balanced sizes")))

design_colours <- c("Original sizes" = "#D55E00",
                    "Balanced sizes" = "#F0E442")

p_pairwise <- ggplot(pairwise_plot_data,
  aes(x = comparison_id,
      y = mean_rate,
      colour = design)) +
  geom_pointrange(aes(ymin = q025,
                      ymax = q975),
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~ panel_label) +
  scale_x_discrete(labels = function(z) unname(pair_labels[z])) +
  scale_colour_manual(values = design_colours,
                      breaks = names(design_colours),
                      name = "Scenario") +
  coord_flip(ylim = c(0, 1)) +
  labs(x = NULL,
       y = "Macro-average correct assignment rate",
       caption = paste("Points: mean across repeated partitions.",
                       "Bars: empirical 2.5th–97.5th percentiles.")) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.box = "vertical")

ggsave(file.path(out_dir, "pairwise_assignment_comparison.png"),
       plot = p_pairwise,
       width = 16,
       height = 8,
       dpi = 300)

# Guardar parámetros y mostrar los resúmenes principales -----------------------
analysis_settings <- list(cv_replicates = cv_replicates,
                          pairwise_replicates = pairwise_replicates,
                          train_fraction = train_fraction,
                          base_seed = base_seed,
                          panels = selected_panels,
                          sample_sizes = sample_sizes,
                          imputation = "Training-column means only",
                          dapc_center = TRUE,
                          dapc_scale = FALSE,
                          max_requested_pcs = 20L,
                          selection_repeated_within_training = FALSE,
                          summary_intervals = "Empirical 2.5th–97.5th percentiles across splits",
                          common_successful_replicates_used = TRUE)

saveRDS(analysis_settings,
        file.path(out_dir, "Analysis_settings.rds"))

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "Session_info.txt"))

cat("\n--- 4 UNIDADES: RESUMEN GLOBAL ---\n")
print(overall_summary |>
        dplyr::filter(comparison == "four_units"),
      n = Inf)

cat("\n--- FALSOS POSITIVOS: SUR ASIGNADO AL NORTE ---\n")
print(ns_false_positive_summary, n = Inf)

cat("\n--- DIMENSIONES RETENIDAS ---\n")
print(dimension_summary, n = Inf)

