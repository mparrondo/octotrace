# ============================================================================ #
# Genome-wide SNP discovery and a reduced diagnostic panel for geographic
# assignment of common octopus (Octopus vulgaris) in the Bay of Biscay and
# adjacent fishing regions
# 
# Data Analysis
# Author: Marina Parrondo Lombardía (parrondomarina@proton.me)
# ============================================================================ #
# ============================================================================ #
# SNP filtering with SNPfiltR // Filtrado de SNPs con SNPfiltR 
#
# Workflow (DeRaad 2022, SNPfiltR vignette):
# A. Quality filters excluding missing data (depth, biallelic, allelic balance,
#    maximum depth, invariant sites) //
#    Filtros de calidad sin missing data (profundidad, bialélicos, balance
#    alélico, profundidad máxima, invariantes)
# B. Exploration of population structure by locality (prior to sample-based
#    filtering) //
#    Exploración de estructura por localidad (antes del filtrado por muestra)
# C. Missing data per sample -> threshold selection -> exclusion of samples
#    with high missingness //
#    Missing data por muestra -> elección del umbral -> exclusión EM
# D. Evaluation of missing data per SNP (PCA/t-SNE using various thresholds) ->
#    threshold selection -> filtering //
#    Evaluación del missing por SNP (PCA/t-SNE con varios umbrales) ->
#    elección del umbral -> filtro
# E. Minor Allele Count (MAC) filtering and assessment of its effect on
#    population structure//
#    Filtro MAC y comprobación de su efecto sobre la estructura
# F. Final checks (missing data per sample, depth, PCA by locality for the
#    Iberian North Atlantic region using the final dataset) //
#    Comprobaciones finales (missing por muestra, profundidad, PCA por localidad
#    en el Atlántico norte ibérico con el dataset final)
# G. Export with versioned filenames + log of counts //
#    Exportación con nombres versionados + log de recuentos
#
# Hardy-Weinberg equilibrium filtering is not applied (as it may remove loci
# under selection; Pearman et al. 2022). LD-based thinning and pruning are
# performed in mPCRselección.
# No se aplica filtro de equilibrio Hardy-Weinberg (puede eliminar loci bajo
# selección; Pearman et al. 2022). El thinning y la poda por LD se hacen en
# mPCRselect.
# ============================================================================ #

suppressPackageStartupMessages({
library(tidyverse)
library(janitor)
library(viridis)
library(vcfR)
library(SNPfiltR)
})

FIG <- "results/figures/"
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

fig <- function(name) file.path(FIG, name)

# Data import and exploration --------------------------------------------------
vcf <- read.vcfR("data/raw/final_octopus_UNO_120803_RAW_SNPs.vcf.recode.vcf")
vcf

# Recording of counts at each step (Filtering table for Methods) ---------------
# Registro de recuentos en cada paso (Tabla de filtrado para Métodos) -------- #
filter_log <- tibble()
log_step <- function(step, v) {
  gt <- extract.gt(v)
  filter_log <<- bind_rows(filter_log,
                           tibble(step = step,
                                  n_snps = nrow(v@fix),
                                  n_ind = ncol(v@gt) - 1,
                                  missing_p = round(100 * mean(is.na(gt)), 2)))
  invisible(v)
}

# Voy a limpiar el archivo de poblaciones porque tiene inconsistencias,
# redundancias y clasificaciones con las que no estoy de acuerdo por imprecisas
pops_raw <- read_tsv("data/processed/pops_complete.txt",
                     col_names = c("id",
                                   "locality_year",
                                   "year",
                                   "locality",
                                   "region",
                                   "pop"),
                     show_col_types = FALSE)

# Paso todo a minusculas
pops <- pops_raw |>
  mutate(across(-id, tolower)) |>
  clean_names()

# Modifico localidades para ajustar precisión. Me baso en emails de 2020-21.
pops <- pops |>
  mutate(locality = case_when(locality == "basquecountry" ~ "pasaia",
                              locality %in% c("puerto_vega",
                                              "puerta_vega",
                                              "puertovega") ~ "puerto_vega",
                              locality == "galicia" ~ "bueu",
                              locality == "portugal" &
                                locality_year == "portugal_2021_olhao" ~ "olhao",
                              locality == "portugal" &
                                locality_year == "portugal_2021_tavira" ~ "tavira",
                              locality == "portugal" ~ "algarve",
                              locality == "barcelona" ~ "deltebre",
                              locality == "turquia" ~ "karaburun",
                              locality == "canarias" ~ "san_andres",
                              locality == "ghana" ~ "tema",
                              .default = locality))
# Modifico regiones
pops <- pops |>
  mutate(region = case_when(region == "cantabrian" &
                              locality_year == "basquecountry_2021" ~ "guipuzkoa",
                            region == "cantabrian" &
                              locality_year == "galicia_2021" ~ "pontevedra",
                            region == "cantabrian" ~ "asturias",
                            region == "atlantic" &
                              locality == "olhao" ~ "algarve",
                            region == "atlantic" &
                              locality == "tavira" ~ "algarve",
                            region == "atlantic" &
                              locality == "algarve" ~ "algarve",
                            region == "mediterranean" &
                              locality == "deltebre" ~ "tarragona",
                            region == "mediterranean" &
                              locality == "karaburun" ~ "izmir",
                            region == "atlantic" &
                              locality == "san_andres" ~ "tenerife",
                            region == "gulg_guinea" ~ "greater_accra",
                            .default = region))

# Modifico poblaciones para que tengan sentido biológico
pops <- pops |>
  mutate(pop = case_when(pop == "east_cantabrian" ~ "northern_iberian_atlantic",
                         pop == "central_cantabrian" ~ "northern_iberian_atlantic",
                         pop == "west_cantabrian" ~ "northern_iberian_atlantic",
                         pop == "central_atlantic" ~ "southern_iberian_atlantic",
                         pop == "mediterranean" &
                           region == "tarragona" ~ "western_mediterranean",
                         pop == "mediterranean" &
                           region == "izmir" ~ "eastern_mediterranean",
                         pop == "south_atlantic" ~ "macaronesia",
                         pop == "gulg_guinea" ~ "gulf_guinea",
                         .default = pop))

# Modifico años, porque son claramente errores. Me baso en emails de 2020-21.
pops <- pops |>
  mutate(year = as.numeric(year),
         year = case_when(year == 2000 &
                            region == "izmir" ~ 2020,
                          year == 2000 &
                            region == "greater_accra" ~ 2021,
                          TRUE ~ year),
         locality_year = paste(locality, year, sep = "_"))

head(pops)

# Guardo
write.csv(pops,
          "data/processed/pops.csv",
          row.names = FALSE)

# Genero archivo popmap. Se trata de un archivo con dos columnas semejante a los
# archivos stacks, las columnas deben nombrarse como "id" y "pop"

# Corro primero el script a nivel de locality para comprobar que no hay estructura
# en el Cantabrico
# popmap <- pops |>
#   dplyr::select(id, pop = locality)

popmap <- pops |>
  dplyr::select(id, pop)

# Verifico
head(popmap)
nrow(popmap)

vcf_inds <- colnames(vcf@gt)[-1]   # -1 porque la primera columna es FORMAT

# Individuos en el VCF que no están en el popmap
setdiff(vcf_inds, popmap$id)

# Individuos en el popmap que no están en el vcf (ind. Olhao 2021)
setdiff(popmap$id, vcf_inds)

# Calculo la intersección de individuos
inds_comunes <- intersect(vcf_inds,
                          popmap$id) # Deberían ser 477 individuos

# Filtro el VCF: Mantengo solo los individuos comunes
# La columna FORMAT se conserva con el índice 1, los individuos empiezan en 2
vcf@gt <- vcf@gt[, c(TRUE, colnames(vcf@gt)[-1] %in% inds_comunes)]

# Filtro el popmap para que coincida exactamente
popmap <- popmap %>%
  filter(id %in% inds_comunes)

# Verifico
ncol(vcf@gt) - 1  # debe ser 477
nrow(popmap)      # debe ser 477

stopifnot(ncol(vcf@gt) - 1 == 477, nrow(popmap) == 477)
log_step("Raw (individuals in popmap)", vcf)

# Backup -----------------------------------------------------------------------
vcf_backup <- vcf
vcf <- vcf_backup

# A. Quality Filters -----------------------------------------------------------
# Implement quality filters that don’t involve missing data. This is because
# removing low data samples will alter percentage/quantile based missing data
# cutoffs, so we wait to implement those until after deciding on our final set
# of samples for downstream analysis.
# ---------------------------------------------------------------------------- #
head(vcf)

# Analisis exploratorio y visualización de las distribuciones de profundidad de
# secuenciación
pdf(fig("01_hard_filter_exploration.pdf"),
    width = 10,
    height = 7)
hard_filter(vcf)
dev.off()

# ---------------------------------------------------------------------------- #
# The VCF file does not include the GQ (Genotype Quality) field in the FORMAT
# for each genotype, so hard_filter() cannot apply that filter. The solution is
# simple: we simply omit the gq argument and apply only the depth filter.
# ---------------------------------------------------------------------------- #
# Hard filter to minimum depth of 5
vcf <- hard_filter(vcf = vcf,
                   depth = 5)
log_step("Genotype depth >= 5", vcf)
vcf

# Eliminamos SNPs multi-alelicos
pdf(fig("02_biallelic.pdf"),
    width = 10,
    height = 7)
vcf <- filter_biallelic(vcf)
dev.off()
log_step("Biallelic SNPs", vcf)
vcf

# Then use this function to filter for allele balance
# ---------------------------------------------------------------------------- #
# Allele balance: a number between 0 and 1 representing the ratio of reads
# showing the reference allele to all reads, considering only reads from
# individuals called as heterozygous, we expect that the allele balance in our
# data (for real loci) should be close to 0.5
# ---------------------------------------------------------------------------- #
pdf(fig("03_allele_balance.pdf"),
    width = 10,
    height = 7)
vcf <- filter_allele_balance(vcf)
dev.off()
log_step("Allele balance (0.25-0.75)", vcf)
vcf

# Visualize and pick appropriate max depth cutoff
# ---------------------------------------------------------------------------- #
# Now we can execute a max depth filter (super high depth loci are likely
# multiple loci stuck together into a single paralogous locus).
# ---------------------------------------------------------------------------- #
pdf(fig("04_max_depth_exploration.pdf"),
    width = 10,
    height = 7)
max_depth(vcf) # dashed line indicates a mean depth across all SNPs of 361.9
dev.off()

# FORMAT/DP: filas = SNPs; columnas = muestras
dp <- vcfR::extract.gt(vcf,
                       element = "DP",
                       as.numeric = TRUE)

# Profundidad media por SNP: misma variable que utiliza max_depth()
dp_snp <- rowMeans(dp, na.rm = TRUE)

# Resumen de la variable
quantile(dp_snp,
         probs = c(0.5, 0.75, 0.9, 0.95, 0.99),
         na.rm = TRUE)

# Percentil 95 de las medias por SNP
max_depth_cutoff <- unname(quantile(dp_snp,
                                    probs = 0.95,
                                    na.rm = TRUE))

# Filtrar y guardar la figura
pdf(fig("05_max_depth_filter.pdf"),
    width = 10,
    height = 7)
vcf <- SNPfiltR::max_depth(vcf, maxdepth = max_depth_cutoff)
dev.off()
log_step(sprintf("Mean SNP depth <= %.0f (95th pct)", max_depth_cutoff), vcf)
vcf

# remove invariant SNPs generated during the genotype filtering steps
pdf(fig("06_invariant_after_quality.pdf"),
    width = 10,
    height = 7)
vcf <- min_mac(vcf, min.mac = 1)
dev.off()
log_step("Remove invariant (MAC >= 1)", vcf)
vcf

# B. Exploration by location (prior to filtering by sample) --------------------

# Backup después de calidad y antes de eliminar individuos
vcf_all_samples <- vcf
popmap_all_samples <- popmap

# Popmap por localidad, ordenado según el VCF
all_ids <- colnames(vcf_all_samples@gt)[-1]

popmap_locality <- pops |>
  dplyr::select(id, pop = locality) |>
  dplyr::filter(id %in% all_ids)

stopifnot(!anyDuplicated(all_ids),
          !anyDuplicated(popmap_locality$id),
          setequal(all_ids, popmap_locality$id),
          !anyNA(popmap_locality$pop))

popmap_locality <- popmap_locality |>
  dplyr::slice(match(all_ids, id))

thresholds <- c(0.65, 0.75, 0.85, 0.95)

pdf(fig("07_pca_locality_all_before_sample_filter.pdf"),
    width = 10,
    height = 7)
pca_locality_before <- assess_missing_data_pca(vcf_all_samples,
                                               popmap = popmap_locality,
                                               thresholds = thresholds,
                                               clustering = FALSE)
dev.off()

# Localidades y colores comunes a los cuatro PCA
localities <- sort(unique(unlist(lapply(pca_locality_before, function(x) x$pop))))

palette_locality <- ggsci::pal_d3("category20")(20)
palette_locality <- palette_locality[c(seq(1, 19, 2),
                                       seq(2, 20, 2))]

stopifnot(length(localities) <= length(palette_locality))

colors_locality <- setNames(palette_locality[seq_along(localities)],
                            localities)

plot_locality_pcas <- function(pca_list, file, title) {
  stopifnot(length(pca_list) == length(thresholds))
  pdf(file,
      width = 10,
      height = 7)
  for (i in seq_along(pca_list)) {
    p <- ggplot(pca_list[[i]],
                aes(PC1, PC2, color = pop)) +
      geom_point(size = 2.3,
                 alpha = 0.9) +
      scale_color_manual(values = colors_locality,
                         limits = intersect(localities, unique(pca_list[[i]]$pop)),
                         name = "Locality") +
      labs(title = title,
           subtitle = paste0("Minimum SNP completeness: ", thresholds[i] * 100, "%"),
           x = "PC1",
           y = "PC2") +
      theme_classic(base_size = 12) +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
    print(p)
  }
  dev.off()
}

## Supplementary Figure S1 -----------------------------------------------------
plot_locality_pcas(pca_locality_before,
                   fig("08_pca_locality_all_before_sample_filter_custom.pdf"),
                   "PCA by locality — before sample filtering")

# Atlántico norte ibérico: exploración inicial (la definitiva va en F)
make_nia_subset <- function(v) {
  ids <- colnames(v@gt)[-1]
  nia <- pops |>
    filter(pop == "northern_iberian_atlantic", id %in% ids) |>
    pull(id)
  v@gt <- v@gt[, c(TRUE, ids %in% nia), drop = FALSE]
  v <- min_mac(v, min.mac = 1)
  ord <- colnames(v@gt)[-1]
  pm <- popmap_locality |>
    filter(id %in% ord) |>
    slice(match(ord, id))
  list(vcf = v, popmap = pm)
}

nia0 <- make_nia_subset(vcf_all_samples)
pdf(fig("09_pca_locality_nia_before_sample_filter.pdf"),
    width = 10,
    height = 7)
pca_locality_nia_before <- assess_missing_data_pca(nia0$vcf,
                                                   popmap = nia0$popmap,
                                                   thresholds = thresholds,
                                                   clustering = FALSE)
dev.off()

# C. Missing data allowed per sample -------------------------------------------
# ---------------------------------------------------------------------------- #
# Missing by sample mostrará si hay individuos con datos muy deficientes que
# conviene eliminar antes de continuar con el filtro de MAF.
# ---------------------------------------------------------------------------- #
pdf(fig("10_missing_by_sample_exploration.pdf"),
    width = 10,
    height = 7)
missing_by_sample(vcf = vcf,
                  popmap = popmap)
dev.off()

ids_antes <- colnames(vcf@gt)[-1]

cutoffs <- c(0.5, 0.4, 0.3, 0.2)

removed <- list()

for (k in seq_along(cutoffs)) {
  pdf(fig(sprintf("%02d_missing_by_sample_cutoff%02d.pdf",
                  10 + k, cutoffs[k] * 100)),
      width = 10,
      height = 7)
  v_tmp <- missing_by_sample(vcf, cutoff = cutoffs[k])
  dev.off()
  removed[[paste0(cutoffs[k] * 100, "%")]] <- setdiff(ids_antes, colnames(v_tmp@gt)[-1])
  print(v_tmp)
}

rm(v_tmp)

removed_summary <- imap_dfr(removed, function(ids, thr)
  popmap |>
    filter(id %in% ids) |>
    count(pop) |>
    mutate(threshold = thr))
print(removed_summary)

p1 <- ggplot(removed_summary,
             aes(pop,
                 n)) +
  geom_col(fill = "steelblue") +
  facet_wrap(~threshold) +
  coord_flip() +
  labs(x = "Population",
       y = "Number of removed individuals") +
  theme_classic()

ggsave(fig("15_removed_individuals_bars.png"),
       p1,
       width = 10,
       height = 10,
       dpi = 300)

p2 <- ggplot(removed_summary,
             aes(threshold,
                 pop,
                 fill = n)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "C") +
  labs(x = "Missing data threshold",
       y = "Population",
       fill = "Removed individuals") +
  theme_classic()

ggsave(fig("16_removed_individuals_heatmap.png"),
       p2,
       width = 10,
       height = 10,
       dpi = 300)

# Uso un cutoff 0.3. Northern Ibernian Atlantic aguanta bien y no hay ninguna
# justificación para sacrificar muestras adicionales de esa población por 0.4%
# de mejora global.

pdf(fig("17_missing_by_sample_final_cutoff30.pdf"),
    width = 10,
    height = 7)
vcf <- missing_by_sample(vcf = vcf,
                         cutoff = 0.3)
dev.off()
log_step("Sample missingness <= 30%", vcf)

muestras_eliminadas <- setdiff(ids_antes,
                               colnames(vcf@gt)[-1])
popmap_all_samples |>
  filter(id %in% muestras_eliminadas) |>
  write_csv("results/tables/removed_samples_missing30.csv")

# Exclusión del Mediterráneo oriental ------------------------------------------
# ---------------------------------------------------------------------------- #
# Tras el filtro por muestra (cutoff 0.3) solo queda 1 de los 3 individuos de
# esta unidad. Se excluye aquí, ANTES de los filtros por sitio (min_mac,
# missing_by_snp), por tres motivos:
# a) Con n = 1 la unidad no admite partición entrenamiento/test en la validación
# cruzada ni una estima de frecuencias alélicas poblacionales.
# b) Un único individuo infla el Fst por locus y condicionaría la selección de
# marcadores.
# c) Excluirlo antes de min_mac evita conservar SNPs cuya variación dependa solo
# de este individuo.
# ---------------------------------------------------------------------------- #
em_ids <- popmap$id[popmap$pop == "eastern_mediterranean" &
                      popmap$id %in% colnames(vcf@gt)[-1]]
print(em_ids)
stopifnot(length(em_ids) == 1)

# Copia con el individuo
vcf_with_em <- vcf
saveRDS(vcf_with_em,
        "data/processed/vcf_with_em_after_sample_filter_4u.rds")

vcf <- vcf[, !(colnames(vcf@gt) %in% em_ids)]
stopifnot(ncol(vcf@gt) - 1 == 452)

writeLines(em_ids, "results/tables/excluded_eastern_mediterranean.txt")
log_step("Exclude Eastern Mediterranean (n = 1)", vcf)

# Subset popmap to only include retained individuals
popmap <- popmap %>%
  filter(id %in% colnames(vcf@gt)[-1])

# Solo deben quedar las 4 unidades de referencia (ni EM, ni Gulf of Guinea, ni
# muestras de mercado, que no forman parte de la línea base)
units_4u <- c("northern_iberian_atlantic",
              "southern_iberian_atlantic",
              "western_mediterranean",
              "macaronesia")

stopifnot(nrow(popmap) == ncol(vcf@gt) - 1,
          setequal(popmap$id,
                   colnames(vcf@gt)[-1]),
          setequal(unique(popmap$pop), units_4u))

print(table(popmap$pop)) # 360 / 45 / 23 / 24

pdf(fig("18_invariant_after_sample_filter.pdf"),
    width = 10,
    height = 7)
vcf <- min_mac(vcf,
               min.mac = 1)
dev.off()
log_step("Remove invariant after sample removal", vcf)
vcf

# D. Missing data allowed per SNP ----------------------------------------------
# Backup después del filtrado por individuo y antes del filtrado por SNP
vcf_before_site_filter <- vcf
popmap_before_site_filter <- popmap

# La evaluación con PCA/t-SNE se hace ANTES de filtrar: se comparan varios
# umbrales sobre el mismo VCF para elegir uno.
snp_thresholds <- c(0.65, 0.75, 0.85, 0.95)

pdf(fig("19_pca_missing_snp_thresholds.pdf"),
    width = 10,
    height = 7)
pca_snp_thr <- assess_missing_data_pca(vcf,
                                       popmap = popmap,
                                       thresholds = snp_thresholds,
                                       clustering = FALSE)
dev.off()

pdf(fig("20_tsne_missing_snp_thresholds.pdf"),
    width = 10,
    height = 7)
tsne_snp_thr <- assess_missing_data_tsne(vcf,
                                         popmap = popmap,
                                         thresholds = snp_thresholds,
                                         clustering = FALSE)
dev.off()

pdf(fig("21_missing_by_snp_exploration.pdf"),
    width = 10,
    height = 7)
missing_by_snp(vcf)
dev.off()

pdf(fig("22_missing_by_snp_cutoff85.pdf"),
    width = 10,
    height = 7)
vcf <- missing_by_snp(vcf,
                      cutoff = 0.85)
dev.off()
log_step("SNP completeness >= 85%", vcf)

# E. Minor allele count --------------------------------------------------------

pdf(fig("23_min_mac_exploration.pdf"),
    width = 10,
    height = 7)
min_mac(vcf)
dev.off()

pdf(fig("24_tsne_before_mac.pdf"),
    width = 10,
    height = 7)
assess_missing_data_tsne(vcf,
                         popmap,
                         clustering = FALSE)
dev.off()

# Por ahora aplicamos un MAC mínimo conservador solo para eliminar singletons y
# errores y después del pipeline de SNPfiltR, el filtro realmente crítico para
# trazabilidad será seleccionar SNPs por Fst entre poblaciones
vcf.mac <- min_mac(vcf,
                   min.mac = 2)
log_step("MAC >= 2", vcf.mac)

pdf(fig("25_tsne_after_mac.pdf"),
    width = 10,
    height = 7)
assess_missing_data_tsne(vcf.mac,
                         popmap,
                         clustering = FALSE)
dev.off()

# F. Comprobaciones finales ----------------------------------------------------
# ---------------------------------------------------------------------------- #
# Finally, we will make sure that the depth look consistent across SNPs and
# samples, following our filtering pipeline.
# ---------------------------------------------------------------------------- #
# Plot depth per snp and per sample
dp <- extract.gt(vcf.mac,
                 element = "DP",
                 as.numeric = TRUE)

pdf(fig("26_depth_heatmap_final.pdf"),
    width = 10,
    height = 7)
heatmap.bp(dp,
           rlabels = FALSE)
dev.off()

# Missing por muestra en el dataset final (para reportar en Resultados)
pdf(fig("27_missing_by_sample_final_dataset.pdf"),
    width = 10,
    height = 7)
missing_by_sample(vcf = vcf.mac,
                  popmap = popmap)
dev.off()

miss_final <- colMeans(is.na(extract.gt(vcf.mac)))
tibble(id = names(miss_final),
       missing = miss_final) |>
  left_join(popmap,
            by = "id") |>
  group_by(pop) |>
  summarise(n = n(),
            mean_missing = mean(missing),
            max_missing = max(missing)) |>
  write_csv("results/tables/missing_per_unit_final.csv")

# PCA por localidad en el Atlántico norte ibérico con el dataset FINAL.
# Es la que respalda en el manuscrito la agrupación de País Vasco, Asturias
# y Galicia en una sola unidad.
nia_f <- make_nia_subset(vcf.mac)

pdf(fig("28_pca_locality_nia_final_dataset.pdf"),
    width = 10,
    height = 7)
pca_locality_nia_final <- assess_missing_data_pca(nia_f$vcf,
                                                  popmap = nia_f$popmap,
                                                  thresholds = NULL,
                                                  clustering = FALSE)
dev.off()



# G. Write out files for downstream analysis -----------------------------------
# Ordeno el popmap según los individuos del VCF final
final_ids <- colnames(vcf.mac@gt)[-1]

stopifnot(!anyDuplicated(final_ids),
          !anyDuplicated(popmap$id),
          setequal(final_ids, popmap$id))

popmap_final <- popmap |>
  dplyr::slice(match(final_ids, id))

# Exporto en formato mPCRselect
popmap_final |>
  dplyr::select(Sample = id,
                Population = pop) |>
  write_csv("data/processed/popmap.csv")

popmap_final |>
  dplyr::select(id) |>
  write_csv("data/processed/ids.csv",
            col_names = FALSE)

# Write out vcf files for downstream analysis ----------------------------------
# Write out vcf with all SNPs
out_vcf <- sprintf("data/processed/octopus.filtered.4u.%dsnps.%dind.vcf.gz",
                   nrow(vcf.mac@fix),
                   ncol(vcf.mac@gt) - 1)
vcfR::write.vcf(vcf.mac, out_vcf)

filter_log <- filter_log |>
  mutate(snps_removed = lag(n_snps) - n_snps,
         ind_removed = lag(n_ind) - n_ind)

print(filter_log, n = Inf)
write_csv(filter_log, "results/tables/snpfiltr_filter_log_4u.csv")

writeLines(c(paste("file:", out_vcf),
             paste("md5:", tools::md5sum(out_vcf)),
             paste("snps:", nrow(vcf.mac@fix)),
             paste("individuals:", ncol(vcf.mac@gt) - 1),
             paste("date:", Sys.time())),
           sub("\\.vcf\\.gz$", ".manifest.txt", out_vcf))

sessionInfo()
