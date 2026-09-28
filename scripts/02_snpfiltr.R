library(tidyverse)
library(janitor)
library(viridis)
library(vcfR)
library(SNPfiltR)

# Data import and exploration --------------------------------------------------
vcf <- read.vcfR("data/raw/final_octopus_UNO_120803_RAW_SNPs.vcf.recode.vcf")
vcf

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

# Individuos en el popmap que no están en el vcf
setdiff(popmap$id, vcf_inds)

# Calculo la intersección de individuos
inds_comunes <- intersect(vcf_inds,
                          popmap$id) # Deberían ser 477 individuos

# Filtro el VCF: Mantengo solo los individuos comunes
# La columna FORMAT se conserva con el índice 1, los individuos empiezan en 2
cols_keep <- c(TRUE, colnames(vcf@gt)[-1] %in% inds_comunes)

vcf@gt <- vcf@gt[, cols_keep]

# Filtro el popmap para que coincida exactamente
popmap <- popmap %>%
  filter(id %in% inds_comunes)

# Verifico
ncol(vcf@gt) - 1  # debe ser 477
nrow(popmap)      # debe ser 477

# Exporto en formato mPCRselect
popmap %>%
  dplyr::select(Sample = id, Population = pop) %>%
  write_csv("data/processed/popmap.csv",
            col_names = TRUE)

popmap %>%
  dplyr::select(id) %>%
  write_csv("data/processed/ids.csv",
            col_names = FALSE)

# Backup -----------------------------------------------------------------------
vcf_backup <- vcf
vcf <- vcf_backup

# Quality Filters --------------------------------------------------------------
# Implement quality filters that don’t involve missing data. This is because
# removing low data samples will alter percentage/quantile based missing data
# cutoffs, so we wait to implement those until after deciding on our final set
# of samples for downstream analysis.
# ---------------------------------------------------------------------------- #
head(vcf)

# Analisis exploratorio y visualización de las distribuciones de profundidad de
# secuenciación
png("results/figures/01_hard_filter_exploration.png",
    width = 1200,
    height = 800,
    res = 300)
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
vcf

# Eliminamos SNPs multi-alelicos
png("results/figures/02_biallelic.png",
    width = 1200,
    height = 800,
    res = 300)
vcf <- filter_biallelic(vcf)
dev.off()
vcf

# Then use this function to filter for allele balance
# ---------------------------------------------------------------------------- #
# Allele balance: a number between 0 and 1 representing the ratio of reads
# showing the reference allele to all reads, considering only reads from
# individuals called as heterozygous, we expect that the allele balance in our
# data (for real loci) should be close to 0.5
# ---------------------------------------------------------------------------- #
png("results/figures/03_allele_balance.png",
     width = 1200,
    height = 800,
    res = 300)
vcf <- filter_allele_balance(vcf)
dev.off()
vcf

# Visualize and pick appropriate max depth cutoff
# ---------------------------------------------------------------------------- #
# Now we can execute a max depth filter (super high depth loci are likely
# multiple loci stuck together into a single paralogous locus).
# ---------------------------------------------------------------------------- #
png("results/figures/04_max_depth_exploration.png",
    width = 1200,
    height = 800,
    res = 300)
max_depth(vcf) # dashed line indicates a mean depth across all SNPs of 400
dev.off()

dp <- extract.gt(vcf,
                 element = "DP",
                 as.numeric = TRUE)

median_dp <- median(dp,
                    na.rm = TRUE)
median_dp

mean_dp <- mean(dp,
                na.rm = TRUE)
mean_dp

quantile(dp,
         probs = c(0.5,
                   0.75,
                   0.9,
                   0.95,
                   0.99),
         na.rm = TRUE)

# percentil-based (95%)
max_depth_cutoff <- quantile(dp, 0.95,
                             na.rm = TRUE)

png("results/figures/05_max_depth.png",
    width = 1200,
    height = 800,
    res = 300)
vcf <- SNPfiltR::max_depth(vcf,
                           maxdepth = max_depth_cutoff)
dev.off()

# remove invariant SNPs generated during the genotype filtering steps
png("results/figures/06_invariant.png",
    width = 1200,
    height = 800,
    res = 300)
vcf <- min_mac(vcf, min.mac = 1)
dev.off()
vcf

# Missing data allowed per sample ----------------------------------------------
# ---------------------------------------------------------------------------- #
# Missing by sample mostrará si hay individuos con datos muy deficientes que
# conviene eliminar antes de continuar con el filtro de MAF.
# ---------------------------------------------------------------------------- #
png("results/figures/07_missing_by_sample_exploration.png",
    width = 1400,
    height = 900,
    res = 300)
missing_by_sample(vcf = vcf,
                  popmap = popmap) # Exploración
dev.off()

ids_antes <- colnames(vcf@gt)[-1]

png("results/figures/08_missing_by_sample_cutoff50.png",
    width = 1400,
    height = 900,
    res = 300)
vcf_50 <- missing_by_sample(vcf,
                            cutoff = 0.5)
dev.off()
vcf_50

png("results/figures/09_missing_by_sample_cutoff40.png",
    width = 1400,
    height = 900,
    res = 300)
vcf_40 <- missing_by_sample(vcf,
                            cutoff = 0.4)
dev.off()
vcf_40

png("results/figures/10_missing_by_sample_cutoff30.png",
    width = 1400,
    height = 900,
    res = 300)
vcf_30 <- missing_by_sample(vcf,
                            cutoff = 0.3)
dev.off()
vcf_30

png("results/figures/11_missing_by_sample_cutoff20.png",
    width = 1400,
    height = 900,
    res = 300)
vcf_20 <- missing_by_sample(vcf,
                            cutoff = 0.2)
dev.off()
vcf_20

removed_inds_50 <- setdiff(colnames(vcf@gt)[-1], colnames(vcf_50@gt)[-1])
removed_inds_40 <- setdiff(colnames(vcf@gt)[-1], colnames(vcf_40@gt)[-1])
removed_inds_30 <- setdiff(colnames(vcf@gt)[-1], colnames(vcf_30@gt)[-1])
removed_inds_20 <- setdiff(colnames(vcf@gt)[-1], colnames(vcf_20@gt)[-1])

popmap %>%
  filter(id %in% removed_inds_50) %>%
  count(pop)

popmap %>%
  filter(id %in% removed_inds_40) %>%
  count(pop)

popmap %>%
  filter(id %in% removed_inds_30) %>%
  count(pop)

popmap %>%
  filter(id %in% removed_inds_20) %>%
  count(pop)

removed_summary <- bind_rows(popmap %>%
                               filter(id %in% removed_inds_20) %>%
                               count(pop) %>%
                               mutate(threshold = "20%"),
                             popmap %>%
                               filter(id %in% removed_inds_30) %>%
                               count(pop) %>%
                               mutate(threshold = "30%"),
                             popmap %>%
                               filter(id %in% removed_inds_40) %>%
                               count(pop) %>%
                               mutate(threshold = "40%"),
                             popmap %>% filter(id %in% removed_inds_50) %>%
                               count(pop) %>%
                               mutate(threshold = "50%"))

p1 <- ggplot(removed_summary,
             aes(x = pop,
                 y = n)) +
  geom_col(fill = "steelblue") +
  facet_wrap(~threshold) +
  coord_flip() +
  labs(x = "Population",
       y = "Number of removed individuals") +
  theme_minimal()

p1

ggsave("results/figures/12_removed_individuals_bars.png",
       plot = p1,
       width = 10,
       height = 10,
       dpi = 300)

p2 <- ggplot(removed_summary,
             aes(x = threshold,
                 y = pop,
                 fill = n)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "C") +
  labs(x = "Missing data threshold",
       y = "Population",
       fill = "Removed individuals") +
  theme_minimal()

p2

ggsave("results/figures/13_removed_individuals_heatmap.png",
       plot = p2,
       width = 10,
       height = 10,
       dpi = 300)

rm(vcf_50,
   vcf_40,
   vcf_30,
   vcf_20)

# Uso un cutoff 0.3. Northern Ibernian Atlantic aguanta bien y no hay ninguna
# justificación para sacrificar 30 muestras adicionales de esa población por 0.4%
# de mejora global.

png("results/figures/14_missing_by_sample_final_cutoff30.png",
    width = 1400,
    height = 900,
    res = 300)
vcf <- missing_by_sample(vcf = vcf,
                         cutoff = 0.3) # cutoff
dev.off()

ids_despues <- colnames(vcf@gt)[-1]

muestras_eliminadas <- setdiff(ids_antes,
                               ids_despues)
print(muestras_eliminadas)

# Subset popmap to only include retained individuals
popmap <- popmap %>%
  filter(id %in% colnames(vcf@gt)[-1])

nrow(popmap)
ncol(vcf@gt) - 1

# Remove invariant sites generated by dropping individuals
png("results/figures/15_min_mac_exploration.png",
    width = 1200,
    height = 800,
    res = 300)
min_mac(vcf)
dev.off()

vcf <- min_mac(vcf,
               min.mac = 1)
vcf
# 0.21% of SNPs fell below a minor allele count of 1 and were removed from the VCF

# Verificar que el missing no crea estructura artificial
assess_missing_data_pca(vcf = vcf,
                        popmap = popmap,
                        thresholds = c(0.5, 0.65, 0.8, 0.95),
                        clustering = FALSE)

# Missing data allowed per SNP -------------------------------------------------
# Decido cuántos SNPs con muchos datos faltantes estoy dispuesta a sacrificar
# para quedarme con un dataset más completo por SNP, sin perder demasiados
# marcadores
png("results/figures/19_missing_by_snp_exploration.png",
    width = 1200,
    height = 800,
    res = 300)
missing_by_snp(vcf) # Exploracion
dev.off()

# Choose a value that retains an acceptable amount of missing data in each
# sample, and maximizes SNPs retained while minimizing overall missing data,
# and filter vcf
png("results/figures/18_missing_by_snp_cutoff50.png",
    width = 1200,
    height = 800,
    res = 300)
vcf <- missing_by_snp(vcf,
                      cutoff = 0.85) # Cutoff
dev.off()
vcf

# Verificar que el missing no crea estructura artificial
assess_missing_data_pca(vcf = vcf,
                        popmap = popmap,
                        thresholds = c(0.65, 0.75, 0.85, 0.95),
                        clustering = FALSE)

# check what t-SNE clustering looks like at an 85% threshold
assess_missing_data_tsne(vcf = vcf,
                         popmap = popmap,
                         thresholds = c(0.65, 0.75, 0.85, 0.95),
                         clustering = FALSE)

vcf

# Effect of a minor allele count (MAC) cutoff ----------------------------------
png("results/figures/19_min_mac_post_snp_filter.png",
    width = 1200,
    height = 800,
    res = 300)
min_mac(vcf) # Exploracion
dev.off()

png("results/figures/20_tsne_before_mac.png",
    width = 1600,
    height = 1000,
    res = 300)
assess_missing_data_tsne(vcf,
                         popmap,
                         clustering = FALSE)
dev.off()

# Por ahora aplicamos un MAC mínimo conservador solo para eliminar singletons y
# errores y después del pipeline de SNPfiltR, el filtro realmente crítico para
# trazabilidad será seleccionar SNPs por Fst entre poblaciones
vcf.mac <- min_mac(vcf,
                   min.mac = 2) # Cutoff

pdf("results/figures/16_tsne_after_mac.pdf",
    width = 1600,
    height = 1000)
assess_missing_data_tsne(vcf.mac,
                         popmap,
                         clustering = FALSE)
dev.off()

# ---------------------------------------------------------------------------- #
# Finally, we will make sure that the depth look consistent across SNPs and
# samples, following our filtering pipeline.
# ---------------------------------------------------------------------------- #
# Plot depth per snp and per sample
dp <- extract.gt(vcf.mac,
                 element = "DP",
                 as.numeric = TRUE)

pdf("results/figures/17_depth_heatmap_final.pdf",
    width = 2000,
    height = 1200)
heatmap.bp(dp,
           rlabels = FALSE)
dev.off()

# Write out vcf files for downstream analysis ----------------------------------
# Write out vcf with all SNPs
vcfR::write.vcf(vcf.mac,
                "data/processed/octopus.filtered.vcf.gz")
