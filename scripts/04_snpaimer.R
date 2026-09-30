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
})

# Write out vcf files for downstream analysis ----------------------------------
# Generar archivos VCF para análisis posteriores ----------------------------- #

# Load the VCF output from mPCRselect
# Cargo el vcf que sale de mPCRselect
vcf_fst_pi <- read.vcfR("data/processed/octopus.fst_pi.vcf.gz")

# Convert to Genind/ Convierte a genind
genind <- vcfR2genind(vcf_fst_pi)

# Assigns populations in the same order as the VCF
# Asigna poblaciones en el mismo orden que el VCF
pop(genind) <- popmap$pop[match(indNames(genind), popmap$id)]

# Verify that there are no NAs (individuals without an assignment)
# Verifica que no hay NAs (individuos sin asignación)
table(is.na(pop(genind)))

# Preselect top SNPs by Fst for snpAIMeR ---------------------------------------
# Preseleccionar los SNP top según el valor Fst para snpAIMeR ---------------- #

# Convert genind to hierfstat format
# Convertir genind a formato hierfstat 
hf_data <- genind2hierfstat(genind)

# Calculate locus-specific statistics (includes locus-specific Fst)
# Calcular estadísticos por locus (incluye Fst por locus)
bs <- basic.stats(hf_data, diploid = TRUE)

# Extract Fst by locus and sort
# Extraer Fst por locus y ordenar
fst_df <- data.frame(locus = locNames(genind),
                     Fst   = bs$perloc$Fst) |>
  dplyr::arrange(dplyr::desc(Fst))

# Top 15 SNPs by Fst
# Top 15 SNPs por Fst
print(head(fst_df, 15))

# Subsett genind to top 15
# Subsetear genind a top 15
top15 <- fst_df$locus[1:15]
genind_top15 <- genind[loc = top15]

# PCA of the top 15 SNPs
# PCA de los top 15 SNPs
x <- tab(genind_top15,
         NA.method="mean")

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

# Export and correct ':' in marker names.
# The first row (header with ind/pop) is removed because read.structure
# does not accept headers with metadata columns—it would count them as extra loci.
#
# Exportar y corregir ':' en nombres de marcadores.
# Se elimina la primera fila (header con ind/pop) porque read.structure
# no acepta encabezados con columnas de metadata — los contaría como loci extra.
genind2structure(genind_top15,
                 file = "data/processed/octopus_top15.str",
                 pops = TRUE)

lines <- readLines("data/processed/octopus_top15.str")

lines[1] <- gsub(":", "_", lines[1])

writeLines(lines[-1],
           "data/processed/octopus_top15_fixed.str")

rm(lines)

# Quick check / Verificación rápida
str_top15 <- read.table("data/processed/octopus_top15_fixed.str",
                        header = FALSE,
                        sep = "\t")
# Rows / Filas (n.ind x 2)
nrow(str_top15)
# Columns / Columnas
ncol(str_top15)
# Markers / Marcadores
ncol(str_top15) - 2

rm(str_top15)

# Configuration / Configuración
config <- list(min_range = 1L,
               max_range = 15L,
               assignment_rate_threshold = 0.9,
               cross_validation_replicates = 1000L,
               working_directory = normalizePath("results/snpAIMeR/"),
               structure_file = normalizePath("data/processed/octopus_top15_fixed.str"),
               number_of_individuals = 453L,
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
snpAIMeR("non-interactive",
         "data/processed/snpAIMeR_config.yaml",
         verbose = TRUE)

all_comb <- read.csv("results/snpAIMeR/All_combinations_assign_rate.csv")
above <- read.csv("results/snpAIMeR/Above_threshold_assign_rate.csv")
panel <- read.csv("results/snpAIMeR/Panel_size_assign_rate.csv")

# Combinations evaluated / Combinaciones evaluadas
nrow(all_comb)
# Combinations above the 0.9 threshold
# Combinaciones sobre umbral 0.9
nrow(above)
# APanel size / Tamaño de panel
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

# Read the header of the original file (before removing the header)
# Leer el header del archivo original (antes de quitar el header)
header <- readLines("data/processed/octopus_top15.str", n = 1)
marker_names <- strsplit(header, "\t")[[1]]
marker_names <- gsub(":", "_", marker_names)  # mismo fix que apliqué antes

# V3=L01, V4=L02, etc. → posiciones 3:17
loci_names <- marker_names[3:17]
names(loci_names) <- paste0("L", sprintf("%02d", 1:15))

best_8 <- all_comb %>%
  filter(panel_size == 8, !is.na(avg_success_rate)) %>%
  slice_max(avg_success_rate,
            n = 1,
            with_ties = FALSE)

# Separate the markers from that combination
# Separar los marcadores de esa combinación
panel_8 <- trimws(strsplit(best_8$marker[[1]], ",", fixed = TRUE)[[1]])

print(panel_8)
print(best_8$avg_success_rate)

# Convert L01...L15 to the original locus names
# Convertir L01...L15 a los nombres originales de los loci
stopifnot(all(panel_8 %in% names(loci_names)))

# Optimal panel / Panel óptimo
print(loci_names[panel_8])

# Subset of the genind object containing the 8 SNPs from the optimal panel
# Subset del genind con los 8 SNPs del panel óptimo
panel_loci <- c("RXHP01002566_1:16165",
                "RXHP01002823_1:48716",
                "RXHP01004691_1:201139",
                "RXHP01005384_1:80616",
                "RXHP01006845_1:6481",
                "RXHP01009923_1:271862",
                "RXHP01018693_1:683",
                "RXHP01028750_1:34890")

genind_panel8 <- genind_top15[, loc = panel_loci]

# PCA
pca_panel8 <- dudi.pca(tab(genind_panel8,
                           NA.method = "mean"),
                       scannf = FALSE,
                       nf = 3)

# Data frame para ggplot
pca_df <- data.frame(PC1 = pca_panel8$li[, 1],
                     PC2 = pca_panel8$li[, 2],
                     pop = pop(genind_panel8))

# Explained variance / Varianza explicada
var_exp <- round(pca_panel8$eig / sum(pca_panel8$eig) * 100, 1)

levels(pop(genind_panel8))

pop_labels <- c("Northern Iberian Atlantic",
                "Southern Iberian Atlantic",
                "Eastern Mediterranean",
                "Western Mediterranean",
                "Macaronesia")

p_pca8 <- ggplot(pca_df,
                 aes(x = PC1,
                     y = PC2,
                     color = pop)) +
  geom_point(alpha = 0.7,
             size = 3) +
  scale_color_okabe_ito(labels = pop_labels,
                        order  = c(5, 1, 7, 3, 6)) +
  labs(title = "8 SNPs",
       x     = paste0("PC1 (", var_exp[1], "%)"),
       y     = paste0("PC2 (", var_exp[2], "%)"),
       color = "Population") +
  theme_classic()

p_pca8

ggsave("results/snpAIMeR/PCA_panel8.png",
       plot   = p_pca8,
       width  = 9,
       height = 6,
       dpi    = 300)

# Function for PCA of a SNP panel/ Función para PCA de un panel de SNPs
run_pca_panel <- function(genind_obj, panel_loci, panel_name, pop_labels) {
  
  # Subset loci
  genind_panel <- genind_obj[, loc = panel_loci]
  
  # PCA
  pca <- dudi.pca(tab(genind_panel, NA.method = "mean"),
                  scannf = FALSE,
                  nf = 3)
  
  # Data frame para ggplot
  pca_df <- data.frame(PC1 = pca$li[, 1],
                       PC2 = pca$li[, 2],
                       pop = pop(genind_panel))
  
  # Explained variance / Varianza explicada
  var_exp <- round(pca$eig / sum(pca$eig) * 100, 1)
  
  # Plot
  p <- ggplot(pca_df,
              aes(x = PC1,
                  y = PC2,
                  color = pop)) +
    geom_point(alpha = 0.7,
               size = 3) +
    scale_color_okabe_ito(labels = pop_labels,
                          order  = c(5, 1, 7, 3, 6)) +
    labs(title = paste0("Panel ", panel_name),
         x     = paste0("PC1 (", var_exp[1], "%)"),
         y     = paste0("PC2 (", var_exp[2], "%)"),
         color = "Population") +
    theme_classic()
  # Save / Guardar
  ggsave(paste0("results/snpAIMeR/PCA_panel", panel_name, ".png"),
         plot = p,
         width = 9,
         height = 6,
         dpi = 300)
  
  return(list(genind = genind_panel,
              pca = pca,
              plot = p,
              variance = var_exp))
}

panel_loci_3 <- c("RXHP01005384_1:80616",
                  "RXHP01021799_1:6426",
                  "RXHP01067935_1:8708")

panel_loci_8 <- c("RXHP01002566_1:16165",
                  "RXHP01002823_1:48716",
                  "RXHP01004691_1:201139",
                  "RXHP01005384_1:80616",
                  "RXHP01006845_1:6481",
                  "RXHP01009923_1:271862",
                  "RXHP01018693_1:683",
                  "RXHP01028750_1:34890")

panel_loci_11 <- c("RXHP01001794_1:112749",
                   "RXHP01002566_1:16165",
                   "RXHP01002823_1:48716",
                   "RXHP01003976_1:15402",
                   "RXHP01004691_1:201139",
                   "RXHP01005384_1:80616",
                   "RXHP01006845_1:6481",
                   "RXHP01009923_1:271862",
                   "RXHP01012430_1:17539",
                   "RXHP01018693_1:683",
                   "RXHP01027781_1:638213")

pca3 <- run_pca_panel(genind_top15,
                      panel_loci_3,
                      "3",
                      pop_labels)

pca8 <- run_pca_panel(genind_top15,
                      panel_loci_8,
                      "8",
                      pop_labels)

pca11 <- run_pca_panel(genind_top15,
                       panel_loci_11,
                       "11",
                       pop_labels)

# Figure 2 ---------------------------------------------------------------------

fig_pca <- (pca3$plot | pca8$plot | pca11$plot) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
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

# Dispersión, matriz de confusión y asignación por pares -----------------------
# Complemento al análisis de snpAIMeR

# Error de la tasa de asignación -----------------------------------------------
# La dispersión que pide Trini ya está en All_combinations_assign_rate.csv:
# panel es la media de avg_success_rate por tamaño de panel, y al agrupar la
# tabla completa se recupera la variabilidad que esa media colapsa

if (!"panel_size" %in% names(all_comb)) {
  all_comb$panel_size <- sapply(strsplit(all_comb$marker, ", "), length)
}

panel_err <- all_comb %>%
  group_by(panel_size) %>%
  summarise(n_comb = n(),
            mean = mean(avg_success_rate),
            sd = sd(avg_success_rate),
            se = sd / sqrt(n_comb),
            q025 = quantile(avg_success_rate, 0.025),
            q975 = quantile(avg_success_rate, 0.975),
            min = min(avg_success_rate),
            max = max(avg_success_rate),
            .groups = "drop")

print(as.data.frame(panel_err),
      digits = 3)

write.csv(panel_err,
          file.path("results/snpAIMeR/Panel_size_dispersion.csv"),
          row.names = FALSE)

# Comprobación: panel_err$mean debe coincidir con panel$avg_success_rate
# Si no coincide, snpAIMeR pondera de otra forma y hay que revisarlo.
print(data.frame(snpAIMeR = panel$avg_success_rate,
                 recalculado = panel_err$mean,
                 dif = panel$avg_success_rate - panel_err$mean))

# Curva con banda de dispersión
p_err <- ggplot(panel_err,
                aes(panel_size, mean)) +
  geom_ribbon(aes(ymin = q025,
                  ymax = q975),
              alpha = 0.18) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_point(data = all_comb %>%
               group_by(panel_size) %>%
               summarise(best = max(avg_success_rate),
                         .groups = "drop"),
             aes(panel_size, best),
             shape = 17,
             size = 2.4) +
  geom_hline(yintercept = 0.9,
             linetype = "dashed",
             linewidth = 0.4) +
  scale_x_continuous(breaks = 1:15) +
  labs(x = "Panel size (number of loci)",
       y = "Correct assignment rate",
       caption = paste("Circles: mean across all combinations of a given size.",
                       "Shaded band: 2.5-97.5% quantiles.",
                       "Triangles: best-performing combination.",
                       "Dashed line: 0.90 threshold.")) +
  theme_classic()

p_err

ggsave(file.path("results/snpAIMeR/Panel_size_with_dispersion.png"),
       p_err,
       width = 8,
       height = 5,
       dpi = 300)

# ¿Es el panel de 11 mejor que el de 8? Comparación explícita
comp <- all_comb %>%
  filter(panel_size %in% c(8, 11))
print(t.test(avg_success_rate ~ panel_size,
             data = comp))
print(wilcox.test(avg_success_rate ~ panel_size,
                  data = comp))
# Ojo: son combinaciones no independientes (comparten loci), así que estos
# contrastes son descriptivos. Lo defendible es reportar medias, sd y solape
# de los rangos, no un p-valor.

# Validación cruzada replicada con DAPC: matriz de confusión -------------------
# Reproduce el esquema de snpAIMeR (75/25 por población, DAPC) pero guarda
# la asignación individuo a individuo, que es lo que snpAIMeR no exporta.

cv_dapc <- function(gi, loci, n_rep = 1000, train_frac = 0.75,
                    n_pca = NULL, seed = 1) {
  
  gi <- gi[, loc = loci]
  set.seed(seed)
  
  pops <- pop(gi)
  n_by_pop <- table(pops)
  
  # Grupos con menos de 4 individuos no admiten un reparto 75/25 informativo
  too_small <- names(n_by_pop)[n_by_pop < 4]
  if (length(too_small)) {
    message("Unidades excluidas de la CV por n < 4: ",
            paste(too_small, collapse = ", "))
  }
  keep <- !(as.character(pops) %in% too_small)
  gi <- gi[keep, ]
  pops <- factor(pop(gi))
  
  if (is.null(n_pca)) n_pca <- min(ncol(tab(gi)) - 1, 20)
  n_da <- max(1, nlevels(pops) - 1)
  
  out <- vector("list", n_rep)
  
  for (r in seq_len(n_rep)) {
    idx_train <- unlist(lapply(levels(pops), function(p) {
      ii <- which(pops == p)
      sample(ii, size = max(2, floor(train_frac * length(ii))))
    }))
    idx_test <- setdiff(seq_along(pops), idx_train)
    if (!length(idx_test)) next
    
    gi_tr <- gi[idx_train, ]
    gi_te <- gi[idx_test, ]
    
    d <- tryCatch(
      dapc(gi_tr, pop = factor(pop(gi_tr)),
           n.pca = min(n_pca, nInd(gi_tr) - 1), n.da = n_da),
      error = function(e) NULL)
    if (is.null(d)) next
    
    pr <- tryCatch(predict.dapc(d, newdata = gi_te),
                   error = function(e) NULL)
    if (is.null(pr)) next
    
    out[[r]] <- data.frame(
      rep      = r,
      ind      = indNames(gi_te),
      true     = as.character(pop(gi_te)),
      assigned = as.character(pr$assign),
      post_max = apply(pr$posterior, 1, max),
      stringsAsFactors = FALSE
    )
  }
  
  loo <- bind_rows(out)
  loo$correct <- loo$true == loo$assigned
  loo
}

loo8 <- cv_dapc(genind_top15,
                panel_loci_8,
                n_rep = 1000)

# Matriz de confusión (proporciones por fila = origen real)
conf <- table(True = loo8$true,
              Assigned = loo8$assigned)

conf_prop <- round(prop.table(conf,
                              margin = 1), 4)

print(conf_prop)

write.csv(as.data.frame.matrix(conf_prop),
          file.path("results/snpAIMeR/Confusion_matrix_panel8.csv"))

# Tasa por unidad con IC binomial (Wilson vía prop.test)
per_unit <- loo8 %>%
  group_by(true) %>%
  summarise(n_tests = n(),
            n_correct = sum(correct),
            .groups = "drop") %>%
  rowwise() %>%
  mutate(rate  = n_correct / n_tests,
         ci_lo = prop.test(n_correct, n_tests)$conf.int[1],
         ci_hi = prop.test(n_correct, n_tests)$conf.int[2]) %>%
  ungroup()

print(as.data.frame(per_unit),
      digits = 4)

write.csv(per_unit,
          file.path("results/snpAIMeR/Per_unit_assignment_panel8.csv"),
          row.names = FALSE)

# Media global (ponderada por individuo) y macro-media (por unidad)
# Ponderada
round(mean(loo8$correct), 4)

# Macro-media no ponderada por tamaño de grupo
round(mean(per_unit$rate), 4)
# El IC binomial aquí es optimista: las réplicas reutilizan los mismos
# individuos, así que las observaciones no son independientes

# Asignación por pares de poblaciones ------------------------------------------
# Reajusta el DAPC dentro de cada par, que es lo que pide Trini

pairwise_assign <- function(gi, loci, n_rep = 500) {
  units <- levels(factor(pop(gi)))
  n_by_pop <- table(pop(gi))
  units <- units[n_by_pop[units] >= 4]
  pares <- combn(units, 2, simplify = FALSE)
  
  bind_rows(lapply(pares, function(p) {
    sub <- gi[as.character(pop(gi)) %in% p, ]
    pop(sub) <- factor(as.character(pop(sub)))
    res <- cv_dapc(sub,
                   loci,
                   n_rep = n_rep)
    n_ok <- sum(res$correct);
    n_t <- nrow(res)
    data.frame(unit_a = p[1],
               unit_b = p[2],
               n_a = sum(pop(sub) == p[1]),
               n_b = sum(pop(sub) == p[2]),
               rate = n_ok / n_t,
               ci_lo = prop.test(n_ok, n_t)$conf.int[1],
               ci_hi = prop.test(n_ok, n_t)$conf.int[2])
  }))
}

pw8 <- pairwise_assign(genind_top15,
                       panel_loci_8,
                       n_rep = 500)

print(as.data.frame(pw8),
      digits = 4)

write.csv(pw8,
          file.path("results/snpAIMeR/Pairwise_assignment_panel8.csv"),
          row.names = FALSE)

# El par que interesa
print(pw8 %>%
        filter(grepl("orth", unit_a) | grepl("orth", unit_b)))

# Mismo cálculo para los paneles de 3 y 11, para la tabla suplementaria
pw3  <- pairwise_assign(genind_top15,
                        panel_loci_3, 
                        n_rep = 500)

pw11 <- pairwise_assign(genind_top15,
                        panel_loci_11,
                        n_rep = 500)

pw_all <- bind_rows(pw3 %>%
                      mutate(panel = "3 loci"),
                    pw8 %>%
                      mutate(panel = "8 loci"),
                    pw11 %>%
                      mutate(panel = "11 loci"))

write.csv(pw_all,
          file.path("results/snpAIMeR/Pairwise_assignment_all_panels.csv"),
          row.names = FALSE)

# Mapa de calor por pares (panel de 8)
p_pw <- ggplot(pw8,
               aes(unit_a,
                   unit_b,
                   fill = rate)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.3f", rate)),
            size = 3.2) +
  scale_fill_viridis_c(limits = c(0.5, 1),
                       option = "mako",
                       direction = -1) +
  labs(x = NULL,
       y = NULL,
       fill = "Correct\nassignment") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1))

p_pw

ggsave(file.path("results/snpAIMeR/Pairwise_assignment_panel8.png"),
       p_pw,
       width = 7.5,
       height = 5.5,
       dpi = 300)

# Efecto del Mediterráneo oriental (n = 1) -------------------------------------
# Comprueba cuánto cambian las tasas globales al excluir esa unidad.

east <- "eastern_mediterranean" # ver levels(pop(genind_top15))
gi_no_east <- genind_top15[as.character(pop(genind_top15)) != east, ]
pop(gi_no_east) <- factor(as.character(pop(gi_no_east)))

loo8_no_east <- cv_dapc(gi_no_east,
                        panel_loci_8,
                        n_rep = 1000)

# Con Mediterráneo oriental:
round(mean(loo8$correct), 4)
# Sin Mediterráneo oriental:
round(mean(loo8_no_east$correct), 4)

# Efecto del desequilibrio de tamaños muestrales en la asignación --------------

EAST <- "eastern_mediterranean"

# Tamaños muestrales:
print(table(pop(genind_top15)))

# cv_dapc con prior uniforme y con submuestreo equilibrado
cv_dapc2 <- function(gi, loci, n_rep = 1000, train_frac = 0.75,
                     prior = c("proportional", "uniform"),
                     balance = FALSE, n_bal = NULL, seed = 1) {
  
  prior <- match.arg(prior)
  gi <- gi[, loc = loci]
  set.seed(seed)
  
  pops <- factor(pop(gi))
  n_by_pop <- table(pops)
  too_small <- names(n_by_pop)[n_by_pop < 4]
  if (length(too_small))
    message("Excluidas por n < 4: ", paste(too_small, collapse = ", "))
  gi <- gi[!(as.character(pops) %in% too_small), ]
  pops <- factor(pop(gi))
  
  # Si balance = TRUE, en cada réplica se submuestrea cada población al
  # tamaño del grupo más pequeño, de modo que el prior sea efectivamente plano
  if (is.null(n_bal)) n_bal <- min(table(pops))
  
  n_pca <- min(ncol(tab(gi)) - 1, 20)
  out <- vector("list", n_rep)
  
  for (r in seq_len(n_rep)) {
    
    if (balance) {
      idx_use <- unlist(lapply(levels(pops), function(p)
        sample(which(pops == p), n_bal)))
    } else {
      idx_use <- seq_along(pops)
    }
    p_use <- factor(pops[idx_use])
    
    idx_train_rel <- unlist(lapply(levels(p_use), function(p) {
      ii <- which(p_use == p)
      sample(ii, size = max(2, floor(train_frac * length(ii))))
    }))
    idx_test_rel <- setdiff(seq_along(idx_use), idx_train_rel)
    if (!length(idx_test_rel)) next
    
    gi_tr <- gi[idx_use[idx_train_rel], ]
    gi_te <- gi[idx_use[idx_test_rel], ]
    p_tr  <- factor(pop(gi_tr))
    
    d <- tryCatch(
      dapc(gi_tr, pop = p_tr,
           n.pca = min(n_pca, nInd(gi_tr) - 1),
           n.da  = max(1, nlevels(p_tr) - 1)),
      error = function(e) NULL)
    if (is.null(d)) next
    
    pr <- tryCatch(predict.dapc(d, newdata = gi_te), error = function(e) NULL)
    if (is.null(pr)) next
    
    # Prior uniforme: se reponderan las posteriores dividiendo por la
    # frecuencia del grupo en el entrenamiento y se renormaliza
    if (prior == "uniform") {
      w <- as.numeric(table(p_tr)[colnames(pr$posterior)])
      post <- sweep(pr$posterior, 2, w, "/")
      post <- post / rowSums(post)
      assigned <- colnames(post)[apply(post, 1, which.max)]
      pmax_    <- apply(post, 1, max)
    } else {
      assigned <- as.character(pr$assign)
      pmax_    <- apply(pr$posterior, 1, max)
    }
    
    out[[r]] <- data.frame(rep = r, ind = indNames(gi_te),
                           true = as.character(pop(gi_te)),
                           assigned = assigned, post_max = pmax_,
                           stringsAsFactors = FALSE)
  }
  
  res <- bind_rows(out)
  res$correct <- res$true == res$assigned
  res
}

# Tres escenarios comparados
esc <- list(proporcional = cv_dapc2(genind_top15,
                                    panel_loci_8,
                                    prior = "proportional"),
            uniforme = cv_dapc2(genind_top15,
                                panel_loci_8,
                                prior = "uniform"),
            equilibrado = cv_dapc2(genind_top15,
                                   panel_loci_8,
                                   prior = "proportional",
                          balance = TRUE))

resumen <- bind_rows(lapply(names(esc), function(nm) {
  d <- esc[[nm]]
  pu <- d %>% group_by(true) %>%
    summarise(rate = mean(correct), .groups = "drop")
  data.frame(escenario = nm,
             ponderada = mean(d$correct),
             macro     = mean(pu$rate))
}))

print(resumen, digits = 4)

# Matriz de confusión por escenario
for (nm in names(esc)) {
  cat("\n---", nm, "---\n")
  print(round(prop.table(table(True = esc[[nm]]$true,
                               Assigned = esc[[nm]]$assigned), 1), 3))
}

# Tasas por unidad en los tres escenarios
per_unit_all <- bind_rows(lapply(names(esc), function(nm) {
  esc[[nm]] %>%
    group_by(true) %>%
    summarise(n_tests = n(),
              n_correct = sum(correct),
              .groups = "drop") %>%
    rowwise() %>%
    mutate(rate = n_correct / n_tests,
           ci_lo = prop.test(n_correct, n_tests)$conf.int[1],
           ci_hi = prop.test(n_correct, n_tests)$conf.int[2],
           escenario = nm) %>%
    ungroup()
}))

print(as.data.frame(per_unit_all),
      digits = 4)

write.csv(per_unit_all,
          file.path("results/snpAIMeR/Per_unit_by_prior.csv"),
          row.names = FALSE)

p_cmp <- ggplot(per_unit_all,
                aes(reorder(true, rate),
                    rate,
                    colour = escenario)) +
  geom_pointrange(aes(ymin = ci_lo,
                      ymax = ci_hi),
                  position = position_dodge(width = 0.45)) +
  geom_hline(yintercept = 0.9,
             linetype = "dashed",
             linewidth = 0.4) +
  coord_flip() +
  ylim(0.4, 1) +
  labs(x = NULL,
       y = "Correct assignment rate",
       colour = NULL) +
  theme_classic()

p_cmp

ggsave(file.path("results/snpAIMeR/Per_unit_by_prior.png"),
       p_cmp,
       width = 8,
       height = 4.5,
       dpi = 300)

# Iberian north-south pair with balanced sizes ---------------------------------
# Par norte-sur ibérico con tamaños equilibrados ----------------------------- #
pair_ns <- c("northern_iberian_atlantic",
             "southern_iberian_atlantic")
gi_ns <- genind_top15[as.character(pop(genind_top15)) %in% pair_ns, ]
pop(gi_ns) <- factor(as.character(pop(gi_ns)))

ns_prop <- cv_dapc2(gi_ns,
                    panel_loci_8,
                    prior = "proportional")
ns_bal  <- cv_dapc2(gi_ns,
                    panel_loci_8,
                    prior = "proportional",
                    balance = TRUE)

# North vs. South, original sizes (360/45)
# Norte vs sur, tamaños originales (360/45)
round(mean(ns_prop$correct), 4)
print(round(prop.table(table(ns_prop$true,
                             ns_prop$assigned), 1), 3))

# North vs. South, downsampled to 45/45
# Norte vs sur, submuestreado a 45/45
round(mean(ns_bal$correct), 4)
print(round(prop.table(table(ns_bal$true, ns_bal$assigned), 1), 3))

# False positive rate for certified origin: individuals from the south
# assigned to the north. This is the relevant figure for traceability.
#
# Tasa de falso positivo para el origen certificado: individuos del sur
# asignados al norte. Es la cifra relevante para trazabilidad.
fp <- ns_bal %>%
  filter(true == pair_ns[2]) %>%
  summarise(n = n(),
            fp = sum(assigned == pair_ns[1])) %>%
  mutate(rate = fp / n,
         ci_lo = prop.test(fp, n)$conf.int[1],
         ci_hi = prop.test(fp, n)$conf.int[2])

print(as.data.frame(fp),
      digits = 4)