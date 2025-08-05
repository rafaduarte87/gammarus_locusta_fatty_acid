# Loading the packages ---------------------------------------------------------

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(performance)
library(emmeans)
library(patchwork)
library(vegan)
library(ggfortify)
library(ggcorrplot)
library(ggtext)
library(pairwiseAdonis)
library(ape)
library(ade4)
library(multcompView)

# Importing the datasets -------------------------------------------------------

# Macroalgae most important FA (more than 1% of relative abundance)

dataset_macroalgae <-
  read_excel("gammarus_dataset.xlsx", sheet = 1,
             range = "A1:T26", na = c("NA", "")) |> 
  mutate(algae = factor(algae,
                        levels = c("ulva", "fucus", "laminaria",
                                   "gracilaria", "chondrus")))

# Amphipods most important FA (more than 1% of relative abundance)

dataset_amphipods <-
  read_excel("gammarus_dataset.xlsx", sheet = 2,
             range = "A1:V31", na = c("NA", "")) |> 
  mutate(algae = factor(algae,
                        levels = c("ulva", "fucus", "laminaria",
                                    "gracilaria", "chondrus", "wild")))

# 1 - Macroalgae ---------------------------------------------------------------

## PERMANOVA -------------------------------------------------------------------

# testing if the FA profiles significantly differ among the five macroalgal species

# separating the numerical data and the factors in different matrices

data_fa_macroalgae <-
  dataset_macroalgae |>
  dplyr::select(`14:0` : `20:5 n-3`)

data_fa_macroalgae_log <-
  dataset_macroalgae |>
  dplyr::select(`14:0` : `20:5 n-3`) |> 
  mutate(across(`14:0` : `20:5 n-3`, ~ log(.x + 1)))

factor_macroalgae <-
  dataset_macroalgae |> 
  dplyr::select(algae : replicate)

# running the test

with(factor_macroalgae, 
     adonis2(data_fa_macroalgae_log ~ algae,
             method = "euc", by = "terms", permutations = 9999,
             data = factor_macroalgae))

# pairwise comparison tests

pairwise.adonis2(
  vegdist(data_fa_macroalgae_log, 
          method = "euclidean") ~ algae,
  data = as.data.frame(dataset_macroalgae))

# SIMPER analysis 

simper_algae <- 
  with(factor_macroalgae, 
       simper(data_fa_macroalgae_log, algae))

summary(simper_algae) 

## PCoA ------------------------------------------------------------------------

# creating a distance matrix

distance_fa_macroalgae <-
  vegdist(data_fa_macroalgae_log, method = "euclidean")

# running the PCoA

pcoa_fa_macroalgae <- 
  pcoa(distance_fa_macroalgae)

pcoa_fa_macroalgae$values

# extracting the PCoA vectors

pcoa_fa_macroalgae_arrows <- 
  compute_arrows(pcoa_fa_macroalgae, data_fa_macroalgae_log)

pcoa_fa_macroalgae_arrows$vectors <- 
  as.data.frame(pcoa_fa_macroalgae_arrows$vectors)

dataset_macroalgae$axis_1 <- pcoa_fa_macroalgae_arrows$vectors[,1]

dataset_macroalgae$axis_2 <- pcoa_fa_macroalgae_arrows$vectors[,2]

# scaling the vectors for a better visualization

arrows_df_macroalgae <- as.data.frame(pcoa_fa_macroalgae_arrows$U * 1)

arrows_df_macroalgae$variable <- rownames(arrows_df_macroalgae)

# plotting the PCoA

plot_pcoa_macroalgae_a <- 
  ggplot() +
  geom_point(data = dataset_macroalgae, 
             aes(axis_1, axis_2, fill = algae), 
             size = 4, shape = 21) +
  labs(x = "**PCO1 (85.5 % of total variation)**",
       y = "**PCO2 (5.8 % of total variation)**",
       fill = "**Macroalgae**")

plot_pcoa_macroalgae_b <-
  plot_pcoa_macroalgae_a +
  scale_x_continuous(breaks = seq(-1, 2, by = 0.5),
                     limits = c(-1, 2)) +
  scale_y_continuous(breaks = seq(-1, 0.75, by = 0.25),
                     limits = c(-1, 0.75)) +
  scale_fill_manual(values = c("#70AD47", "#FFCC00", "#993300",
                               "#FF6666", "#FF0000"),
                    labels = c("*Ulva*", "*Fucus*", "*Laminaria*", 
                               "*Gracilaria*", "*Chondrus*")) 

plot_pcoa_macroalgae_c <-
  plot_pcoa_macroalgae_b +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                                        margin = 
                                        margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12, 
                                        margin = 
                                        margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "inside",
        legend.position.inside = c(0.23, 0.88),
        legend.title = element_markdown(size = 10),
        legend.text = element_markdown(size = 8),
        legend.spacing.x = unit(0.25, "cm"),
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot_pcoa_macroalgae_d <-
  plot_pcoa_macroalgae_c +
  guides(fill = guide_legend(nrow = 2),
                theme = theme(legend.byrow = TRUE)) + 
  geom_segment(data = as.data.frame(pcoa_fa_macroalgae_arrows$U * 1),
               x = 0, y = 0, 
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(2, "mm"))) +
  geom_text(data = arrows_df_macroalgae, 
            aes(x = Axis.1, y = Axis.2, label = variable),
            position = position_jitter(), 
            size = 4)

# saving the graph

ggsave("PCoA FA macroalgae.png", plot_pcoa_macroalgae_d,
       height = 12.5, width = 17.5, units = "cm", dpi = 600)

## PCA -------------------------------------------------------------------------

pca_FA_macroalgae <-
  prcomp(data_fa_macroalgae_log,
         center = TRUE, scale. = TRUE)

summary(pca_FA_macroalgae)

# extract PC axes

pca_FA_macroalgae_values <-
  data.frame(dataset_macroalgae$algae,
             dataset_macroalgae$replicate,
             pca_FA_macroalgae$x)

colnames(pca_FA_macroalgae_values)[1] <-"algae"

colnames(pca_FA_macroalgae_values)[2] <-"replicate"

pca_FA_macroalgae_values <-
  pca_FA_macroalgae_values  |> 
  as_tibble()

plot_pca_FA_macroalgae <-
  autoplot(pca_FA_macroalgae, data = dataset_macroalgae, 
           fill = 'algae', shape = 21, size = 4,
           loadings = TRUE, loadings.colour = 'black',
           loadings.label = TRUE,
           loadings.label.size = 3,
           loadings.label.colour = "black") +
  labs(x = "**PC1 (61.39 % of total variation)**",
       y = "**PC2 (18.90 % of total variation)**",
       fill = "**Macroalgae**") +
  scale_fill_manual(values = c("#70AD47", "#FFCC00", "#993300",
                               "#FF6666", "#FF0000"),
                    labels = c("*Ul*", "*Fu*", "*La*", "*Gr*", "*Ch*")) +
  scale_x_continuous(breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4),
                     limits = c(-0.3, 0.4)) +
  scale_y_continuous(breaks = seq(from = -0.4, to = 0.3, by = 0.1),
                     limits = c(-0.4, 0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                                        margin = 
                                        margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12, 
                                        margin = 
                                        margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "inside",
        legend.position.inside = c(0.15, 0.9),
        legend.title = element_markdown(size = 10),
        legend.text = element_markdown(size = 8),
        legend.spacing.x = unit(0.25, "cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(nrow = 2,
                             theme = theme(legend.byrow = TRUE)))

# saving the graph

ggsave("PCA FA macroalgae.png", plot_pca_FA_macroalgae,
       height = 12.5, width = 17.5, units = "cm", dpi = 600)

## One-way ANOVAs --------------------------------------------------------------

# testing if the FA classes and ratios differ among the five macroalgal species

### SFA ------------------------------------------------------------------------

model_macroalgae_SFA <- 
  lm(log(SFA) ~ algae, data = dataset_macroalgae)

check_model(model_macroalgae_SFA, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_macroalgae_SFA)

anova(model_macroalgae_SFA)

pairs(emmeans(model_macroalgae_SFA, "algae"))

# defining the letters that represent the significant differences of the Tukey test

aov_macroalgae_SFA <- 
  aov(log(SFA) ~ algae, data = dataset_macroalgae)

tukey_macroalgae_SFA <- TukeyHSD(aov_macroalgae_SFA)

multcompLetters3("algae", "SFA", 
                 tukey_macroalgae_SFA$algae[,"p adj"], 
                 as.data.frame(dataset_macroalgae))

### MUFA ------------------------------------------------------------------------

model_macroalgae_MUFA <- 
  lm(log(MUFA) ~ algae, data = dataset_macroalgae)

check_model(model_macroalgae_MUFA, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_macroalgae_MUFA)

anova(model_macroalgae_MUFA)

pairs(emmeans(model_macroalgae_MUFA, "algae"))

# defining the letters that represent the significant differences of the Tukey test

aov_macroalgae_MUFA <- 
  aov(log(MUFA) ~ algae, data = dataset_macroalgae)

tukey_macroalgae_MUFA <- TukeyHSD(aov_macroalgae_MUFA)

multcompLetters3("algae", "MUFA", 
                 tukey_macroalgae_MUFA$algae[,"p adj"], 
                 as.data.frame(dataset_macroalgae))

### PUFA ------------------------------------------------------------------------

model_macroalgae_PUFA <- 
  lm(PUFA ~ algae, data = dataset_macroalgae)

check_model(model_macroalgae_PUFA, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_macroalgae_PUFA)

anova(model_macroalgae_PUFA)

pairs(emmeans(model_macroalgae_PUFA, "algae"))

# defining the letters that represent the significant differences of the Tukey test

aov_macroalgae_PUFA <- 
  aov(PUFA ~ algae, data = dataset_macroalgae)

tukey_macroalgae_PUFA <- TukeyHSD(aov_macroalgae_PUFA)

multcompLetters3("algae", "PUFA", 
                 tukey_macroalgae_PUFA$algae[,"p adj"], 
                 as.data.frame(dataset_macroalgae))

### n-3 / n-6 ------------------------------------------------------------------

model_macroalgae_n3_n6 <- 
  lm(`n-3/n-6` ~ algae, data = dataset_macroalgae)

check_model(model_macroalgae_n3_n6, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_macroalgae_n3_n6)

anova(model_macroalgae_n3_n6)

pairs(emmeans(model_macroalgae_n3_n6, "algae"))

dataset_macroalgae |> 
  group_by(algae) |> 
  summarise(mean_n_3_n_6 = mean(`n-3/n-6`),
            sd_n_3_n_6 = sd(`n-3/n-6`))

# defining the letters that represent the significant differences of the Tukey test

aov_macroalgae_n3_n6 <- 
  aov(`n-3/n-6` ~ algae, data = dataset_macroalgae)

tukey_macroalgae_n3_n6 <- TukeyHSD(aov_macroalgae_n3_n6)

multcompLetters3("algae", "n-3/n-6", 
                 tukey_macroalgae_n3_n6$algae[,"p adj"], 
                 as.data.frame(dataset_macroalgae))

### PUFA / SFA ------------------------------------------------------------------

model_macroalgae_PUFA_SFA <- 
  lm(log(`PUFA/SFA`) ~ algae, data = dataset_macroalgae)

check_model(model_macroalgae_PUFA_SFA, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_macroalgae_PUFA_SFA)

anova(model_macroalgae_PUFA_SFA)

### Calculating mean and sd for some FAs ---------------------------------------

dataset_macroalgae |> 
  group_by(algae) |> 
  summarise(mean_14 = mean(`14:0`),
            sd_14 = sd(`14:0`),
            mean_16 = mean(`16:0`),
            sd_16 = sd(`16:0`),
            mean_18 = mean(`18:0`),
            sd_18 = sd(`18:0`),
            mean_18_1_n_9 = mean(`18:1 n-9`),
            sd_18_1_n_9 = sd(`18:1 n-9`))

## Graphics for macroalgae -----------------------------------------------------

# selecting only the columns of the FA classes and transforming them into a new factor

dataset_macroalgae_classes <-
  dataset_macroalgae |> 
  select(algae, replicate, SFA, MUFA, PUFA) |> 
  pivot_longer(cols = c(SFA, MUFA, PUFA),
               names_to = "class",
               values_to = "abundance") |> 
  mutate(class = factor(class,
                        levels = c("SFA", "MUFA", "PUFA")))

# calculating mean and standard errors

summary_macroalgae_FA_classes <-
  dataset_macroalgae_classes |> 
  group_by(class, algae) |> 
  summarise(mean_abundance = mean(abundance),
            sd_abundance = sd(abundance),
            sample_size = n()) |> 
  mutate(se_abundance = sd_abundance / sqrt(sample_size))

# plotting

plot_macroalgae_FA_classes_a <-
  ggplot(summary_macroalgae_FA_classes,
         aes(x = algae, y = mean_abundance, fill = algae)) +
  geom_bar(stat = "identity", position = position_dodge(),
           colour = "black", width = 0.5) + 
  geom_errorbar(aes(ymin = mean_abundance,
                    ymax = mean_abundance + se_abundance),
                width = 0.2,
                position = position_dodge(0.9)) +
  facet_wrap(~ class,
             ncol = 1) +
  labs(x = "**Species**", 
       y = expression(bold("\u03bc"*g~"FA"~"/"~"mg of DW")), 
       fill = "**Species**")

plot_macroalgae_FA_classes_b <-
  plot_macroalgae_FA_classes_a +
  scale_x_discrete(labels = c("*Ulva*", "*Fucus*", "*Laminaria*", 
                              "*Gracilaria*", "*Chondrus*")) +
  scale_y_continuous(breaks = seq(0, 4, by = 1),
                     limits = c(0, 4),
                     expand = c(0, 0)) + 
  scale_fill_manual(values = c("#70AD47", "#FFCC00", "#993300",
                               "#FF6666", "#FF0000"))

plot_macroalgae_FA_classes_c <-
  plot_macroalgae_FA_classes_b +
  theme_bw() +
  theme(axis.text.x = element_markdown(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                                        margin = 
                                        margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_text(size = 12, 
                                    margin = 
                                    margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "none",
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.background = element_rect(colour = "black", linewidth = 0.5),
        panel.spacing.y = unit(0.25, "cm"))

ggsave("Macroalgae FA classes.png",
       plot_macroalgae_FA_classes_c,
       height = 22.5, width = 15, units = "cm", dpi = 600)

# 2 - Amphipods ---------------------------------------------------------------

## PERMANOVA -------------------------------------------------------------------

# testing if the FA profiles of neonates significantly differed according to the different algal diets of parental amphipods

# separating the numerical data and the factors in different matrices

data_fa_amphipods <-
  dataset_amphipods |>
  dplyr::select(`16:0` : `22:6 n-3`)

data_fa_amphipods_log <-
  dataset_amphipods |>
  dplyr::select(`16:0` : `22:6 n-3`) |> 
  mutate(across(`16:0` : `22:6 n-3`, ~ log(.x + 1)))

factor_amphipods <-
  dataset_amphipods |> 
  dplyr::select(algae : replicate)

# running the test

with(factor_amphipods, 
     adonis2(data_fa_amphipods_log ~ algae,
             method = "euc", by = "terms",
             data = factor_amphipods))

# pairwise comparison tests

pairwise.adonis2(
  vegdist(data_fa_amphipods_log, 
          method = "euclidean") ~ algae,
  data = as.data.frame(dataset_amphipods))

# SIMPER analysis 

simper_amphipods <- 
  with(factor_amphipods, 
       simper(data_fa_amphipods_log, algae))

summary(simper_amphipods) 

## PCoA ------------------------------------------------------------------------

# creating a distance matrix

distance_fa_amphipods <-
  vegdist(data_fa_amphipods_log, method = "euclidean")

# running the PCoA

pcoa_fa_amphipods <- 
  pcoa(distance_fa_amphipods)

pcoa_fa_amphipods$values

# extracting the PCoA vectors

pcoa_fa_amphipods_arrows <- 
  compute_arrows(pcoa_fa_amphipods, data_fa_amphipods_log)

pcoa_fa_amphipods_arrows$vectors <- 
  as.data.frame(pcoa_fa_amphipods_arrows$vectors)

dataset_amphipods$axis_1 <- pcoa_fa_amphipods_arrows$vectors[,1]

dataset_amphipods$axis_2 <- pcoa_fa_amphipods_arrows$vectors[,2]

# scaling the vectors for a better visualization

arrows_df_amphipods <- as.data.frame(pcoa_fa_amphipods_arrows$U * 1)

arrows_df_amphipods$variable <- rownames(arrows_df_amphipods)

# plotting the PCoA

plot_pcoa_amphipods_a <- 
  ggplot() +
  geom_point(data = dataset_amphipods, 
             aes(axis_1, axis_2, fill = algae), 
             size = 4, shape = 21) +
  labs(x = "**PCO1 (56.5 % of total variation)**",
       y = "**PCO2 (27.5 % of total variation)**",
       fill = "**Diet**")

plot_pcoa_amphipods_b <-
  plot_pcoa_amphipods_a +
  scale_x_continuous(breaks = seq(-1, 1.5, by = 0.5),
                     limits = c(-1, 1.5)) +
  scale_y_continuous(breaks = seq(-0.75, 0.75, by = 0.25),
                     limits = c(-0.75, 0.75)) +
  scale_fill_manual(values = c("#70AD47", "#FFCC00", "#993300",
                               "#FF6666", "#FF0000", "#4472C4"),
                    labels = c("*Ulva*", "*Fucus*", "*Laminaria*", 
                               "*Gracilaria*", "*Chondrus*", "Wild")) 

plot_pcoa_amphipods_c <-
  plot_pcoa_amphipods_b +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                                        margin = 
                                        margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12, 
                                        margin = 
                                        margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "inside",
        legend.position.inside = c(0.22, 0.9),
        legend.title = element_markdown(size = 10),
        legend.text = element_markdown(size = 8),
        legend.spacing.x = unit(0.25, "cm"),
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot_pcoa_amphipods_d <-
  plot_pcoa_amphipods_c +
  guides(fill = guide_legend(nrow = 2),
                theme = theme(legend.byrow = TRUE)) + 
  geom_segment(data = as.data.frame(pcoa_fa_amphipods_arrows$U * 1),
               x = 0, y = 0, 
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(2, "mm"))) +
  geom_text(data = arrows_df_amphipods, 
            aes(x = Axis.1, y = Axis.2, label = variable),
            position = position_jitter(), 
            size = 4)

# saving the graph

ggsave("PCoA FA amphipods.png", plot_pcoa_amphipods_d,
       height = 12.5, width = 17.5, units = "cm", dpi = 600)

## PCA -------------------------------------------------------------------------

pca_FA_amphipods <-
  prcomp(data_fa_amphipods_log,
         center=TRUE,scale.=TRUE)

summary(pca_FA_amphipods)

# extract PC axes

pca_FA_amphipods_values <-
  data.frame(dataset_amphipods$algae,
             dataset_amphipods$replicate,
             pca_FA_amphipods$x)

colnames(pca_FA_amphipods_values)[1] <-"algae"

colnames(pca_FA_amphipods_values)[2] <-"replicate"

pca_FA_amphipods_values <-
  pca_FA_amphipods_values  |> 
  as_tibble()

plot_pca_FA_amphipods <-
  autoplot(pca_FA_amphipods, data = dataset_amphipods, 
           fill = 'algae', shape = 21, size = 4,
           loadings = TRUE, loadings.colour = 'black',
           loadings.label = TRUE,
           loadings.label.size = 3,
           loadings.label.colour = "black") +
  labs(x = "**PC1 (58.04 % of total variation)**",
       y = "**PC2 (21.67 % of total variation)**",
       fill = "**Diet**") +
  scale_fill_manual(values = c("#70AD47", "#FFCC00", "#993300",
                               "#FF6666", "#FF0000", "#4472C4"),
                    labels = c("*Ul*", "*Fu*", "*La*", "*Gr*", "*Ch*", "Wd")) +
  scale_x_continuous(breaks = seq(from = -0.4, to = 0.3, by = 0.1),
                     limits = c(-0.3, 0.4)) +
  scale_y_continuous(breaks = seq(from = -0.4, to = 0.6, by = 0.2),
                     limits = c(-0.4, 0.6)) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                                        margin = 
                                        margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12, 
                                        margin = 
                                        margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "inside",
        legend.position.inside = c(0.15, 0.9),
        legend.title = element_markdown(size = 10),
        legend.text = element_markdown(size = 8),
        legend.spacing.x = unit(0.25, "cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(nrow = 2,
                             theme = theme(legend.byrow = TRUE)))

# saving the graph

ggsave("PCA FA amphipods.png", plot_pca_FA_amphipods,
       height = 12.5, width = 17.5, units = "cm", dpi = 600)

## One-way ANOVAs ---------------------------------------------------------------

# testing if the dry weight and the number of neonates per brood differ according to the parental diet 

### Dry weight -----------------------------------------------------------------

model_amphipods_dry_weight <- 
  lm(log(dry_weight) ~ algae, data = dataset_amphipods)

check_model(model_amphipods_dry_weight, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_dry_weight)

anova(model_amphipods_dry_weight)

pairs(emmeans(model_amphipods_dry_weight, "algae"))

# defining the letters to put over the barplot to show significant differences

aov_dry_weight <- 
  aov(log(dry_weight) ~ algae, data = dataset_amphipods)

tukey_dry_weight <- TukeyHSD(aov_dry_weight)

multcompLetters3("algae", "dry_weight", 
                 tukey_dry_weight$algae[,"p adj"], 
                 as.data.frame(dataset_amphipods))

dataset_amphipods |> 
  group_by(algae) |> 
  summarise(mean_dw = mean(dry_weight),
            sd_dw = sd(dry_weight),
            sample_size = n()) |> 
  mutate(se_dw = sd_dw / sqrt(sample_size)) |> 
  as.data.frame()

### Number of newborns per brood -----------------------------------------------

model_amphipods_n_brood <- 
  lm(n_brood ~ algae, data = dataset_amphipods)

check_model(model_amphipods_n_brood, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_n_brood)

anova(model_amphipods_n_brood)

pairs(emmeans(model_amphipods_n_brood, "algae"))

# defining the letters to put over the barplot to show significant differences

aov_n_brood <- 
  aov(n_brood ~ algae, data = dataset_amphipods)

tukey_n_brood <- TukeyHSD(aov_n_brood)

multcompLetters3("algae", "n_brood", 
                 tukey_n_brood$algae[,"p adj"], 
                 as.data.frame(dataset_amphipods))

dataset_amphipods |> 
  group_by(algae) |> 
  summarise(mean_brood = mean(n_brood),
            sd_brood = sd(n_brood),
            sample_size = n()) |> 
  mutate(se_brood = sd_brood / sqrt(sample_size)) |> 
  as.data.frame()

### 16:0 -----------------------------------------------------------------------

model_amphipods_16_0 <- 
  lm(log(`16:0`) ~ algae, data = dataset_amphipods)

check_model(model_amphipods_16_0, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_16_0)

anova(model_amphipods_16_0)

dataset_amphipods |> 
  group_by(algae) |> 
  summarise(mean_16 = mean(`16:0`),
            sd_16 = sd(`16:0`))

### 18:0 -----------------------------------------------------------------------

model_amphipods_18_0 <- 
  lm(log(`18:0`) ~ algae, data = dataset_amphipods)

check_model(model_amphipods_18_0, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_18_0)

anova(model_amphipods_18_0)

pairs(emmeans(model_amphipods_18_0, "algae"))

dataset_amphipods |> 
  group_by(algae) |> 
  summarise(mean_18 = mean(`18:0`),
            sd_18 = sd(`18:0`))

# defining the letters to put over the barplot to show significant differences

aov_18 <- 
  aov(log(`18:0`) ~ algae, data = dataset_amphipods)

tukey_18 <- TukeyHSD(aov_18)

multcompLetters3("algae", "18:0", 
                 tukey_18$algae[,"p adj"], 
                 as.data.frame(dataset_amphipods))

### 18:1 n-9 -------------------------------------------------------------------

model_amphipods_18_1_n_9 <- 
  lm(`18:1 n-9` ~ algae, data = dataset_amphipods)

check_model(model_amphipods_18_1_n_9, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_18_1_n_9)

anova(model_amphipods_18_1_n_9)

dataset_amphipods |> 
  group_by(algae) |> 
  summarise(mean_18_1_n_9 = mean(`18:1 n-9`),
            sd_18_1_n_9 = sd(`18:1 n-9`))

### 18:1 n-7 -------------------------------------------------------------------

model_amphipods_18_1_n_7 <- 
  lm(`18:1 n-7` ~ algae, data = dataset_amphipods)

check_model(model_amphipods_18_1_n_7, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_18_1_n_7)

anova(model_amphipods_18_1_n_7)

### 18:2 n-6 -------------------------------------------------------------------

model_amphipods_18_2_n_6 <- 
  lm(`18:2 n-6` ~ algae, data = dataset_amphipods)

check_model(model_amphipods_18_2_n_6, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_18_2_n_6)

anova(model_amphipods_18_2_n_6)

### 18:3 n-3 -------------------------------------------------------------------

model_amphipods_18_3_n_3 <- 
  lm(`18:3 n-3` ~ algae, data = dataset_amphipods)

check_model(model_amphipods_18_3_n_3, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_18_3_n_3)

anova(model_amphipods_18_3_n_3)

### 20:4 n-6 -------------------------------------------------------------------

model_amphipods_20_4_n_6 <- 
  lm(log(`20:4 n-6`) ~ algae, data = dataset_amphipods)

check_model(model_amphipods_20_4_n_6, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_20_4_n_6)

anova(model_amphipods_20_4_n_6)

### 20:5 n-3 (EPA) -------------------------------------------------------------

model_amphipods_20_5_n_3 <- 
  lm(`20:5 n-3` ~ algae, data = dataset_amphipods)

check_model(model_amphipods_20_5_n_3, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_20_5_n_3)

anova(model_amphipods_20_5_n_3)

dataset_amphipods |> 
  group_by(algae) |> 
  summarise(mean_20_5_n_3 = mean(`20:5 n-3`),
            sd_20_5_n_3 = sd(`20:5 n-3`))

### 22:5 n-3 -------------------------------------------------------------------

model_amphipods_22_5_n_3 <- 
  lm(`22:5 n-3` ~ algae, data = dataset_amphipods)

check_model(model_amphipods_22_5_n_3, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_22_5_n_3)

anova(model_amphipods_22_5_n_3)

### 22:6 n-3 (DHA) -------------------------------------------------------------

model_amphipods_22_6_n_3 <- 
  lm(`22:6 n-3` ~ algae, data = dataset_amphipods)

check_model(model_amphipods_22_6_n_3, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_22_6_n_3)

anova(model_amphipods_22_6_n_3)

pairs(emmeans(model_amphipods_22_6_n_3, "algae"))

dataset_amphipods |> 
  group_by(algae) |> 
  summarise(mean_22_6_n_3 = mean(`22:6 n-3`),
            sd_22_6_n_3 = sd(`22:6 n-3`))

# defining the letters to put over the barplot to show significant differences

aov_22_6_n_3 <- 
  aov(`22:6 n-3` ~ algae, 
      data = dataset_amphipods)

tukey_22_6_n_3 <- TukeyHSD(aov_22_6_n_3)

multcompLetters3("algae", "22:6 n-3", 
                 tukey_22_6_n_3$algae[,"p adj"], 
                 as.data.frame(dataset_amphipods), reversed = TRUE)

### SFA ------------------------------------------------------------------------

model_amphipods_SFA <- 
  lm(log(SFA) ~ algae, data = dataset_amphipods)

check_model(model_amphipods_SFA, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_SFA)

anova(model_amphipods_SFA)

pairs(emmeans(model_amphipods_SFA, "algae"))

dataset_amphipods |> 
  group_by(algae) |> 
  summarise(mean_SFA = mean(SFA),
            sd_SFA = sd(SFA))

# defining the letters to put over the barplot to show significant differences

aov_SFA <- 
  aov(log(SFA) ~ algae, 
      data = dataset_amphipods)

tukey_SFA <- TukeyHSD(aov_SFA)

multcompLetters3("algae", "SFA", 
                 tukey_SFA$algae[,"p adj"], 
                 as.data.frame(dataset_amphipods), reversed = TRUE)

### MUFA ------------------------------------------------------------------------

model_amphipods_MUFA <- 
  lm(MUFA ~ algae, data = dataset_amphipods)

check_model(model_amphipods_MUFA, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_MUFA)

anova(model_amphipods_MUFA)

dataset_amphipods |> 
  group_by(algae) |> 
  summarise(mean_MUFA = mean(MUFA),
            sd_MUFA = sd(MUFA))

### PUFA ------------------------------------------------------------------------

model_amphipods_PUFA <- 
  lm(PUFA ~ algae, data = dataset_amphipods)

check_model(model_amphipods_PUFA, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_PUFA)

anova(model_amphipods_PUFA)

dataset_amphipods |> 
  group_by(algae) |> 
  summarise(mean_PUFA = mean(PUFA),
            sd_PUFA = sd(PUFA))

### n-3 / n-6 ------------------------------------------------------------------

model_amphipods_n3_n6 <- 
  lm(`n-3/n-6` ~ algae, data = dataset_amphipods)

check_model(model_amphipods_n3_n6, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_n3_n6)

anova(model_amphipods_n3_n6)

dataset_amphipods |> 
  group_by(algae) |> 
  summarise(mean_n3_n6 = mean(`n-3/n-6`),
            sd_n3_n6 = sd(`n-3/n-6`))

### EPA / DHA ------------------------------------------------------------------

model_amphipods_EPA_DHA <- 
  lm(`EPA/DHA` ~ algae, data = dataset_amphipods)

check_model(model_amphipods_EPA_DHA, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_EPA_DHA)

anova(model_amphipods_EPA_DHA)

pairs(emmeans(model_amphipods_EPA_DHA, "algae"))

dataset_amphipods |> 
  group_by(algae) |> 
  summarise(mean_EPA_DHA = mean(`EPA/DHA`),
            sd_EPA_DHA = sd(`EPA/DHA`))

# defining the letters to put over the barplot to show significant differences

aov_EPA_DHA <- 
  aov(`EPA/DHA` ~ algae, 
      data = dataset_amphipods)

tukey_EPA_DHA <- TukeyHSD(aov_EPA_DHA)

multcompLetters3("algae", "EPA/DHA", 
                 tukey_EPA_DHA$algae[,"p adj"], 
                 as.data.frame(dataset_amphipods))

### PUFA / SFA ------------------------------------------------------------------

model_amphipods_PUFA_SFA <- 
  lm(`PUFA/SFA` ~ algae, data = dataset_amphipods)

check_model(model_amphipods_PUFA_SFA, 
            check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_amphipods_PUFA_SFA)

anova(model_amphipods_PUFA_SFA)

pairs(emmeans(model_amphipods_PUFA_SFA, "algae"))

dataset_amphipods |> 
  group_by(algae) |> 
  summarise(mean_PUFA_SFA = mean(`PUFA/SFA`),
            sd_PUFA_SFA = sd(`PUFA/SFA`))

# defining the letters to put over the barplot to show significant differences

aov_PUFA_SFA <- 
  aov(`PUFA/SFA` ~ algae, 
      data = dataset_amphipods)

tukey_PUFA_SFA <- TukeyHSD(aov_PUFA_SFA)

multcompLetters3("algae", "PUFA/SFA", 
                 tukey_PUFA_SFA$algae[,"p adj"], 
                 as.data.frame(dataset_amphipods))

## Graphics for neonates -------------------------------------------------------

### Dry weight + Number per brood ----------------------------------------------

# calculating mean and standard errors

summary_neonates_traits <-
  dataset_amphipods |> 
  group_by(algae) |> 
  summarise(across(c(dry_weight, n_brood),
                   list(mean = mean,
                        sd = sd))) |> 
  mutate(dry_weight_se = dry_weight_sd / sqrt(5),
         n_brood_se = n_brood_sd / sqrt(5))

# adding letters representing significant differences in the table  

df_letters_n_brood <-
  as.data.frame.list(letters_n_brood$algae) |> 
  rownames_to_column(var = "algae") |> 
  as_tibble() |> 
  select(algae, Letters) |> 
  mutate(algae = factor(algae,
                        levels = c("ulva", "fucus", "laminaria", 
                                   "gracilaria", "chondrus", "wild"))) |> 
  arrange(algae) |> 
  select(Letters) |> 
  pull()

df_letters_dry_weight <-
  as.data.frame.list(letters_dry_weight$algae) |> 
  rownames_to_column(var = "algae") |> 
  as_tibble() |> 
  select(algae, Letters) |> 
  mutate(algae = factor(algae,
                        levels = c("ulva", "fucus", "laminaria", 
                                   "gracilaria", "chondrus", "wild"))) |> 
  arrange(algae) |> 
  select(Letters) |> 
  pull()

# number of newborns per brood

plot_amphipoda_n_brood_a <-
  ggplot(summary_neonates_traits,
         aes(x = algae, y = n_brood_mean, fill = algae)) +
  geom_bar(stat = "identity", position = position_dodge(),
           colour = "black", width = 0.5) + 
  geom_errorbar(aes(ymin = n_brood_mean,
                    ymax = n_brood_mean + n_brood_se),
                width = 0.2,
                position = position_dodge(0.9)) +
  labs(x = "**Diet**", y = "**# Offspring / brood**", fill = "**Diet**")

plot_amphipoda_n_brood_b <-
  plot_amphipoda_n_brood_a +
  scale_x_discrete(labels = c("*Ulva*", "*Fucus*", "*Laminaria*", 
                              "*Gracilaria*", "*Chondrus*", "Wild")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20),
                     limits = c(0, 100),
                     expand = c(0, 0)) + 
  scale_fill_manual(values = c("#70AD47", "#FFCC00", "#993300",
                               "#FF6666", "#FF0000", "#4472C4"))

plot_amphipoda_n_brood_c <-
  plot_amphipoda_n_brood_b +
  theme_bw() +
  theme(axis.text.x = element_markdown(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                                        margin = 
                                        margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12, 
                                        margin = 
                                        margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "none",
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# dry weight plot

plot_amphipoda_dw_a <-
  ggplot(summary_neonates_traits,
       aes(x = algae, y = dry_weight_mean, fill = algae)) +
  geom_bar(stat = "identity", position = position_dodge(),
           colour = "black", width = 0.5) + 
  geom_errorbar(aes(ymin = dry_weight_mean,
                    ymax = dry_weight_mean + dry_weight_se),
                width = 0.2,
                position = position_dodge(0.9)) +
  labs(x = "**Diet**", y = "**Dry Weight (mg)**", fill = "**Diet**")

plot_amphipoda_dw_b <-
  plot_amphipoda_dw_a +
  scale_x_discrete(labels = c("*Ulva*", "*Fucus*", "*Laminaria*", 
                              "*Gracilaria*", "*Chondrus*", "Wild")) +
  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5),
                     limits = c(0, 2.5),
                     expand = c(0, 0)) + 
  scale_fill_manual(values = c("#70AD47", "#FFCC00", "#993300",
                             "#FF6666", "#FF0000", "#4472C4"))

plot_amphipoda_dw_c <-
  plot_amphipoda_dw_b +
  theme_bw() +
  theme(axis.text.x = element_markdown(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                                        margin = 
                                        margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12, 
                                        margin = 
                                        margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "none",
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# joining the plots

plot_amphipoda_traits <-
  plot_amphipoda_n_brood_c +
  plot_amphipoda_dw_c +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = 'bold'))

# remove x title and text for plot 1

plot_amphipoda_traits[[1]] <-
  plot_amphipoda_traits[[1]] +
  theme(
    plot.margin = margin(0.5, 0.5, 0.25, 0.5, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank())

# change margins for plot 2

plot_amphipoda_traits[[2]] <-
  plot_amphipoda_traits[[2]] +
  theme(
    plot.margin = margin(0.25, 0.5, 0.5, 0.5, "cm"))

ggsave("Neonate traits.png",
       plot_amphipoda_traits,
       height = 20, width = 15, units = "cm", dpi = 600)

### FA classes -----------------------------------------------------------------

# selecting only the columns of the FA classes and transforming them into a new factor

dataset_amphipods_classes <-
  dataset_amphipods |> 
  select(algae, replicate, SFA, MUFA, PUFA) |> 
  pivot_longer(cols = c(SFA, MUFA, PUFA),
               names_to = "class",
               values_to = "abundance") |> 
  mutate(class = factor(class,
                        levels = c("SFA", "MUFA", "PUFA")))

# calculating mean and standard errors

summary_amphipods_FA_classes <-
  dataset_amphipods_classes |> 
  group_by(class, algae) |> 
  summarise(mean_abundance = mean(abundance),
            sd_abundance = sd(abundance),
            sample_size = n()) |> 
  mutate(se_abundance = sd_abundance / sqrt(sample_size))

# plotting

plot_amphipoda_FA_classes_a <-
  ggplot(summary_amphipods_FA_classes,
         aes(x = algae, y = mean_abundance, fill = algae)) +
  geom_bar(stat = "identity", position = position_dodge(),
           colour = "black", width = 0.5) + 
  geom_errorbar(aes(ymin = mean_abundance,
                    ymax = mean_abundance + se_abundance),
                width = 0.2,
                position = position_dodge(0.9)) +
  facet_wrap(~ class,
             ncol = 1) +
  labs(x = "**Diet**", 
       y = expression(bold("\u03bc"*g~"FA"~"/"~"mg of DW")), 
       fill = "**Diet**")

plot_amphipoda_FA_classes_b <-
  plot_amphipoda_FA_classes_a +
  scale_x_discrete(labels = c("*Ulva*", "*Fucus*", "*Laminaria*", 
                              "*Gracilaria*", "*Chondrus*", "Wild")) +
  scale_y_continuous(breaks = seq(0, 8, by = 2),
                     limits = c(0, 8),
                     expand = c(0, 0)) + 
  scale_fill_manual(values = c("#70AD47", "#FFCC00", "#993300",
                               "#FF6666", "#FF0000", "#4472C4"))

plot_amphipoda_FA_classes_c <-
  plot_amphipoda_FA_classes_b +
  theme_bw() +
  theme(axis.text.x = element_markdown(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                                        margin = 
                                        margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_text(size = 12, 
                                        margin = 
                                        margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "none",
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.background = element_rect(colour = "black", linewidth = 0.5),
        panel.spacing.y = unit(0.25, "cm"))

ggsave("Neonate FA classes.png",
       plot_amphipoda_FA_classes_c,
       height = 22.5, width = 15, units = "cm", dpi = 600)

### FA differences to wild "control" -------------------------------------------

# selecting only the columns of the most common FAs and transforming them into a new factor

dataset_amphipods_FA <-
  dataset_amphipods |> 
  select(algae, replicate, `16:0`:`22:6 n-3`) |> 
  pivot_longer(cols = `16:0`:`22:6 n-3`,
               names_to = "FAME",
               values_to = "abundance") |> 
  mutate(FAME = factor(FAME))

# calculating mean 

summary_amphipods_FA <-
  dataset_amphipods_FA |> 
  group_by(FAME, algae) |> 
  summarise(mean_abundance = mean(abundance))

# calculating the absolute difference between the mean FAs of the macroalgal diets to the wild baseline 

difference_amphipods_FA_wild <-
  summary_amphipods_FA |> 
  mutate(across(mean_abundance, ~.x-.x[algae == "wild"])) |> 
  rename(dif_wild = mean_abundance)

# plotting

plot_amphipoda_FA_differences_a <-
  ggplot(difference_amphipods_FA_wild |> 
         filter(algae != "wild"),
         aes(x = FAME, y = dif_wild, fill = algae)) +
  geom_bar(stat = "identity", position = position_dodge(0.75),
           colour = "black", width = 0.5) +
  labs(x = "**Fatty Acids**", 
       y = expression(
           bold("Difference in FA abundance "*"(\u03bc"*g~"FA"~"/"~"mg of DW)")), 
       fill = "**Diet**")

plot_amphipoda_FA_differences_b <-
  plot_amphipoda_FA_differences_a +
  scale_x_discrete(labels = c("16:0", "18:0", "18:1*n*-7", "18:1*n*-9",
                              "18:2*n*-6", "18:3*n*-3", "20:4*n*-6",
                              "20:5*n*-3", "22:5*n*-3", "22:5*n*-6")) +
  scale_y_continuous(breaks = seq(-2.5, 0.5, by = 0.5),
                     limits = c(-2.5, 0.5),
                     expand = c(0, 0)) + 
  scale_fill_manual(values = c("#70AD47", "#FFCC00", "#993300",
                               "#FF6666", "#FF0000"),
                    labels = c("*Ulva*", "*Fucus*", "*Laminaria*", 
                               "*Gracilaria*", "*Chondrus*")) +
  geom_hline(yintercept = 0,  linetype = "dashed", linewidth = 1,
             colour = "#4472C4")

plot_amphipoda_FA_differences_c <-
  plot_amphipoda_FA_differences_b +
  theme_bw() +
  theme(axis.text.x = element_markdown(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                                        margin = 
                                        margin(t = 20, r = 20, b = 10, l = 20)),
        axis.title.y = element_text(size = 12, 
                                    margin = 
                                    margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.15),
        legend.title = element_markdown(size = 10),
        legend.text = element_markdown(size = 10),
        legend.spacing.x = unit(0.5, "cm"),
        legend.key.spacing.y = unit(0.5, "cm"),
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(nrow = 2,
                             theme = theme(legend.byrow = TRUE)))

ggsave("Neonate FA differences.png",
       plot_amphipoda_FA_differences_c,
       height = 17.5, width = 22.5, units = "cm", dpi = 600)
