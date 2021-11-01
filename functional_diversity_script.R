# 1: DATA PREPARATION ####


# R version used: 4.0.3.

# Load the needed libraries.
library(readxl)
library(glmmTMB) # Version used: 1.1.2.2
library(tidyverse) # 1.3.1
library(DHARMa) # Version used: 0.4.3
library(MuMIn) # Version used: 1.43.17
library(Information) # Version used: 0.0.9
library(Hmisc) # Version used: 4.4-2
library(corrplot) # Version used: 0.84
library(ltm) # Version used: 1.1-1 
library(FD) # Version used: 1.0-12
library(xlsx) # Version used: 0.6.5
library(coRanking) # Version used: 0.2.1
library(betapart) # Version used: 1.5.4
library(BAT) # Version used: 2.7.0

# Read in the source files, which are the community matrices and the trait
# matrix.
read.csv2(file = "functional diversity/data/1998_breeding_all_species.csv", 
          header = T) -> all_species_1998
read.csv2(file = "functional diversity/data/2018_breeding_all_species.csv", 
          header = T) -> all_species_2018

# Ready the data frames. These data frames have all the species wetland species 
# that have ever been detected in the study area. Therefore, we need to filter 
# out those that do not occur in our dataset.
all_species_1998[,1] -> rownames(all_species_1998)
all_species_2018[,1] -> rownames(all_species_2018)
all_species_1998[,1] <- NULL ; all_species_2018[,1] <- NULL
all_species_1998[rowSums(all_species_1998)!=0, ] -> study_comm_1998
all_species_2018[rowSums(all_species_2018)!=0, ] -> study_comm_2018

# Read in the trait matrix.
read_excel(path = "functional diversity/data/functional_trait_matrix.xlsx", 
           sheet = 1) -> trait_matrix_raw
class(trait_matrix_raw)
glimpse(trait_matrix_raw)


# 2: FUNCTIONAL TRAITS ####


## 2.1: BRAIN MASS RESIDUALS ####

# Create a new column by using the residuals from the regression brain 
# mass ~ body mass. The data in the "brain_mass_residuals_va" is raw brain mass
# data (in grams). Build a generalized linear model (GLM) and obtain the 
# residuals.
trait_matrix_raw %>% 
  filter(!is.na(brain_mass_residuals_va)) %>%
  dplyr::summarise(min = min(brain_mass_residuals_va), 
                   max = max(brain_mass_residuals_va))
# We have a positive continuous response variable, so we will use a GLM with 
# gamma distributed errors.
m_brain_mass <- glmmTMB(data = trait_matrix_raw, 
                       formula = brain_mass_residuals_va ~ body_mass_wi_du,
                       family = Gamma(link = "log"))
summary(m_brain_mass) # The effect is significant. Let's check the model.

# Model validation #
# 1) Check the distribution of the residuals.
boxplot(residuals(m_brain_mass)); abline(h = 0, col = "red", lty = 2)
# There is an outlier.
boxplot(residuals(m_brain_mass), ylim = c(-10,10))
abline(h = 0, col = "red", lty = 2) # The distribution looks okay.
# 2) Plot the residuals against the only predictor.
plot(trait_matrix_raw$
         brain_mass_residuals_va
       [!is.na(trait_matrix_raw$brain_mass_residuals_va)],
     residuals(m_brain_mass), ylim = c(-10,10))
abline(h = 0, col = "red", lty = 2) # There is a pattern in the plot. We may try
# to log-transform the body mass and re-model the data.
log(trait_matrix_raw$body_mass_wi_du) -> body_mass_log_trans
m_brain_mass <- glm(formula = trait_matrix_raw$
                            brain_mass_residuals_va ~ body_mass_log_trans,
                       family = Gamma(link = "log"))
summary(m_brain_mass) # The effect persists. 
# Check the residuals.
plot(trait_matrix_raw$
         brain_mass_residuals_va
       [!is.na(trait_matrix_raw$brain_mass_residuals_va)],
     residuals(m_brain_mass), ylim = c(-10,10))
abline(h = 0, col = "red", lty = 2)
# The residuals look much better now. Proceed with this model.
# 1) Check the distribution of the residuals.
boxplot(residuals(m_brain_mass)); abline(h = 0, col = "red", lty = 2)
# There is an outlier.
boxplot(residuals(m_brain_mass), ylim = c(-10,10))
abline(h = 0, col = "red", lty = 2) # The distribution looks okay.
# 2) The second step has already been conducted. See above.
# 3) Check the model for autocorrelation.
par(mfrow = c(1, 2))
acf(residuals(m_brain_mass)); pacf(residuals(m_brain_mass)) # No problem here.
par(mfrow = c(1, 1))
# 4) Do the DHARMa tests by simulating residuals from the model.
set.seed(13)
simulateResiduals(m_brain_mass, n = 1000) -> sim_brain_model
set.seed(13); testResiduals(sim_brain_model)
# No problem here. We can assume that the model is a valid one.
# We can now use the residuals from this as a trait.
residuals(m_brain_mass, type = "response") -> resid_brain

# Calculate an R-squared value for this model to report in the model.
r.squaredGLMM(m_brain_mass) # r2 = 0.94. 

# Deal with the NAs and include the residuals into the data frame.
trait_matrix_raw %>% 
  filter(!is.na(brain_mass_residuals_va)) %>% 
  mutate(brain_mass_residuals_va = resid_brain) -> resid_matrix
trait_matrix_raw %>% 
  filter(is.na(brain_mass_residuals_va)) %>% 
  bind_rows(resid_matrix) %>% 
  arrange(species) -> trait_matrix # This is the new trait matrix that has the
# the brain mass residuals instead of raw brain mass data.

# Now, check the data classes of each trait.
glimpse(trait_matrix)
# The categorical trait, early developmental mode, should be converted to 
# factor.
trait_matrix %>% 
  mutate(early_developmental_mode_st = 
           as.factor(early_developmental_mode_st )) -> trait_matrix
glimpse(trait_matrix)


## 2.2: REMOVING CORRELATED TRAITS ####


# Firstly, there are some numeric traits (e.g., diet traits) that need to be 
# binarized.
func_binary <- function(x){
  Information::is.binary(as.matrix(trait_matrix[, x]))
} # To get the binaries.
sapply(X = (1:88), func_binary) -> indexer
trait_matrix[, indexer] -> binaries

# Remember, there is only one categorical variable and it's the early 
# developmental mode. 
trait_matrix[, !indexer] -> non_binaries
non_binaries[, -c(1,12)] -> non_binaries # Removing the categorical one and the
# species name column.
# There are some binaries in this tibble, and that's because they have NAs. 
# Remove them manually and add them to the binaries tibble.
non_binaries[, 28:39]
bind_cols(non_binaries[, 28:39], binaries) -> binaries
non_binaries[, -c(28:39)] -> non_binaries

# Now that we are ready, we can start by checking the correlation between the 
# continuous variables.
corr_matrix_cont <- rcorr(as.matrix(non_binaries))

# Plot the correlation matrix to check the values manually and to select between
# the correlated traits. 
png(filename = "functional diversity/figs/corr_plot_cont.png", width = 10, 
    height = 10, units = "in", res = 800, type = "cairo")
corrplot(corr = corr_matrix_cont$r, type = "upper", diag = F, outline = T, 
         tl.col = "black", method = "number", p.mat = corr_matrix_cont$P, 
         tl.cex = 0.75, insig = "blank", number.cex = 0.4, cl.cex = 1)
dev.off()
# Remove one from each "correlation couple". While removing, pick the one that
# is significantly correlated with more of the other traits. We set the
# threshold correlation coefficient (Spearman's r) level as 0.70 (inclusive).

# Here, we chose to remove body length, bill length, wing length, incubation 
# duration, broods per year, and generation length (a total of six traits).
non_binaries %>% 
  dplyr::select(-body_length_st, -bill_length_st, -wing_length_st, 
         -incubation_duration_st, -generation_length_bl,
         -broods_per_year_st) -> non_binaries_uncorr
# Next, check the correlation between the binaries by using Spearman's Rho.
corr_matrix_bins <- rcorr(as.matrix(binaries), type = "spearman")
# Plot the matrix and pick out the correlated ones.
 png(filename = "functional diversity/figs/corr_plot_bins.png", width = 10, 
     height = 10, units = "in", res = 800, type = "cairo")
 corrplot(corr = corr_matrix_bins$r, type = "upper", diag = F, outline = T, 
         tl.col = "black", method = "number", p.mat = corr_matrix_bins$P, 
         tl.cex = 0.5, insig = "blank", number.cex = 0.3, cl.cex = 0.75)
dev.off()

# Here, we chose to remove nesting substrate: ground-open, nesting behavior:
# solitary, breeding habitats urban suburban, shrubland, and woodland, and daily
# activity period: crepuscular (a total of six traits).
binaries %>%
  dplyr::select(-nesting_substrate_ground_open_va, 
                -nesting_behavior_solitary_st, 
                -daily_activity_period_crepuscular_pe,
         -breeding_habitat_urban_suburban_pe, -breeding_habitat_shrubland_pe,
         -breeding_habitat_woodland_pe) -> binaries_uncorr
# Unite the binaries and non-binaries data frames.
bind_cols(non_binaries_uncorr, trait_matrix$early_developmental_mode_st, 
          binaries_uncorr) -> trait_matrix_uncorr
"early_developmental_mode_st" -> colnames(trait_matrix_uncorr)[23]

# Lastly, check if any of the binaries are correlated with continuous ones. Use 
# point biserial correlation coefficient for the purpose.
as.data.frame(trait_matrix_uncorr) -> trait_matrix_uncorr_df
for(a in 1:22){
  for(b in 24:75){
    if(abs(biserial.cor(x = trait_matrix_uncorr_df[, a], 
                        y = trait_matrix_uncorr_df[, b], 
                        use = "complete.obs")) > 0.7){
      print(c(colnames(trait_matrix_uncorr_df)[a], 
              colnames(trait_matrix_uncorr_df)[b]))
    }
  }
}
# It looks like the problematic trait is the daily activity plasticity, so let's
# remove it, too.
trait_matrix_uncorr_df[, -9] -> trait_matrix_uncorr_df
# At last, we are done with the correlated trait removal, and, in the end, we 
# are left with 74 traits. Rename the resulting data frames.
trait_matrix_uncorr_df -> trait_matrix_final; remove(trait_matrix_uncorr_df)
# Our trait matrix is free of correlated traits now, and it is ready to be used
# in functional diversity calculations. 


# 3: FUNCTIONAL DIVERSITY ####


# Ready the community matrices so that both data frames have the same row names, 
# that is, the species names.
Negate("%in%") -> "%ni%" # Opposite of %in%.
rownames(study_comm_1998) %ni% rownames(study_comm_2018) -> row_indexer_1998
rownames(study_comm_1998)[row_indexer_1998] -> 
  rownames(study_comm_2018[98:124, ])
# Populate the newly added rows with zeroes.
0 -> study_comm_2018[98:124, ]
# Repeat the same steps for the 1998's data.
rownames(study_comm_2018) %ni% rownames(study_comm_1998) -> row_indexer_2018
rownames(study_comm_2018)[row_indexer_2018] -> 
  rownames(study_comm_1998[121:124, ])
0 -> study_comm_1998[121:124, ]
study_comm_1998[order(rownames(study_comm_1998)), ] -> study_comm_1998
study_comm_2018[order(rownames(study_comm_2018)), ] -> study_comm_2018
# One last check.
rownames(study_comm_1998) == rownames(study_comm_2018)
colnames(study_comm_1998) == colnames(study_comm_2018)

# Read in the weights file. See the supplementary material for information on
# how the weights were calculated.
read_xlsx(path = "functional diversity/data/trait_weights.xlsx") -> weights
# Check the order of the traits.
trait_matrix_final[, order(colnames(trait_matrix_final))] -> trait_matrix_final
weights[order(weights$TRAIT), ] -> weights
weights$TRAIT == colnames(trait_matrix_final)


## 3.1: ALPHA & GAMMA FUNCTIONAL DIVERSITY ####


# Firstly, find out which species are present in 1998 and 2018 in the whole 
# basin.
data.frame(X1998 = rowSums(study_comm_1998), 
           X2018 = rowSums(study_comm_2018)) -> basin_incidences
# Binarize this data frame.
func_binarize <- function(df){
  if(class(df) != "data.frame"){
    stop("df must be a data frame")
  }
  for(i in 1:nrow(df))
    for(j in 1:ncol(df)){
      if(df[i, j] > 0){
        1 -> df[i, j]
      }
    }
  return(df)
}
func_binarize(basin_incidences) -> basin_incidences
# Transpose for use in dbFD function.
as.data.frame(t(basin_incidences)) -> basin_incidences

# Matching of species names in the community matrix and traits matrix can be a
# little problematic in dbFD. Therefore, we will assign them simple names.
colnames(basin_incidences) == trait_matrix$species
trait_matrix$species -> rownames(trait_matrix_final)
rownames(trait_matrix_final) == colnames(basin_incidences) # No problem here.

# We can now proceed with the functional diversity calculations at gamma scale.
dbFD(x = trait_matrix_final, a = basin_incidences, w = weights$WEIGHT, 
     corr = "cailliez", calc.CWM = T, calc.FDiv = F, 
     CWM.type = "all" , print.pco = T, m = 12) -> FD_incidences
FD_incidences$FRic
100 - (3.606616e-10/1.054840e-09)*100 # The decline in gamma-scale Functional 
# Richness (FRic) is 65.8%, which is 3.4 times larger than the decline in 
# taxonomic richness (19.2%).
# For this calculation, we could use a maximum of 12 PCoA axes, which explained
# 49.3% of the total variation in the data.

# Extract the community weighted means (CWMs) for each trait and export it as a
# .xlsx file for further inspections. Also, export the PCoA axes.
FD_incidences$CWM -> CWMs
write.xlsx(x = CWMs, file = "functional diversity/output/CWMs.xlsx")
FD_incidences$x.axes -> pcoa_axes

# Also, calculate the AUC for this dimensionality reduction.
gowdis(x = trait_matrix_final, w = weights$WEIGHT) -> gower_dists
coranking(Xi = gower_dists, X = pcoa_axes[, 1:12], input_Xi = "dist", 
          input_X = "data", use = "C") -> x1
R_NX(x1) -> x2
AUC_ln_K(x2) # AUC: 0.75 for the dimensionality reduction. A critical approach 
# to the dimensionality reduction, as well as the application of the methods 
# used in Mouillot et al. (2021), is given in section 3.3.

# Next, calculate alpha-scale FD metrics for the study years and squares. To do
# that, we will have to merge the two community matrices into one.
study_comm_1998 -> study_comm_1998_merge
paste(colnames(study_comm_1998), "1998", 
      sep = "_") -> colnames(study_comm_1998_merge)
study_comm_2018 -> study_comm_2018_merge
paste(colnames(study_comm_2018), "2018", 
      sep = "_") -> colnames(study_comm_2018_merge)
study_comm_1998_merge[, order(colnames(study_comm_1998_merge))] -> 
  study_comm_1998_merge
study_comm_2018_merge[, order(colnames(study_comm_2018_merge))] -> 
  study_comm_2018_merge
cbind.data.frame(study_comm_1998_merge, study_comm_2018_merge) -> study_comms

# Check the species richness in the squares because we need to satisfy s>t,
# which will define the maximum number of PCoA axes that can be used in 
# alpha-scale FD calculations.
hist(colSums(study_comms), breaks = seq(from = 0, to = 120, by = 5), 
     xlim = c(0, 85), ylim = c(0, 25), labels = T)
abline(v = 5, col = "red", lty = 2)
hist(colSums(study_comms), breaks = seq(from = 0, to = 120, by = 1), 
     xlim = c(0, 85), ylim = c(0, 25), labels = T)
abline(v = 5, col = "red", lty = 2)

# Ready the data frame and do the final checks.
as.data.frame(t(study_comms)) -> study_comms
colnames(study_comms) == rownames(trait_matrix_final)

# It looks like setting m = 4 seems to be a good compromise between the number 
# of grid squares that can be included in the study and the quality of the 
# dimensionality reduction. The same number of PCoA axes is the maximum number 
# of dimensions that can be included in the beta diversity calculations (see the
# section 3.2).
which(rowSums(study_comms) <= 4)
length(which(rowSums(study_comms) <= 4))
# 12 such communities in 2018, and 1 in 1998. Remove them and save as a new data
# frame.
which(rowSums(study_comms) <= 4) -> indexer_S4
study_comms[-indexer_S4, ] -> study_comms_S4

# We should compare the same exact squares, because removal of some of the 
# squares result in an uneven case for the study years. Therefore, use GLMMs to
# make the comparison.
which(rowSums(study_comms) <= 4) - 35 -> indexer_shared
indexer_shared[- c(1, 2)] -> indexer_shared
study_comms_S4[-indexer_shared, ] -> study_comms_shared

# Check if any of the species does not occur in the any community.
which(colSums(study_comms_shared) == 0) -> zero_species
study_comms_shared[, -zero_species] -> study_comms_shared
# Remove the related columns from the trait matrix as well.
trait_matrix_final[-zero_species, ] -> trait_matrix_shared

# Finally, calculate FD at the grid square scale.
dbFD(x = trait_matrix_shared, a = study_comms_shared, w = weights$WEIGHT, 
     corr = "cailliez", calc.CWM = T, CWM.type = "all", print.pco = T,
     m = 4, stand.FRic = T) -> FD_alpha_shared
# The representation of the trait variation is 31.3%. 
# Calculate AUC.
gowdis(x = trait_matrix_shared, 
       w = weights$WEIGHT) -> gower_dists_shared
coranking(Xi = gower_dists_shared, X = pcoa_axes[-zero_species, 1:4], 
          input_Xi = "dist", input_X = "data", use = "C") -> x1
R_NX(x1) -> x2
AUC_ln_K(x2) # AUC: 0.54 for the shared dataset.

# Extract FRic and FDis for the years.
FD_alpha_shared$FRic[1:22] -> FRic_1998
FD_alpha_shared$FRic[23:length(FD_alpha_shared$FRic)] -> FRic_2018	
data.frame(year = as.factor(c(rep("1998", 22), rep("2018", 22))),
           FRic = as.vector(c(FRic_1998	, FRic_2018))) -> FRic_df		
FD_alpha_shared$FDis[1:22] -> FDis_1998			
FD_alpha_shared$FDis[23:length(FD_alpha_shared$FDis)] -> FDis_2018			
data.frame(year = as.factor(c(rep("1998", 22), rep("2018", 22))),
           FDis = as.vector(c(FDis_1998, FDis_2018))) -> FDis_df		
# Plot the data.
ggplot() +
  geom_boxplot(data = FRic_df, aes(x = year, y = FRic, fill = year),
               alpha = 0.8, notch = T) +
  scale_fill_manual(values = c("deepskyblue", "deepskyblue4")) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  ylab(label = "Functional Richness (FRic)") + xlab(label = NULL) +
  theme(axis.title.x = element_text(size = 12)) +
  guides(fill = guide_legend("Year")) -> p_fric_alpha
p_fric_alpha
ggplot() +
  geom_boxplot(data = FDis_df, aes(x = year, y = FDis, fill = year), 
               alpha = 0.8, notch = T) +
  scale_fill_manual(values = c("deepskyblue", "deepskyblue4")) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  ylab(label = "Functional Dispersion (FDis)") + xlab(label = NULL) +
  theme(axis.title.x = element_text(size = 12)) +
  guides(fill = guide_legend("Year")) -> p_fdis_alpha
p_fdis_alpha
# The differences look more striking now. Make the formal comparison by using
# GLMMs.
as.factor(rep((paste("square", 1:22, sep = "_")), 2)) -> FRic_df$ID
m_fric_alpha <- glmmTMB(data = FRic_df, formula = FRic ~ year, 
                        family = gaussian(link = "identity"))
summary(m_fric_alpha) # The difference is significant with 1998's communities
# having higher FRic (p = 0.0167*). Validate the model.

# Model validation #
# 1) Check the distribution of the residuals. 
boxplot(residuals(m_fric_alpha)); abline(h = 0, col = "red", lty = 2)
# Nearly perfect.
# 2) Plot the residuals against the only predictor.
boxplot(residuals(m_fric_alpha) ~ FRic_df$year)
abline(h = 0, col = "red", lty = 2)
# The model performs only marginally worse at 2018 level. No problem here.
# 3) Check the model for autocorrelation.
par(mfrow = c(1, 2))
acf(residuals(m_fric_alpha)); pacf(residuals(m_fric_alpha))
par(mfrow = c(1, 1))
# No problem here except for a random spike in PACF plot, which can be ignored.
# Also check for spatial autocorrelation.
read.csv2("functional diversity/data/relative_coords.csv", 
          header = T) -> rel_coords
rep(rel_coords$x, 2) -> FRic_df$x
rep(rel_coords$y, 2) -> FRic_df$y
set.seed(13)
simulateResiduals(m_fric_alpha, n = 1000) -> sim_fric_alpha_2
set.seed(13); recalculateResiduals(sim_fric_alpha_2, 
                     group = FRic_df$ID) -> sim_fric_alpha_2
set.seed(13); testSpatialAutocorrelation(sim_fric_alpha_2, 
                           x = FRic_df[1:22,]$x,
                           y = FRic_df[1:22,]$y)
# The residuals have no spatial autocorrelation.
# 4) Use the DHARMa's utilities to check the model. 
set.seed(13); simulateResiduals(m_fric_alpha, n = 1000) -> sim_fric_alpha
set.seed(13); testResiduals(sim_fric_alpha)
# The model is valid and the FRic of communities of 1998 are significantly 
# higher (p = 0.0167*).

# Now, model the FDis data. 
as.factor(rep((paste("square", 1:22, sep = "_")), 2)) -> FDis_df$ID
m_fdis_alpha <- glmmTMB(data = FDis_df, formula = FDis ~ year, 
                        family = gaussian(link = "identity"))
summary(m_fdis_alpha) # The difference is not significant.


## 3.2: FUNCTIONAL DIVERSITY & COMPOSITIONAL STABILITY ####


# Calculate various FD metrics for the 1998's communities to check whether FD
# provides any compositional stability.
study_comm_1998[order(colnames(study_comm_1998))] -> study_comm_1998_ordered
which(colSums(study_comm_1998_ordered) < 3) # Because FD metrics we will use 
# require a minimum of three species.
study_comm_1998_ordered[, -6] -> study_comm_1998_for_FD
min(colSums(study_comm_1998_for_FD))
hist(colSums(study_comm_1998_for_FD), breaks = seq(from = 0, to = 85, by = 1),
     labels = T, ylim = c(0, 5))
abline(v = 12, col = "red", lty = 2); abline(v = 5, col = "blue", lty = 2)
# There are six squares (out of 33) that have fewer than 12 (the number of PCoA
# axes used in gamma FD) species and three squares that have fewer than four 
# (the number of PCoA axes used in alpha and beta FD) species.
# Try the first 12 PCoA axes first.
which(colSums(study_comm_1998_for_FD) <= 12) -> indexer_S12
study_comm_1998_for_FD[, -indexer_S12] -> study_comm_1998_for_FD
which(rowSums(study_comm_1998_for_FD) == 0) -> zero_species_1998_FD
study_comm_1998_for_FD[-zero_species_1998_FD, ] -> study_comm_1998_for_FD
trait_matrix_final[-zero_species_1998_FD, ] -> trait_matrix_1998_FD
as.data.frame(t(study_comm_1998_for_FD)) -> study_comm_1998_for_FD
study_comm_1998_for_FD[order(rownames(study_comm_1998_for_FD)), ] ->
  study_comm_1998_for_FD

# Calculate the FD metrics.
dbFD(x = trait_matrix_1998_FD, a = study_comm_1998_for_FD, w = weights$WEIGHT,
     corr = "cailliez", m = 12, print.pco = T, 
     CWM.type = "all", calc.CWM = T) -> FD_1998_for_stability
# The percent representation of the variation after the reduction is 49.8%. 
# Calculate AUC for this dimensionality reduction.
FD_1998_for_stability$x.axes -> pcoa_axes_1998_FD
gowdis(x = trait_matrix_1998_FD, 
       w = weights$WEIGHT) -> gower_dists_1998_FD
coranking(Xi = gower_dists_1998_FD, X = pcoa_axes_1998_FD[, 1:12], 
          input_Xi = "dist", input_X = "data", use = "C") -> x1
R_NX(x1) -> x2
AUC_ln_K(x2) # AUC: 0.75 for this dimensionality reduction.

# Read in the pairwise taxonomic beta diversity that was calculated in the other
# script.
read.csv2("taxonomic diversity/output/self_beta.csv", ) -> self_tax_beta
self_tax_beta[, -1] -> self_tax_beta
NA -> self_tax_beta[self_tax_beta == "NaN"]
self_tax_beta[-c(6, indexer_S12), ] -> self_tax_beta_1998_FD
data.frame(Bsim = as.numeric(self_tax_beta_1998_FD$beta.sim),
           Bsne = as.numeric(self_tax_beta_1998_FD$beta.sne),
           Bsor = as.numeric(self_tax_beta_1998_FD$beta.sor),
           FRic = as.numeric(FD_1998_for_stability$FRic),
           FDis = as.numeric(FD_1998_for_stability$FDis),
           FDiv = as.numeric(FD_1998_for_stability$FDiv),
           FEve = as.numeric(FD_1998_for_stability$FEve)) -> stability_df
rcorr(as.matrix(stability_df)) -> stability_corr_matrix
corrplot(corr = stability_corr_matrix$r, type = "upper", diag = F, outline = T, 
         tl.col = "black", method = "number", p.mat = stability_corr_matrix$P, 
         tl.cex = 0.75, insig = "blank", number.cex = 1, cl.cex = 1)
# FEve (rho = -0.46) and FDis (rho = -0.42) seem to provide a slight immunity to
# compositional change and provide some compositional stability.

# Check the relationship by using a GLM to make a formal conclusion about the 
# relationship.
stability_df$Bsor - 0.00001 -> stability_df$Bsor
m_fdis_stability <- glmmTMB(data = stability_df, formula = Bsor ~ FDis, 
                            family = beta_family(link = "logit"))
summary(m_fdis_stability) # The relationship is significant (p = 0.0355*). Let's
# check the model.

# Model validation #
# 1) Check the distribution of the residuals.
boxplot(residuals(m_fdis_stability)); abline(h = 0, col = "red", lty = 2)
# The distribution is not perfect, but it may be acceptable considering the low
# sample size.
# 2) Plot the residuals against the predictor.
plot(stability_df$FDis, residuals(m_fdis_stability))
abline(h = 0, col = "red", lty = 2)
# No problem here.
# 3) Check for residual autocorrelation.
par(mfrow = c(1, 2))
acf(residuals(m_fdis_stability)); pacf(residuals(m_fdis_stability))
par(mfrow = c(1, 1))
# No dependence issues.
# 4) Use the DHARMa utilities.
set.seed(13); simulateResiduals(m_fdis_stability, 
                                n = 1000) -> sim_fdis_stability
set.seed(13); testResiduals(sim_fdis_stability)
# We conclude that the model is valid with FDis providing significant 
# compositional stability to the communities (p = 0.00569**). 

# Build another model for FRic.
m_fric_stability <- glmmTMB(data = stability_df, formula = Bsor ~ scale(FRic), 
                            family = beta_family(link = "logit"))
summary(m_fric_stability)
# The effect is not significant.


## 3.3: BETA FUNCTIONAL DIVERSITY ####


# We will now calculate pairwise functional beta diversity for the study years.
# Start with calculating pairwise functional beta diversity.


### 3.3.1: PAIRWISE BETA FUNCTIONAL DIVERSITY ####


# We will use the study_comms_shared data frame for the calculations, but we 
# need to split the data frame first.
study_comms_shared[1:21 ,] -> comms_1998_shared
study_comms_shared[22:42 ,] -> comms_2018_shared
FD_alpha_shared$x.axes -> pcoa_axes_shared
functional.betapart.core.pairwise(x = comms_1998_shared, 
                                  traits = pcoa_axes_shared[1:4], 
                                  return.details = F, parallel = T, 
                                  progress = T) -> FD_beta_1998
functional.betapart.core.pairwise(x = comms_2018_shared, 
                                  traits = pcoa_axes_shared[1:4],
                                  return.details = F, parallel = T, 
                                  progress = T) -> FD_beta_2018
# We can now use the functional.beta.pair() function to calculate the turnover, 
# nestedness components of functional beta diversity, as well as the total beta
# functional diversity.
functional.beta.pair(FD_beta_1998, index.family = "sorensen") -> FD_beta_1998
functional.beta.pair(FD_beta_2018, index.family = "sorensen") -> FD_beta_2018

# Plot the components of functional beta diversity. To do that, ready the data
# in "tidy" format.
data.frame(year = as.factor(c(rep("1998", 630), rep("2018", 630))),
           value = c(FD_beta_1998$funct.beta.sim, FD_beta_1998$funct.beta.sne,
                     FD_beta_1998$funct.beta.sor, FD_beta_2018$funct.beta.sim,
                     FD_beta_2018$funct.beta.sne, FD_beta_2018$funct.beta.sor),
           index = as.factor(c(rep("FBsim", 210), rep("FBsne", 210), 
                               rep("FBsor", 210), rep("FBsim", 210),
                               rep("FBsne", 210), 
                               rep("FBsor", 210)))) -> beta_data
# Generate the plot.
 ggplot() +
  geom_boxplot(data = beta_data, aes(x = index, y = value, fill = year), 
               alpha = 0.8, notch = T) +
  scale_fill_manual(values = c("deepskyblue", "deepskyblue4")) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  ylab(label = "Index value") + xlab(label = "Index") +
  theme(axis.title.x = element_text(size = 12)) +
  scale_y_continuous(position = "right") +
  guides(fill = guide_legend("Year")) -> p_functional_beta
p_functional_beta
# Save the plot (optional).
 ggsave(p_functional_beta, filename ="functional diversity/figs/func_beta.tiff",
       units = "mm", dpi = 400, width = 120, height = 75)


# Compare the components of beta functional diversity by using Wilcoxon signed
# rank test because of the reasons listed in section 3.2 of the taxonomic
# diversity script.
wilcox.test(x = as.vector(FD_beta_1998$funct.beta.sim),  
            y = as.vector(FD_beta_2018$funct.beta.sim),
            alternative = "two.sided")
# 2018's communities have significantly higher FBsim (p = 8.242e-11).
wilcox.test(x = as.vector(FD_beta_1998$funct.beta.sne),  
            y = as.vector(FD_beta_2018$funct.beta.sne),
            alternative = "two.sided")
# The difference is not significant (p = 0.2093).
wilcox.test(x = as.vector(FD_beta_1998$funct.beta.sor),  
            y = as.vector(FD_beta_2018$funct.beta.sor),
            alternative = "two.sided")
# 2018's communities have significantly higher FBsor (p = 1.665e-14).


### 3.3.2: DISTANCE DECAY MODELS ####


# Get the distances first.
dist(x = FRic_df[1:21, 4:5], method = "euclidean") -> utm_distances
as.matrix(utm_distances) -> utm_distances
utm_distances # Just to check.
as.dist(m = utm_distances, diag = F, upper = T) -> utm_distances

# Check the beta functional diversity data.
range(FD_beta_1998$funct.beta.sim) # There are zeroes and the data is 
# continuous. However, we will need to do something about the zeroes because we
# will use log link functions in the distance decay models. 
func_log_trans <- function(x){
  (x * 32 + 0.5)/ 33
}
lapply(X = FD_beta_1998, FUN = func_log_trans) -> FD_beta_1998
lapply(X = FD_beta_2018, FUN = func_log_trans) -> FD_beta_2018

# Build the distance decay models now. For that purpose, define a function that
# This time, exponential models are prioritized in cases where delta AIC < 2.
func_decay_select <- function(input){
  if(class(input) != "list"){
    stop("Error: input must be a output from the beta.pair function")
    # An extra step of not-so-useful precaution.
  }
    power_models <- list(); exp_models <- list(); selected_models <- list()
    # Define some empty lists that will be populated later.
    for(i in 1:length(input)){
      set.seed(169)
      if(i == 1){
        index_type <- "sim"
      } else if(i == 2){
        index_type <- "sne"
      } else {
        index_type <- "sor" # The indices are saved in this order after the
      }                     # beta.pair() function is run.
              decay.model(input[[i]], utm_distances, model.type = "power", 
        y.type = "dissimilarities", perm = 1000) -> power_models[[i]]
              decay.model(input[[i]], utm_distances, model.type = "exponential", 
        y.type = "dissimilarities", perm = 1000) -> exp_models[[i]]
              if(power_models[[i]]$model$aic - exp_models[[i]]$model$aic > 2){
                mod_type <- "exponential"
                print(paste("p =",exp_models[[i]]$p.value, "for",index_type,
                            "model, type =",mod_type))
                selected_models[[i]] <- exp_models[[i]]
              } else {
                mod_type <- "power"
                print(paste("p =",power_models[[i]]$p.value, "for",index_type,
                            "model, type =",mod_type))
                selected_models[[i]] <- power_models[[i]]
                # This chunk builds and selects between the power and 
                # exponential models and stores them in selected_models object.
                # The selection is AIC-based. 
              }
    }
    return(selected_models)
}

# Run this function on 1998's and 2018's pairwise beta diversity matrices.
set.seed(13); func_decay_select(FD_beta_1998) -> ddms_1998
set.seed(13); func_decay_select(FD_beta_2018) -> ddms_2018
# None of the distance decay models (ddms) except for FBsim (slope = 0.028) and
# FBsor (slope = 0.025) for 2018 has significant slopes. Let's check if the 
# model with a significant slope is a valid model.

# Model validation #
# 1) Check the distribution of the residuals around zero.
boxplot(ddms_1998[[3]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# The distribution is not ideal, but it may be acceptable.
# 2) Plot the residuals against the only predictor in the model: distances to 
# see if there are any patterns in the plot. 
plot(log(utm_distances), 
      ddms_1998[[3]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# There seems to be no significant patterns here. We will assume that the model
# is valid.

# Produce a graph to be used in the manuscript (optional).
 tiff(filename = "functional diversity/figs/functional_ddms.tiff", 
      width = 230, height = 80, units = "mm", res = 600)
 par(mfrow = c(1, 3))
 for(i in 1:3){
   plot.decay(ddms_1998[[i]], remove.dots = T, 
             col = "deepskyblue", ylim = c(0,1), lwd = 3.5)
    plot.decay(ddms_2018[[i]], remove.dots = T, col = "deepskyblue4", 
               ylim = c(0,1), lwd = 3.5, add = T)
 }
 dev.off()

# Use the bootstrapped model parameters to compare the slopes of the FBsor 
# ddms.
set.seed(13); boot.coefs.decay(ddms_1998[[1]], R = 10000) -> boot_ddm_1998_sim
set.seed(13); boot.coefs.decay(ddms_2018[[1]], R = 10000) -> boot_ddm_2018_sim
set.seed(13); boot.coefs.decay(ddms_1998[[3]], R = 10000) -> boot_ddm_1998_sor
set.seed(13); boot.coefs.decay(ddms_2018[[3]], R = 10000) -> boot_ddm_2018_sor

# Plot the model parameters for FBsim and FBsor side-by-side.
adjustcolor(col = "deepskyblue", alpha.f = 0.9) -> blue_tr
adjustcolor(col = "deepskyblue4", alpha.f = 0.9) -> dark_blue_tr
range(c(boot_ddm_1998_sim$boot.coefs[,2], 
        boot_ddm_2018_sim$boot.coefs[,2])) -> current_range
boxplot(boot_ddm_2018_sim$boot.coefs[,2], col = dark_blue_tr,
        ylim = current_range,
        main = "FBsim Slope Bootstrapped - Blue:1998, Dark blue:2018")
boxplot(boot_ddm_1998_sim$boot.coefs[,2], col = blue_tr,
        ylim = current_range, add = T)
range(c(boot_ddm_1998_sor$boot.coefs[,2], 
        boot_ddm_2018_sor$boot.coefs[,2])) -> current_range
boxplot(boot_ddm_2018_sor$boot.coefs[,2], col = dark_blue_tr,
        ylim = current_range,
        main = "FBsor Slope Bootstrapped - Blue:1998, Dark blue:2018")
boxplot(boot_ddm_1998_sor$boot.coefs[,2], col = blue_tr,
        ylim = current_range, add = T)

# Plot a proper box plot now to be included in the paper.
data.frame(slope = c(boot_ddm_1998_sim$boot.coefs[, 2], 
                     boot_ddm_2018_sim$boot.coefs[, 2]), 
           year = as.factor(c(rep("1998", 10000), 
                              rep("2018", 10000)))) -> FBsim_slopes
ggplot() +
  geom_boxplot(data = FBsim_slopes, 
               aes(x = year, y = slope, fill = year), alpha = 0.8) +
  scale_fill_manual(values = c("deepskyblue", "deepskyblue4")) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  theme(axis.title.x = element_text(size = 12)) + 
  ylab("Bootstrapped DDM slope") + xlab("Year") +
  guides(fill = guide_legend("Year")) -> p_ddm_slope_FBsim
p_ddm_slope_FBsim
# Save the plot (optional). 
# ggsave(filename = "functional diversity/figs/p_ddm_slope_FBsim.tiff", 
#       device = "tiff", dpi = 600, height = 4, width = 3.4)
data.frame(slope = c(boot_ddm_1998_sor$boot.coefs[, 2], 
                     boot_ddm_2018_sor$boot.coefs[, 2]), 
           year = as.factor(c(rep("1998", 10000), 
                              rep("2018", 10000)))) -> FBsor_slopes
ggplot() +
  geom_boxplot(data = FBsor_slopes, 
               aes(x = year, y = slope, fill = year), alpha = 0.8) +
  scale_fill_manual(values = c("deepskyblue", "deepskyblue4")) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  theme(axis.title.x = element_text(size = 12)) + 
  ylab("Bootstrapped DDM slope") + xlab("Year") +
  guides(fill = guide_legend("Year")) -> p_ddm_slope_FBsor
p_ddm_slope_FBsor
# Save the plot (optional). 
# ggsave(filename = "functional diversity/figs/p_ddm_slope_FBsor.tiff", 
#       device = "tiff", dpi = 600, height = 4, width = 3.4)

# Make the formal comparison of slopes by calculating p-value as the proportion 
# of times one of the slopes is larger than the other.
FBsim_slopes$slope[1:10000] - FBsim_slopes$slope[10001:20000] -> slope_dif_FBsim
length(which(slope_dif_FBsim < 0)) 
length(which(slope_dif_FBsim > 0)) # Smaller one.
(length(which(slope_dif_FBsim > 0)) / 10000)*2
# p-value for this comparison is 0.11 and the difference is not significant.
FBsor_slopes$slope[1:10000] - FBsor_slopes$slope[10001:20000] -> slope_dif_FBsor
length(which(slope_dif_FBsor < 0)) 
length(which(slope_dif_FBsor > 0)) # Smaller one.
(length(which(slope_dif_FBsor > 0)) / 10000)*2
# p-value for this comparison is 0.25 and the difference is not significant.

# Do the same for the intercepts.
data.frame(intercept = c(boot_ddm_1998_sim$boot.coefs[, 1], 
                     boot_ddm_2018_sim$boot.coefs[, 1]), 
           year = as.factor(c(rep("1998", 10000), 
                              rep("2018", 10000)))) -> FBsim_intercept
FBsor_intercept$intercept[1:10000] - FBsim_intercept$intercept[10001:20000] -> 
  intercept_dif_FBsim
length(which(intercept_dif_FBsim < 0)) # Smaller one.
length(which(intercept_dif_FBsim > 0)) 
(length(which(intercept_dif_FBsim < 0)) / 10000)*2
# The difference is not significant (p = 0.64).
data.frame(intercept = c(boot_ddm_1998_sim$boot.coefs[, 1], 
                     boot_ddm_2018_sim$boot.coefs[, 1]), 
           year = as.factor(c(rep("1998", 10000), 
                              rep("2018", 10000)))) -> FBsor_intercept
FBsor_intercept$intercept[1:10000] - FBsor_intercept$intercept[10001:20000] -> 
  intercept_dif_FBsor
length(which(intercept_dif_FBsor < 0)) # Smaller one.
length(which(intercept_dif_FBsor > 0)) 
(length(which(intercept_dif_FBsor < 0)) / 10000)*2
# The difference is not significant (p = 0.44).


## 3.3: QUALITY OF DIMENSIONALITY REDUCTION ####


# Write a function that will calculate the AUC value for each PCoA axes included
# along with a net benefit (as described in Mouillot et al. (2021)) to calculate
# the elbow inflexion point for the SHARED dataset.
seq(from = 0, to = 1, by = 1/20) -> diag_line_points
func_AUC <- function(max_PCoA){
  if(class(max_PCoA) != "numeric"){
    stop("max_PCoa must be a number")
  }
  vector() -> AUC_values_shared; vector() -> elbow_AUC_shared
  data.frame(AUC = rep(NA, 20), elbow = rep(NA, 20)) -> AUC_df
  for(i in 1:max_PCoA){
  coranking(Xi = gower_dists_shared, X = pcoa_axes_shared[, 1:i], 
            input_Xi = "dist", input_X = "data", use = "C") -> x1
  R_NX(x1) -> x2
  AUC_ln_K(x2) -> AUC_values_shared[i]
  AUC_values_shared[i] - diag_line_points[i] -> elbow_AUC_shared[i]
  AUC_values_shared[i] -> AUC_df[i, 1]
  elbow_AUC_shared[i] -> AUC_df[i, 2]
  }
  return(AUC_df) 
}
# Try it on the first 20 PCoA axes.
func_AUC(20) -> AUC_first_20

# Add the number of PCoA axes included to the data frame.
1:20-> AUC_first_20$PCoA_no

# Now produce a beautiful plot to be included in the supplementary material
# (optional).
ggplot() +
  geom_line(data = AUC_first_20, aes(x = PCoA_no, y = AUC)) +
  geom_point(data = AUC_first_20, aes(x = PCoA_no, y = AUC), size = 2) +
  geom_line(data = AUC_first_20, aes(x = PCoA_no, y = elbow), col = "gray45") +
  geom_point(data = AUC_first_20, aes(x = PCoA_no, y = elbow), col = "gray45",
             size = 2) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  geom_point(data = data.frame(x = 5, y = AUC_first_20$AUC[5]), aes(x, y), 
             pch = 1, size = 6, stroke = 0.75) +
  geom_point(data = data.frame(x = 8, y = AUC_first_20$AUC[8]), aes(x, y), 
             pch = 1, size = 6, stroke = 0.75) +
  geom_point(data = data.frame(x = 4, y = AUC_first_20$AUC[4]), aes(x, y), 
             pch = 1, size = 6, colour = "firebrick1", stroke = 0.75) +
  geom_segment(aes(x = 5, y = 0, xend = 5, yend = AUC_first_20$AUC[5]), 
               lty = 2, col = "gray25") +
  geom_segment(aes(x = 8, y = 0, xend = 8, yend = AUC_first_20$AUC[8]), 
               lty = 2, col = "gray25") +
  geom_segment(aes(x = 4, y = 0, xend = 4, yend = AUC_first_20$AUC[4]), 
               lty = 2, col = "firebrick1") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 20), 
                     breaks = seq(from = 0, to = 20, by = 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  annotate(geom = "text", x = 3.5, y = 0.60, label = "AUC = 0.54", 
           col = "firebrick2") +
  annotate(geom = "text", x = 4.8, y = 0.66, label = "AUC = 0.61") +
  annotate(geom = "text", x = 7.9, y = 0.76, label = "AUC = 0.71") +
  xlab("Number of PCoA Axes") -> p_dim_red_shared
p_dim_red_shared

# Save the plot (optional).
 ggsave("functional diversity/figs/dim_reduction_shared.tiff", device = "tiff", 
       dpi = 400, width = 8, height = 4, units = "in")

# Repeat all of the above for the complete dataset, where we also used the first
# 12 PCoA axes.
func_AUC <- function(max_PCoA){
  if(class(max_PCoA) != "numeric"){
    stop("max_PCoa must be a number")
  }
  vector() -> AUC_values; vector() -> elbow_AUC
  data.frame(AUC = rep(NA, 20), elbow = rep(NA, 20)) -> AUC_df
  for(i in 1:max_PCoA){
  coranking(Xi = gower_dists, X = pcoa_axes[, 1:i], 
            input_Xi = "dist", input_X = "data", use = "C") -> x1
  R_NX(x1) -> x2
  AUC_ln_K(x2) -> AUC_values[i]
  AUC_values[i] - diag_line_points[i] -> elbow_AUC[i]
  AUC_values[i] -> AUC_df[i, 1]
  elbow_AUC[i] -> AUC_df[i, 2]
  }
  return(AUC_df)
}
# Try it on the first 20 PCoA axes.
func_AUC(20) -> AUC_first_20

# Add the number of PCoA axes included to the data frame.
1:20-> AUC_first_20$PCoA_no

ggplot() +
  geom_line(data = AUC_first_20, aes(x = PCoA_no, y = AUC)) +
  geom_point(data = AUC_first_20, aes(x = PCoA_no, y = AUC), size = 2) +
  geom_line(data = AUC_first_20, aes(x = PCoA_no, y = elbow), col = "gray45") +
  geom_point(data = AUC_first_20, aes(x = PCoA_no, y = elbow), col = "gray45",
             size = 2) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  geom_point(data = data.frame(x = 5, y = AUC_first_20$AUC[5]), aes(x, y), 
             pch = 1, size = 6, stroke = 0.75) +
  geom_point(data = data.frame(x = 8, y = AUC_first_20$AUC[8]), aes(x, y), 
             pch = 1, size = 6, stroke = 0.75) +
  geom_point(data = data.frame(x = 12, y = AUC_first_20$AUC[12]), aes(x, y), 
             pch = 1, size = 6, colour = "green4", stroke = 0.75) +
  geom_segment(aes(x = 5, y = 0, xend = 5, yend = AUC_first_20$AUC[5]), 
               lty = 2, col = "gray25") +
  geom_segment(aes(x = 8, y = 0, xend = 8, yend = AUC_first_20$AUC[8]), 
               lty = 2, col = "gray25") +
  geom_segment(aes(x = 12, y = 0, xend = 12, yend = AUC_first_20$AUC[12]), 
               lty = 2, col = "green4") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 20), 
                     breaks = seq(from = 0, to = 20, by = 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  annotate(geom = "text", x = 12, y = 0.8, label = "AUC = 0.75", 
           col = "green4") +
  annotate(geom = "text", x = 4.8, y = 0.67, label = "AUC = 0.61") +
  annotate(geom = "text", x = 7.9, y = 0.76, label = "AUC = 0.71") +
  xlab("Number of PCoA Axes") -> p_dim_red
p_dim_red

# Save the plot (optional).
# ggsave("functional diversity/figs/dim_reduction.tiff", device = "tiff", 
#       dpi = 400, width = 8, height = 4, units = "in")
 
# Also, calculate the correlation between the distances calculated by using the
# complete dataset and 12 and 4 PCoA axes.
gowdis(x = trait_matrix_shared, w = weights$WEIGHT) -> gower_dists_full
gowdis(x = pcoa_axes_shared[, 1:12]) -> gower_dists_12 # First 12 PCoA axes
gowdis(x = pcoa_axes_shared[, 1:4]) -> gower_dists_4   # First 4 PCoA axes

# Now, calculate the correlation between the matrices.
mantel(gower_dists_full, gower_dists_12, method = "spearman", 
       permutations = 10000)
# The distance matrices are strongly and positively correlated (rho = 0.83, 
# p< 0.001***)
mantel(gower_dists_full, gower_dists_4, method = "spearman", 
       permutations = 10000)
# The distance matrices are strongly and positively correlated (rho = 0.85, 
# p< 0.001***)
# According to the Mantel tests' results, both of the matrices are strongly and
# significantly correlated with the original distance matrix, which was
# calculated by using all of the traits.

## 3.4: FUNCTIONAL ORIGINALITY ####


# Remove the communities that have zero species in 2018, because we aim to 
# compare Functional Originality (FOri) of the extinct and extant species in the
# communities of each grid square.
which(colSums(study_comm_2018) == 0) -> zeroes_in_2018
study_comm_1998[, -zeroes_in_2018] -> no_zero_comms_1998
study_comm_2018[, -zeroes_in_2018] -> no_zero_comms_2018

# Define a function to calculate FOri with the first 12 PCoA axes used in gamma
# scale FD calculations. See 3.1 and 3.3 for the quality of this dimensionality 
# reduction.
func_fori_12PCoA <- function(x){
  as.data.frame(t(x)) -> x
  originality(comm = x, tree = pcoa_axes[, 1:12], distance = gower_dists, 
              abund = F, relative = T)
}
# Apply the function.
as.data.frame(apply(X = no_zero_comms_1998, MARGIN = 2, 
                    FUN = func_fori_12PCoA)) -> fori_1998_comms

# Locate each extinct species' position.
which(no_zero_comms_1998 == 1 & no_zero_comms_2018 == 0) -> extinct_pos
matrix(ncol = ncol(fori_1998_comms), 
       nrow = nrow(fori_1998_comms)) -> extinct_matrix_pos
1 -> extinct_matrix_pos[extinct_pos]
colnames(fori_1998_comms) -> colnames(extinct_matrix_pos)
rownames(study_comm_1998) -> rownames(extinct_matrix_pos)
as.data.frame(extinct_matrix_pos) -> extinct_matrix_pos

# Get the mean ranks of the extinct and extant species in each community, and
# learn whether the functionally the most original species in each of the grid
# squares have gone extinct or not. 
for(i in 1:ncol(fori_1998_comms)){
  colnames(extinct_matrix_pos)[i] -> square_name
  which(extinct_matrix_pos[, i] == 1) -> extinct_indexer
  rank(fori_1998_comms[, i], na.last = "keep")[extinct_indexer] -> extinct_ranks
  rank(fori_1998_comms[, i], 
     na.last = NA) %ni% extinct_ranks -> extant_rank_indexer
  rank(fori_1998_comms[, i], 
     na.last = NA)[extant_rank_indexer] -> extant_ranks
  mean(extinct_ranks, na.rm = T) -> mean_extinct_rank
  mean(extant_ranks, na.rm = T) -> mean_extant_rank
  # Mean ranks of the extinct and extant species in each square.
  min(rank(fori_1998_comms[, i], na.last = NA), na.rm = T) -> min_rank
  # The highest rank in square.
  if(min_rank %in% extinct_ranks == T) {
  "The most functionally original species is among the extinct ones" -> max_rslt
    } else {
  "The most functionally original species lives" -> max_rslt
      }
  print(paste("Mean rank of extinct species in", square_name, "is", 
              round(mean_extinct_rank, 2), sep = " "))
  print(paste("Mean rank of extant species in", square_name, "is", 
              round(mean_extant_rank, 2), sep = " "))
  print(max_rslt)
}
# We can see that 23 out of 28 have lost the most original species. In addition, 
# extinct species had on average 15% higher ranks (in terms of their FOri) when 
# compared to the extant species.




