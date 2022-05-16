# 1: DATA PREPARATION ####


# R version used: 4.0.3.

# Load the needed libraries. 
library(tidyverse) # Version used: 1.3.1
library(glmmTMB) # Version used: 1.1.3
library(DHARMa) # Version used: 0.4.5
library(MASS) # Version used: 7.3-54
library(MuMIn) # Version used: 1.46.0
library(betapart) # Version used: 1.5.6
library(ggpubr) # Version used: 0.4.0
library(performance) # Version used: 0.9.0

# Read in the source files. The files given here include A-coded breeding 
# records (i.e., possible records) too. For more information on the breeding 
# codes see the article's supplementary material.
read.csv2(file = "taxonomic diversity/data/1998_breeding_all_species.csv", 
          header = T) -> all_species_1998
read.csv2(file = "taxonomic diversity/data/2018_breeding_all_species.csv", 
          header = T) -> all_species_2018

# Ready the data frames. These data frames have all the species wetland species 
# that have ever been detected in the study area. Therefore, we need to filter 
# out those that do not occur in our dataset.
all_species_1998[,1] -> rownames(all_species_1998)
all_species_2018[,1] -> rownames(all_species_2018)
all_species_1998[,1] <- NULL ; all_species_2018[,1] <- NULL
all_species_1998[rowSums(all_species_1998)!=0, ] -> study_comm_1998
all_species_2018[rowSums(all_species_2018)!=0, ] -> study_comm_2018
# Check the row numbers of these data frames. There were 120 breeding birds in
# the study area, Konya Closed Basin (KCB) in 1998 and 97 in 2018.

# Define a function that's reverse of %in% for the purpose of getting the names 
# of the species that have gone extinct and have colonized the KCB. 
Negate("%in%") -> "%ni%"
rownames(study_comm_1998)[(rownames(study_comm_1998) %ni% 
rownames(study_comm_2018))] -> present_in_1998_absent_in_2018
rownames(study_comm_2018)[(rownames(study_comm_2018) %ni% 
rownames(study_comm_1998))] -> present_in_2018_absent_in_1998
present_in_1998_absent_in_2018 # The 27 species that have gone extinct.
present_in_2018_absent_in_1998 # The 4 colonizers/recolonizers, which means we
# have a total of 124 species in our study. 

# The following chunk gets the names of the all 124 species in both data frames
# as row names. 
rownames(study_comm_1998[c(121, 122, 123, 124), ]) <- 
  c("Cattle egret", "Great cormorant", "White-tailed lapwing", 
    "Yellow-legged gull")
rep(0, 34) -> study_comm_1998[c(121, 122, 123, 124), ]  # Populating the newly 
# added rows with zeroes (they are absent in 1998).
study_comm_1998[order(row.names(study_comm_1998)), ] -> study_comm_1998
present_in_1998_absent_in_2018 ->
  rownames(study_comm_2018[c(seq(from = 98, to = 124, by = 1)), ])
rep(0, 34) -> study_comm_2018[c(seq(from = 98, to = 124, by = 1)), ]
study_comm_2018[order(row.names(study_comm_2018)), ] -> study_comm_2018
rownames(study_comm_1998) == rownames(study_comm_2018) # Final check to see if 
# the species list is the same in both data frames.
remove(all_species_1998, all_species_2018)

# Transpose the data frames for the use in beta.multi function (see below).
as.data.frame(t(study_comm_1998)) -> study_comm_1998
as.data.frame(t(study_comm_2018)) -> study_comm_2018


# 2: TAXONOMIC ALPHA & GAMMA DIVERSITY: SPECIES RICHNESS ####


# The decline in gamma-scale species richness was 120-97=23 species, which is a
23/120*100 # 19.2% decline.

# Calculate the species richness for the squares to get the alpha-scale richness
# values.
colSums(study_comm_1998) -> richness_1998
colSums(study_comm_2018) -> richness_2018
c(richness_1998, richness_2018) -> richness_both 
# Unite the data in a single data frame that is in a "tidy" format.
data.frame(year = c(rep("1998",34), rep("2018", 34)), 
           square = rep(colnames(study_comm_1998), 2),
           richness = richness_both, 
           index = rep("Species richness", 68)) -> alpha_richness_df
# Convert the year column into a factor to ready the data frame for 
# visualization with ggplot and for statistical modeling (see below).
as.factor(alpha_richness_df$year) -> alpha_richness_df$year
# Sort the data frame according to the square names and years.
tibble(alpha_richness_df) -> alpha_richness_df
alpha_richness_df %>% 
  arrange(year, square) -> alpha_richness_df

# Also read in the area data and relative coordinates.
read.csv2(file = "taxonomic diversity/data/areas.csv", 
          header = T) -> areas
rep(as.numeric(areas$area), 2) -> alpha_richness_df$area
read.csv2("taxonomic diversity/data/relative_coords.csv", 
          header = T) -> rel_coords
rel_coords$square[1:34] == alpha_richness_df$square[1:34]
rep(rel_coords$x, 2) -> alpha_richness_df$x
rep(rel_coords$y, 2) -> alpha_richness_df$y
# Divide by 50,000 to transform the x and y coordinates.
alpha_richness_df$x / 50000 -> alpha_richness_df$xdiv
alpha_richness_df$y / 50000 -> alpha_richness_df$ydiv

# Transpose the community matrices.
as.data.frame(t(study_comm_1998)) -> study_comm_1998
as.data.frame(t(study_comm_2018)) -> study_comm_2018
  
## Data exploration ##
# Now explore the data. Our task is to see whether the study years differ in 
# terms of their species richness. Other possible covariates are area and the
# spatial coordinates.

# Check for zero inflation in both data frames.
which(rowSums(study_comm_1998)==0) ; which(rowSums(study_comm_2018)==0)
# No square has zero breeding species in 1998's atlas data, but 6 squares have 
# zero breeding species in 2018's atlas data. It is unlikely that we have a zero
# inflation problem. Still, this will be checked with DHARMa's utilities.

# Check for the outliers.
dotchart(alpha_richness_df$richness, 
         labels = alpha_richness_df$square, bg = "red")
# There are some high values, but they should be kept as there is nothing wrong 
# with these data points. We can check the Cook distances post-hoc to see if 
# there are any influential observations.
dotchart(alpha_richness_df$area, 
         labels = alpha_richness_df$square, bg = "green")
# No outliers.

# Plot the response against the possible covariates and random effects. Start
# with the year. Also check for any possible interactions and homogeneity in 
# variances.
ggplot() +
  geom_boxplot(data = alpha_richness_df, 
               aes(x = index, y = richness, fill = year), 
               alpha = 0.8, notch = T) +   ylab(label = "Species Richness") +
  scale_fill_manual(values = c("deepskyblue", "deepskyblue4")) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x.bottom = element_blank()) +
  guides(fill = guide_legend("Year")) -> p_alpha_richness
p_alpha_richness # There seems to be a year effect, there's a slight variance
                 # but it is not too serious. 
# Save the plot (optional).
# ggsave(filename = "taxonomic diversity/figs/p_alpha_richness.tiff", 
#       device = "tiff", dpi = 500, height = 3, width = 2.9)

# Plot species richness vs. the area.
ggplot(data = alpha_richness_df) +
  geom_point(aes(x = area, y = richness)) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  xlab("Area") + ylab("Richness") + guides(col = guide_legend("Year")) +
  geom_smooth(aes(x = area, y = richness)) # Pooled data.
ggplot(data = alpha_richness_df) +
  geom_point(aes(x = area, y = richness, col = year)) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  xlab("Area") + ylab("Richness") + guides(col = guide_legend("Year")) +
  geom_smooth(aes(x = area, y = richness, col = year)) # Separate years.
# Obviously, there is an area effect, but there don't seem to be an interaction 
# with the year covariate.

# Plot richness vs. square IDs for each year. 
par(mfrow = c(2, 1)); par(mar = c(2, 2, 2, 2))
boxplot(alpha_richness_df$richness[1:34] ~ alpha_richness_df$square[1:34], 
        xlab = "Squares", ylab = "Richness")
boxplot(alpha_richness_df$richness[35:68] ~ alpha_richness_df$square[35:68], 
        xlab = "Squares", ylab = "Richness")
par(mfrow = c(1, 1))
# Huge variation. Our test should be a paired test anyway, so this variable will
# be defined as a random variable since we are not after its effect.

# Plot the response against x and y coordinates to look for spatial patterns.
ggplot() +
  geom_point(data = alpha_richness_df[1:34, ], aes(x = x, y = y, col =richness), 
             size = 5) +
  scale_color_continuous(type = "gradient") +
  guides(col = guide_legend("Richness")) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  xlab("Easting") + ylab("Northing") -> p_check1 # For 1998's data.
ggplot() +
  geom_point(data = alpha_richness_df[35:68, ], aes(x = x, y = y, col = richness), 
             size = 5) +
  scale_color_continuous(type = "gradient") +
  guides(col = guide_legend("Richness")) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  xlab("Easting") + ylab("Northing") -> p_check2  # For 2018's data
ggarrange(p_check1, p_check2, common.legend = T, nrow = 1, ncol = 2)
# The smaller, marginal squares seem to have lower species richness as expected.
# Keep this in mind. This suggest we may consider including area or coordinates
# in our model.

# Summary of the data exploration:
# (i)   No zero-inflation, no outliers.
# (ii)  There seems to be a year effect
# (iii) There is a relationship between area and species richness.
# (iv)  Lots of variation among richness in squares.
# (v)   Marginal squares have smaller areas and lower species richness.

# Now, compare the species richness values between the years by using 
# generalized linear mixed effect models (GLMMs).
# Let's scale the area covariate before we start.
scale(alpha_richness_df$area) -> alpha_richness_df$area
as.factor(alpha_richness_df$square) -> alpha_richness_df$square
# Build the model. Define area and square as random effects and see if there is
# there is spatial autocorrelation. If that's the case try to include spatial
# coordinates in the model or switch to a GAM, simply.
m_alpha_richness <- glmmTMB(data = alpha_richness_df, na.action = "na.fail",
                            formula = richness ~ year + (1|area) + (1|square), 
                            family = nbinom1(link = "log"))
summary(m_alpha_richness) # The year effect seems to be significant. Let's 
                          # check the model.

## Model validation ##
# 1) Plot residuals vs. the fitted values.
par(mfrow = c(1,1))
residuals(m_alpha_richness) -> resids
fitted(m_alpha_richness) -> fits
plot(fits, resids, pch = 16); abline(h = 0, col = "red", lty = 2)
# This looks okayish.
# 2) Plot the residuals vs. the predictor.
boxplot(resids ~ alpha_richness_df$year); abline(h = 0, col = "red", lty = 2)
# They look good.
# 3) Check normality of the random effects.
check_model(m_alpha_richness, check = "reqq")
# No problem here.
# 4) Plot the residuals against spatial coordinates to look for spatial patterns
# and also test to spatial autocorrelation.
as.data.frame(resids) -> resids
as.data.frame(resids[order(as.numeric(rownames(resids))), ]) -> resids
resids$`resids[order(as.numeric(rownames(resids))), ]` -> 
  alpha_richness_df$resids
ggplot(data = alpha_richness_df[1:34, ]) +
  geom_point(aes(x = x, y = y, fill = resids), 
             col = "black", pch = 21, size = 6) +  
  scale_fill_gradient2(low = "green", mid = "yellow", high = "red") +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  theme(legend.position = "None") +
  xlab("Easting") + ylab("Northing") # This should be checked (1998).
ggplot(data = alpha_richness_df[35:68, ]) +
  geom_point(aes(x = x, y = y, fill = resids), 
             col = "black", pch = 21, size = 6) +  
  scale_fill_gradient2(low = "green", mid = "yellow", high = "red") +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  theme(legend.position = "None") +
  xlab("Easting") + ylab("Northing") # This should be checked too (2018).
# Test for spatial autocorrelation now.
simulateResiduals(m_alpha_richness, 1000) -> sim_alpha_richness
set.seed(13); recalculateResiduals(sim_alpha_richness, 
                     group = alpha_richness_df$square) -> sim_alpha_richness_2
set.seed(13); testSpatialAutocorrelation(sim_alpha_richness_2, 
                           x = alpha_richness_df[1:34,]$x,
                           y = alpha_richness_df[1:34,]$y)
# No significant autocorrelation, which means there is no need to include the
# coordinates in the model.
# 5) Check the acf and pacf plots.
par(mfrow = c(1, 2))
acf(resids); pacf(resids) # There are random spikes in both plots.
par(mfrow = c(1, 1))
# There are no patterns in acf and pacf plots.
# 6) Use the DHARMa utilities to check the model for overdispersion and 
# outliers.
testDispersion(sim_alpha_richness)
testZeroInflation(sim_alpha_richness)
testOutliers(sim_alpha_richness)
# No problem here.
# 7) Check whether the model performs better than the null model.
dredge(m_alpha_richness)
# The model is better than the null model.
# The model is valid and the year 1998 has significantly higher species 
# richness (p = 4.51e-11***).


# 3: TAXONOMIC BETA DIVERSITY ####


# Calculate pairwise beta diversity for each pair of squares for 1998 and 2018.
# Order the community matrices by their row names. This is done for the distance
# decay relationship that we will try to obtain (see below).
study_comm_1998[order(rownames(study_comm_1998)), ] -> study_comm_1998
study_comm_2018[order(rownames(study_comm_2018)), ] -> study_comm_2018
beta.pair(study_comm_1998, index.family = "sorensen") -> pairwise_1998
beta.pair(study_comm_2018, index.family = "sorensen") -> pairwise_2018


## 3.1: DISTANCE-DECAY MODELS ####


# Get the distances.
dist(x = alpha_richness_df[1:34, 6:7], method = "euclidean") -> utm_distances
as.matrix(utm_distances) -> utm_distances # Just to check.
as.dist(m = utm_distances, diag = F, upper = T) -> utm_distances

# Build the distance decay models now. For that purpose, define a function that
# builds and selects distance decay models (ddms).
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
              if(exp_models[[i]]$model$aic - power_models[[i]]$model$aic > 2){
                mod_type <- "power"
                print(paste("p =",power_models[[i]]$p.value, "for",index_type,
                            "model, type =",mod_type))
                selected_models[[i]] <- power_models[[i]]
              } else {
                mod_type <- "exponential"
                print(paste("p =",exp_models[[i]]$p.value, "for",index_type,
                            "model, type =", mod_type))
                selected_models[[i]] <- exp_models[[i]]
                # This chunk builds and selects between the power and 
                # exponential models and stores them in selected_models object.
                # The selection is AIC-based. 
              }
    }
    return(selected_models)
}

# Run this function on 1998's and 2018's pairwise beta diversity matrices.
set.seed(13); func_decay_select(pairwise_1998) -> ddms_1998
set.seed(13); func_decay_select(pairwise_2018) -> ddms_2018
# All models, except for the distance decay model (ddm) for Bsim and Bsne in 
# 2018, had significant slopes. Check the validity of these models now.

# Model validation #
# 1) Check the normality and distribution of the residuals.
boxplot(ddms_1998[[3]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# They look pretty good. No problem here.
# 2) Plot the residuals against the only predictor in the model: distances to 
# see if there are any patterns in the plot. 
plot(utm_distances,
      ddms_1998[[3]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# There seems to be no significant patterns here.
# 3) Check fitted vs. residuals graph.
plot(ddms_1998[[3]]$model$fitted.values, 
     ddms_1998[[3]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# The model is valid.
# 1) Check the distribution of the residuals around zero.
boxplot(ddms_2018[[1]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# They look pretty good.
# 2) Plot the residuals against the only predictor in the model: distances to 
# see if there are any patterns in the plot. 
plot(utm_distances,
     ddms_2018[[1]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# 3) Plot the residuals against the fitted values, as it is easier than 
# filtering the distances matrix.
plot(ddms_2018[[1]]$model$fitted.values, 
      ddms_2018[[1]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# The pattern here is because of the nature of the covariate (distance), so 
# there seems to be no significant patterns here. The model is valid.
# 1) Check the distribution of the residuals around zero.
boxplot(exp_mod_Bsor_2018$model$residuals); abline(h = 0, col = "red", lty = 2)
# The distribution of the residuals is suboptimal.
# 2) Plot the residuals against the fitted values, as it is easier than 
# filtering the distances matrix.
plot(ddms_2018[[3]]$model$fitted.values, 
      ddms_2018[[3]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# There seems to be no significant patterns here. The residuals' distribution is
# not perfect, but we conclude that is is acceptable and the model is valid.

# Visualize the distance decay models now.
par(mfrow = c(1, 3))
lapply(ddms_1998, FUN = plot.decay); lapply(ddms_2018, FUN = plot.decay)
par(mfrow = c(1, 1))

# Define two transparent colors to be used in the plots.
adjustcolor(col = "deepskyblue", alpha.f = 0.9) -> blue_tr
adjustcolor(col = "deepskyblue4", alpha.f = 0.9) -> dark_blue_tr

# Calculate bootstrapped parameters for these models now compare the slopes of 
# the models between the study years. Start with the species turnover component.
set.seed(13); boot.coefs.decay(ddms_1998[[1]], 
                               resamples = 10000) -> boot_ddm_1998_sim
set.seed(13); boot.coefs.decay(ddms_2018[[1]], 
                               resamples = 10000) -> boot_ddm_2018_sim


### 3.1.1: COMPARISON OF THE SLOPES ####


# Plot the model parameters side-by-side.
data.frame(slope = c(boot_ddm_1998_sim$boot.coefs[, 2], 
                     boot_ddm_2018_sim$boot.coefs[, 2]), 
           year = as.factor(c(rep("1998", 10000), 
                              rep("2018", 10000)))) -> Bsim_slopes
ggplot() +
  geom_boxplot(data = Bsim_slopes, 
               aes(x = year, y = slope, fill = year), alpha = 0.8) +
  scale_fill_manual(values = c("deepskyblue", "deepskyblue4")) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  theme(axis.title.x = element_text(size = 12)) + 
  labs(fill = "Year") +
  ylab("Bootstrapped DDM slope for Bsim") + xlab("Year") -> p_ddm_slope_Bsim
p_ddm_slope_Bsim
# Save the plot (optional).
# ggsave(filename = "taxonomic diversity/figs/p_ddm_slope_Bsim.tiff", 
#       device = "tiff", dpi = 600, height = 4, width = 3.4)

# Make the formal comparison by calculating p-value as the proportion of times
# one of the slopes is larger than the other.
Bsim_slopes$slope[1:10000] - Bsim_slopes$slope[10001:20000] -> slope_dif_Bsim
length(which(slope_dif_Bsim < 0))
length(which(slope_dif_Bsim > 0)) # Smaller one.
(length(which(slope_dif_Bsim > 0)) / 10000)*2
# 2018's communities have steeper slope for the decay of Bsim (p = 0.0146*).

# Comparisons for Bsne models won't be made as the best Bsne model for 2018 does
# not have a significant slope.

# Now compare the slopes of the Bsor models. 
set.seed(13)
boot.coefs.decay(ddms_1998[[3]], resamples = 10000) -> boot_ddm_1998_sor
set.seed(13)
boot.coefs.decay(ddms_2018[[3]], resamples = 10000) -> boot_ddm_2018_sor

# Plot the model parameters side-by-side.
range(c(boot_ddm_1998_sor$boot.coefs[,2], 
        boot_ddm_2018_sor$boot.coefs[,2])) -> current_range
boxplot(notch = T, boot_ddm_2018_sor$boot.coefs[,2], col = dark_blue_tr,
        ylim = current_range,
        main = "Slope Bsor Bootstrapped - Blue:1998, Dark blue:2018")
boxplot(notch = T, boot_ddm_1998_sor$boot.coefs[,2], col = blue_tr,
        ylim = current_range, add = T)
# The slopes look quite different.
data.frame(slope = c(boot_ddm_1998_sor$boot.coefs[, 2], 
                     boot_ddm_2018_sor$boot.coefs[, 2]), 
           year = as.factor(c(rep("1998", 10000), 
                              rep("2018", 10000)))) -> Bsor_slopes
ggplot() +
  geom_boxplot(data = Bsor_slopes, 
               aes(x = year, y = slope, fill = year), alpha = 0.8) +
  scale_fill_manual(values = c("deepskyblue", "deepskyblue4")) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  theme(axis.title.x = element_text(size = 12)) + 
  labs(fill = "Year") +
  ylab("Bootstrapped DDM slope for BSor") + xlab("Year") -> p_ddm_slope_Bsor
p_ddm_slope_Bsor
# Save the plot (optional).
# ggsave(filename = "taxonomic diversity/figs/p_ddm_slope_Bsor.tiff", 
#       device = "tiff", dpi = 600, height = 4, width = 3.5)

# Make the formal comparison.
Bsor_slopes$slope[1:10000] - Bsor_slopes$slope[10001:20000] -> slope_dif_Bsor
length(which(slope_dif_Bsor < 0))
length(which(slope_dif_Bsor > 0)) # Smaller one.
(length(which(slope_dif_Bsor > 0)) / 10000)*2
# p-value for this comparison is zero (p = 0) and and the ddm of 2018 has a 
# significantly steeper slope.

# Produce a graph to visualize all the models side-by-side (optional).
# tiff(filename = "taxonomic diversity/figs/taxonomic_ddms.tiff", 
#      width = 240, height = 100, units = "mm", res = 500)
# par(mfrow = c(1, 3))
# for(i in 1:3){
#   plot.decay(ddms_1998[[i]], remove.dots = T, 
#             col = "deepskyblue", ylim = c(0,1), lwd = 3.5)
#    plot.decay(ddms_2018[[i]], remove.dots = T, col = "deepskyblue4", 
#               ylim = c(0,1), lwd = 3.5, add = T)
# }
# dev.off()



### 3.1.2: COMPARISON OF THE INTERCEPTS ####


# Check the intercepts now, but bootstrap the coefficients for the nestedness
# models as well.
set.seed(13)
boot.coefs.decay(ddms_1998[[2]], resamples = 10000) -> boot_ddm_1998_sne
set.seed(13)
boot.coefs.decay(ddms_2018[[2]], resamples = 10000) -> boot_ddm_2018_sne
data.frame(intercept = c(boot_ddm_1998_sim$boot.coefs[, 1], 
                     boot_ddm_2018_sim$boot.coefs[, 1]), 
           year = as.factor(c(rep("1998", 10000), 
                              rep("2018", 10000)))) -> Bsim_intercepts
data.frame(intercept = c(boot_ddm_1998_sne$boot.coefs[, 1], 
                         boot_ddm_2018_sne$boot.coefs[, 1]), 
           year = as.factor(c(rep("1998", 10000), 
                              rep("2018", 10000)))) -> Bsne_intercepts
data.frame(intercept = c(boot_ddm_1998_sor$boot.coefs[, 1], 
                     boot_ddm_2018_sor$boot.coefs[, 1]), 
           year = as.factor(c(rep("1998", 10000), 
                              rep("2018", 10000)))) -> Bsor_intercepts
# Check the box plots.
ggplot() +
  geom_boxplot(data = Bsim_intercepts, notch = T,
               aes(x = year, y = intercept, fill = year), alpha = 0.8) +
  scale_fill_manual(values = c("deepskyblue", "deepskyblue4")) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  theme(axis.title.x = element_text(size = 12)) + 
  ylab("Bootstrapped DDM intercept for Bsim") + xlab("Year") +
  guides(fill = guide_legend("Year")) -> p_ddm_intercept_Bsim
p_ddm_intercept_Bsim # The intercepts look significant and different.
# Save the plot (optional).
# ggsave(filename = "taxonomic diversity/figs/p_ddm_intercept_Bsim.tiff", 
#       device = "tiff", dpi = 600, height = 4, width = 3.5)
ggplot() +
  geom_boxplot(data = Bsne_intercepts, notch = T,
               aes(x = year, y = intercept, fill = year), alpha = 0.8) +
  scale_fill_manual(values = c("deepskyblue", "deepskyblue4")) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  theme(axis.title.x = element_text(size = 12)) + 
  ylab("Bootstrapped DDM intercept for Bsne") + xlab("Year") +
  guides(fill = guide_legend("Year")) -> p_ddm_intercept_Bsne
p_ddm_intercept_Bsne # They don't look too different.
ggplot() +
  geom_boxplot(data = Bsor_intercepts, notch = T, 
               aes(x = year, y = intercept, fill = year), alpha = 0.8) +
  scale_fill_manual(values = c("deepskyblue", "deepskyblue4")) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  theme(axis.title.x = element_text(size = 12)) + 
  ylab("Bootstrapped DDM intercept for Bsor") + xlab("Year") +
  guides(fill = guide_legend("Year")) -> p_ddm_intercept_Bsor
p_ddm_intercept_Bsor # They look different and significant.
# Save the plot (optional).
# ggsave(filename = "taxonomic diversity/figs/p_ddm_intercept_Bsor.tiff", 
#       device = "tiff", dpi = 600, height = 4, width = 3.5)

# Calculate the p-values now.
Bsim_intercepts$intercept[1:10000] - Bsim_intercepts$intercept[10001:20000] -> 
  intercept_dif_Bsim
length(which(intercept_dif_Bsim < 0)) 
length(which(intercept_dif_Bsim > 0)) # Smaller one.
(length(which(intercept_dif_Bsim > 0)) / 10000)*2
# The difference in intercepts is significant (p = 6e-04***) with 2018's 
# communities having a higher intercept for Bsim ddm.

Bsne_intercepts$intercept[1:10000] - Bsne_intercepts$intercept[10001:20000] -> 
  intercept_dif_Bsne
length(which(intercept_dif_Bsne < 0)) 
length(which(intercept_dif_Bsne > 0)) # Smaller one.
(length(which(intercept_dif_Bsne > 0)) / 10000)*2
# The difference in intercepts is not significant (p = 0.669).
Bsor_intercepts$intercept[1:10000] - Bsor_intercepts$intercept[10001:20000] -> 
  intercept_dif_Bsor
length(which(intercept_dif_Bsor < 0)) 
length(which(intercept_dif_Bsor > 0)) # Smaller one.
(length(which(intercept_dif_Bsor > 0)) / 10000)*2
# The difference in intercepts is significant (p = 0) with 2018's communities 
# having a higher intercept for Bsor ddm.


## 3.2: COMPONENTS OF BETA DIVERSITY ####


# Create a data frame that will have all the data.
data.frame(index = c(rep("Bsim", 1122), rep("Bsne", 1122), rep("Bsor", 1122)), 
           year = as.factor(rep(c(rep("1998", 561), rep("2018", 561)), 3)),
           identifier = as.factor(rep(rep(paste("X", 1:561, sep = ""), 2), 3)),
           value = c(as.vector(pairwise_1998$beta.sim),
                     as.vector(pairwise_2018$beta.sim),
                     as.vector(pairwise_1998$beta.sne),
                     as.vector(pairwise_2018$beta.sne),
                     as.vector(pairwise_1998$beta.sor),
                     as.vector(pairwise_2018$beta.sor))) -> beta_data

# We will use Wilcoxon signed-rank test, which is a non-parametric test, to 
# compare the dissimilarities in the two matrices.
wilcox.test(x = as.vector(pairwise_1998$beta.sim),  
            y = as.vector(pairwise_2018$beta.sim),
            alternative = "two.sided")
# 2018's communities have significantly higher Bsim values (p < 2.2e-16***).
wilcox.test(x = as.vector(pairwise_1998$beta.sne),  
            y = as.vector(pairwise_2018$beta.sne),
            alternative = "two.sided")
# 1998's communities have significantly higher Bsne values (p = 5.075e-08***).
wilcox.test(x = as.vector(pairwise_1998$beta.sor),  
            y = as.vector(pairwise_2018$beta.sor),
            alternative = "two.sided")
# 2018's communities have significantly higher Bsor values (p < 2.2e-16***).

# Create a plot for the taxonomic beta diversity components.
ggplot() +
  geom_boxplot(data = beta_data, aes(x = index, y = value, fill = year), 
               alpha = 0.8, notch = T) +
  scale_fill_manual(values = c("deepskyblue", "deepskyblue4")) +
  theme(panel.border = element_rect(colour = "gray35", fill = NA, size = 0.5)) +
  ylab(label = "Index value") + xlab(label = "Index") +
  theme(axis.title.x = element_text(size = 12)) +
  guides(fill = guide_legend("Year")) -> p_beta_taxonomic
p_beta_taxonomic
# Save the plot.
# ggsave(filename = "taxonomic diversity/figs/p_taxonomic_beta.tiff", dpi = 400, 
#        units = "mm", width = 120, height = 75)




