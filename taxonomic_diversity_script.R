# 1: DATA PREPARATION ####


# R version used: 4.0.3.

# Load the needed libraries. 
library(ggplot2) # Version used: 3.3.5
library(glmmTMB) # Version used: 1.1.2.2
library(DHARMa) # Version used: 0.4.3
library(MASS) # Version used: 7.3-53
library(MuMIn) # Version used: 1.43.17
library(betapart) # Version used: 1.5.4

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

# Define a function that's reverse of %in% with the ultimate goal of getting
# the names of the species that have gone extinct and have colonized the KCB. 
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
rowSums(study_comm_1998) -> richness_1998
rowSums(study_comm_2018) -> richness_2018
c(richness_1998, richness_2018) -> richness_both 
# Unite the data in a single data frame that is in "tidy" format.
data.frame(year = c(rep("1998",34), rep("2018", 34)), 
           square = rep(rownames(study_comm_1998), 2),
           richness = richness_both, 
           index = rep("Species richness", 68)) -> alpha_richness_df
# Convert the year column into a factor to ready the data frame for 
# visualization with ggplot and for statistical modeling (see below).
as.factor(alpha_richness_df$year) -> alpha_richness_df$year

# Plot a box plot to visualize the alpha-scale richness data for both years.
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
p_alpha_richness
# Save the plot (optional).
 ggsave(filename = "taxonomic diversity/figs/p_alpha_richness.tiff", 
       device = "tiff", dpi = 500, height = 3, width = 2.9)

# Check if any of the squares have zero species richness in both data frames.
which(rowSums(study_comm_1998)==0) ; which(rowSums(study_comm_2018)==0)
# No square has zero breeding species in 1998's atlas data, but 6 squares have 
# zero breeding species in 2018's atlas data.

# Now, compare the species richness values between the years by using 
# generalized linear models (GLMs).
m_alpha_richness <- glmmTMB(data = alpha_richness_df, formula = richness ~ year, 
                        family = nbinom1(link = "log"))
summary(m_alpha_richness) # Year 2018 seems to have significantly lower species 
# richness (p = 9.97e-06***). Check the model assumptions and validity.
# Start checking the model.

# Model validation #
# 1) Plot the model's residuals to check their distribution around zero.
boxplot(residuals(m_alpha_richness)) ; abline(h = 0, col = "red", lty = 2)
# The distribution is not perfectly symmetrical around zero, but it seems 
# acceptable.
# 2) Check the model's performance at different levels of the explanatory 
# variable, year, by plotting the residuals against the year variable.
boxplot(residuals(m_alpha_richness) ~ alpha_richness_df$year)
abline(h = 0, col = "red", lty =2)
# The model seems to perform a marginally worse at predicting for the 1998 level
# but the difference is at acceptable levels. 
# 3) Plot the residuals against the square ID variable, which was not included
# in the model.
boxplot(residuals(m_alpha_richness) ~ alpha_richness_df$square)
abline(h = 0, col = "red", lty =2)
# The model's performance at the levels of the square ID variable is quite
# variable. We can now include this variable in the model as a random variable 
# to control for the effect of the square ID.
as.factor(alpha_richness_df$square) -> alpha_richness_df$square
m_alpha_richness_2 <- glmmTMB(data = alpha_richness_df, 
                            formula = richness ~ year + (1 |square), 
                        family = nbinom1(link = "log"))
summary(m_alpha_richness_2) # The effect persists (p = 7.68e-11***).
AICc(m_alpha_richness, m_alpha_richness_2) # The GLMM has a way lower AICc, 
# which means this model is the better model. Let's check the new model.
# 1) Plot the model's residuals to check their distribution around zero.
boxplot(residuals(m_alpha_richness_2)) ; abline(h = 0, col = "red", lty = 2)
# The distribution is not perfectly symmetrical around zero, but it seems 
# barely acceptable.
# 2) Check the model's performance at different levels of the explanatory 
# variable, year, by plotting the residuals against the year variable.
boxplot(residuals(m_alpha_richness_2) ~ alpha_richness_df$year)
abline(h = 0, col = "red", lty =2)
# The model still  seems to perform a marginally worse at predicting for the 
# 1998 level, but we think the difference is at acceptable levels. 
# 3) Plot the residuals against the square ID variable, which is now included
# in the model.
boxplot(residuals(m_alpha_richness_2) ~ alpha_richness_df$square)
abline(h = 0, col = "red", lty =2)
# The model's performance at the levels of the square ID variable looks way 
# better and closer to the h=0 line, but it is not perfect. We will assume it is
# acceptable.
# 4) Test the model for spatial autocorrelation as the squares that are close to
# each other might have similar richness values and richness changes. To do that
# add the coordinates (relative UTM coordinates in our case) to the data frame.
read.csv2("taxonomic diversity/data/relative_coords.csv", 
          header = T) -> rel_coords
alpha_richness_df[order(alpha_richness_df$year,
                        alpha_richness_df$square), ] -> alpha_richness_df
rel_coords$square[1:34] == alpha_richness_df$square[1:34]
rep(rel_coords$x, 2) -> alpha_richness_df$x
rep(rel_coords$y, 2) -> alpha_richness_df$y
set.seed(13)
simulateResiduals(m_alpha_richness_2, 1000) -> sim_alpha_richness_2
set.seed(13); recalculateResiduals(sim_alpha_richness_2, 
                     group = alpha_richness_df$square) -> sim_alpha_richness_2
set.seed(13); testSpatialAutocorrelation(sim_alpha_richness_2, 
                           x = alpha_richness_df[1:34,]$x,
                           y = alpha_richness_df[1:34,]$y)
# The residuals have no spatial autocorrelation.
# 5) Check the model residuals for autocorrelation and spatial autocorrelation.
par(mfrow = c(1, 2))
acf(residuals(m_alpha_richness_2)); pacf(residuals(m_alpha_richness_2))
par(mfrow = c(1, 1))
# The deviations are barely significant seem random, so we will assume that 
# there is no autocorrelation problem in the residuals. 
# 6) Use the utilities of the DHARMa package to test the residuals for 
# uniformity and over/underdispersion.
set.seed(13)
simulateResiduals(m_alpha_richness_2, 1000) -> sim_alpha_richness_2
set.seed(13); testResiduals(sim_alpha_richness_2)
# No problem in the results. 
# We can conclude that the model is a valid one and the that the 2018's 
# communities have significantly lower species richness (p = 7.68e-11***). 


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
dist(x = alpha_richness_df[1:34, 5:6], method = "euclidean") -> utm_distances
as.matrix(utm_distances) -> utm_distances # Just to check.
as.dist(m = utm_distances, diag = F, upper = T) -> utm_distances

# Some of the beta diversity values are zeros, which cannot be used with a log
# link function (see below) in distance decay modeling. Therefore, we need to
# transform the data.
func_log_trans <- function(x){
   (x * 32 + 0.5)/ 33
}
lapply(X = pairwise_1998, FUN = func_log_trans) -> pairwise_1998
lapply(X = pairwise_2018, FUN = func_log_trans) -> pairwise_2018

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
set.seed(13); func_decay_select(pairwise_1998) -> ddms_1998
set.seed(13); func_decay_select(pairwise_2018) -> ddms_2018
# All models, except for the distance decay model (ddm) for Bsim and Bsne in 
# 2018, had significant slopes. Check the validity of these models now.

# Model validation #
# 1) Check the distribution of the residuals around zero.
boxplot(ddms_1998[[1]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# They look pretty good. No problem here.
# 2) Plot the residuals against the only predictor in the model: distances to 
# see if there are any patterns in the plot. 
plot(log(utm_distances), 
      ddms_1998[[1]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# There seems to be no significant patterns other than the model performing a
# little bit worse for the higher distances. The model looks valid. Repeat the
# procedure for the other five models.
# 1) Check the distribution of the residuals around zero.
boxplot(ddms_1998[[2]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# The distribution looks decent. 
# 2) Plot the residuals against the only predictor in the model: distances to 
# see if there are any patterns in the plot. 
plot(log(utm_distances), 
      ddms_1998[[2]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# There seems to be no significant patterns. The model looks good.
# 1) Check the distribution of the residuals around zero.
boxplot(ddms_1998[[3]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# They look pretty good. No problem here.
# 2) Plot the residuals against the only predictor in the model: distances to 
# see if there are any patterns in the plot. 
plot(log(utm_distances),
      ddms_1998[[3]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# There seems to be no significant patterns here. The model looks good.
# 1) Check the distribution of the residuals around zero.
boxplot(ddms_2018[[1]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# They look pretty good. No problem here.
# 2) Plot the residuals against the fitted values, as it is easier than 
# filtering the distances matrix.
plot(ddms_2018[[1]]$model$fitted.values, 
      ddms_2018[[1]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# There seems to be no significant patterns other than the model performing a
# little bit worse for lower fitted values. Still, the model looks okay.
# There seems to be no significant patterns here. The model looks good.
# 1) Check the distribution of the residuals around zero.
boxplot(ddms_2018[[2]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# The distribution of the residuals do not look perfect.
# 2) Plot the residuals against the fitted values, as it is easier than 
# filtering the distances matrix.
plot(ddms_2018[[2]]$model$fitted.values, 
      ddms_2018[[2]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# There seems to be no significant patterns. The model does not have a 
# significant slope, and the residuals are not perfect. There is no valid
# distance-decay relationship here.
# There seems to be no significant patterns here. The model looks good.
# 1) Check the distribution of the residuals around zero.
boxplot(ddms_2018[[3]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# The distribution of the residuals is not perfect. The model may not be valid.
# 2) Plot the residuals against the fitted values, as it is easier than 
# filtering the distances matrix.
plot(ddms_2018[[3]]$model$fitted.values, 
      ddms_2018[[3]]$model$residuals); abline(h = 0, col = "red", lty = 2)
# There seems to be no significant patterns here. The residuals' distribution is
# not perfect, but we conclude that is is acceptable.

# Visualize the distance decay models now.
par(mfrow = c(1, 3))
lapply(ddms_1998, FUN = plot.decay); lapply(ddms_2018, FUN = plot.decay)
par(mfrow = c(1, 1))

# Define two transparent colors to be used in the plots.
adjustcolor(col = "deepskyblue", alpha.f = 0.9) -> blue_tr
adjustcolor(col = "deepskyblue4", alpha.f = 0.9) -> dark_blue_tr

# Calculate bootstrapped parameters for these models now compare the slopes of 
# the models between the study years. Start with the species turnover component.
set.seed(13); boot.coefs.decay(ddms_1998[[1]], R = 10000) -> boot_ddm_1998_sim
set.seed(13); boot.coefs.decay(ddms_2018[[1]], R = 10000) -> boot_ddm_2018_sim

# Plot the model parameters side-by-side.
range(c(boot_ddm_1998_sim$boot.coefs[,2], 
        boot_ddm_2018_sim$boot.coefs[,2])) -> current_range
boxplot(notch = T, boot_ddm_2018_sim$boot.coefs[,2], col = dark_blue_tr,
        ylim = current_range,
        main = "Bsim Slope Bootstrapped - Blue:1998, Dark blue:2018")
boxplot(notch = T, boot_ddm_1998_sim$boot.coefs[,2], col = blue_tr,
        ylim = current_range, add = T)

# Produce a proper plot.
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
  ylab("Bootstrapped DDM slope") + xlab("Year") -> p_ddm_slope_Bsim
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
# 2018's communities have significantly steeper slopes (p = 0.007). 

# Repeat the procedure for the total beta diversity, but because the best Bsor 
# model for 2018 was the exponential model, we need to either get the 
# exponential model for 1998 or get the power model for 2018 to be able to 
# compare the slopes. Let's start with trying the exponential model for 1998.
set.seed(13)
decay.model(pairwise_1998$beta.sor, utm_distances, model.type = "exponential", 
            perm = 1000, y.type = "dissimilarities") -> exp_mod_Bsor_1998
# Check if the AIC scores are comparable (i.e., delta AIC<2).
ddms_1998[[3]]$model$aic; exp_mod_Bsor_1998$model$aic
# They are really close. The reason why the power model was chosen was because 
# the func_decay_select function automatically was defined to picks the power
# model unless the exponential model is significantly better. Since the two
# models are not significantly different, we can use the exponential model to 
# make the parameter comparison, but before that, let's check the model.

# Model validation #
# 1) Check the distribution of the residuals around zero.
boxplot(exp_mod_Bsor_1998$model$residuals); abline(h = 0, col = "red", lty = 2)
# They look nearly perfect.
# 2) Plot the residuals against the fitted values, as it is easier than 
# filtering the distances matrix.
plot(exp_mod_Bsor_1998$model$fitted.values, 
      exp_mod_Bsor_1998$model$residuals); abline(h = 0, col = "red", lty = 2)
# No patterns here. The model is valid and we can proceed.
set.seed(13)
boot.coefs.decay(exp_mod_Bsor_1998, R = 10000) -> boot_ddm_1998_sor
set.seed(13)
boot.coefs.decay(ddms_2018[[3]], R = 10000) -> boot_ddm_2018_sor

# Plot the model parameters side-by-side.
range(c(boot_ddm_1998_sor$boot.coefs[,2], 
        boot_ddm_2018_sor$boot.coefs[,2])) -> current_range
boxplot(notch = T, boot_ddm_2018_sor$boot.coefs[,2], col = dark_blue_tr,
        ylim = current_range,
        main = "Slope Bsor Bootstrapped - Blue:1998, Dark blue:2018")
boxplot(notch = T, boot_ddm_1998_sor$boot.coefs[,2], col = blue_tr,
        ylim = current_range, add = T)
# The slopes look quite different.

# Produce a graph to be used in the manuscript.
# ddms_1998[[3]] -> spare
# exp_mod_Bsor_1998 -> ddms_1998[[3]]
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
# spare -> ddms_1998[[3]]

# Now, let's compare the slopes of the Bsor models. Put the data in the "tidy"
# format to use it in a proper box plot to be included in the supplementary
# material.

# Repeat the same for the Bsor models.
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
  ylab("Bootstrapped DDM slope") + xlab("Year") +
  guides(fill = guide_legend("Year")) -> p_ddm_slope_Bsor
p_ddm_slope_Bsor
# ggsave(filename = "taxonomic diversity/figs/p_ddm_slope_Bsor.tiff", 
#       device = "tiff", dpi = 600, height = 4, width = 3.5)

# Make the formal comparison.
Bsor_slopes$slope[1:10000] - Bsor_slopes$slope[10001:20000] -> slope_dif_Bsor
length(which(slope_dif_Bsor < 0))
length(which(slope_dif_Bsor > 0)) # Smaller one.
(length(which(slope_dif_Bsor > 0)) / 10000)*2
# p-value for this comparison is zero (p = 0) and and the ddm of 2018 
# has a significantly steeper slope.

# Comparisons for Bsne models won't be made as the best Bsne model for 2018 does
# not have a significant slope.

# Check the intercepts now.
data.frame(intercept = c(boot_ddm_1998_sim$boot.coefs[, 1], 
                     boot_ddm_2018_sim$boot.coefs[, 1]), 
           year = as.factor(c(rep("1998", 10000), 
                              rep("2018", 10000)))) -> Bsim_intercepts
data.frame(intercept = c(boot_ddm_1998_sor$boot.coefs[, 1], 
                     boot_ddm_2018_sor$boot.coefs[, 1]), 
           year = as.factor(c(rep("1998", 10000), 
                              rep("2018", 10000)))) -> Bsor_intercepts
Bsim_intercepts$intercept[1:10000] - Bsim_intercepts$intercept[10001:20000] -> 
  intercept_dif_Bsim
length(which(intercept_dif_Bsim < 0)) # Smaller one.
length(which(intercept_dif_Bsim > 0)) 
(length(which(intercept_dif_Bsim < 0)) / 10000)*2
# The difference in intercepts is significant (p = 0.0306) with 2018's 
# communities having a higher intercept for Bsim ddm.
Bsor_intercepts$intercept[1:10000] - Bsor_intercepts$intercept[10001:20000] -> 
  intercept_dif_Bsor
length(which(intercept_dif_Bsor < 0)) 
length(which(intercept_dif_Bsor > 0)) # Smaller one.
(length(which(intercept_dif_Bsor > 0)) / 10000)*2
# The difference in intercepts is significant (p = 0) with 2018's communities 
# having a higher intercept for Bsor ddm.


## 3.2: COMPONENTS OF BETA DIVERSITY ####


# Compare the beta diversity and its components between the study years by using
# GLMMs. Start with visualizing the data.
int1 <- 1 # Intermediary objects for the while loops (see below).
int2 <- 1
col_id <- c(); row_id <- c()
# The loops aim to create column ID and row ID columns to be used as random 
# effects in the GLMM.
while(int1 <= 33){
  col_id <- c(col_id, rep(int1, 34-int1))
  int1 <- int1 + 1
}
col_id
while(int2 <= 33){
  row_id <- c(row_id, int2:33)
  int2 <- int2 + 1
}
row_id
# Create a data frame that will have all the data.
data.frame(index = c(rep("Bsim", 1122), rep("Bsne", 1122), rep("Bsor", 1122)), 
           year = as.factor(rep(c(rep("1998", 561), rep("2018", 561)), 3)),
           identifier = as.factor(rep(rep(paste("X", 1:561, sep = ""), 2), 3)),
           value = c(as.vector(pairwise_1998$beta.sim),
                     as.vector(pairwise_2018$beta.sim),
                     as.vector(pairwise_1998$beta.sne),
                     as.vector(pairwise_2018$beta.sne),
                     as.vector(pairwise_1998$beta.sor),
                     as.vector(pairwise_2018$beta.sor)),
           col_ID = rep(col_id, 6),
           row_ID = rep(row_id, 6)) -> beta_data
# Duplicate the data frame and transform the beta diversity data so that there 
# are no zeroes left. This is because a log link will be used within gamma 
# family. 
beta_data_transformed <- beta_data
beta_data_transformed$value <- beta_data_transformed$value + 0.001
# Build the model.
m_Bsim <- glmmTMB(data = beta_data[1:1122, ], 
             formula = value ~ year + 
               (1 | row_ID) + (1 | col_ID), family = Gamma(link = "log"))
summary(m_Bsim) # The effect is significant. Validated the model.
# 1) Check the distribution of the residuals around zero.
boxplot(residuals(m_Bsim)); abline(h = 0, lty = 2, col = "red")
# Nearly perfect.
# 2) Plot the residuals against the explanatory variable.
boxplot(residuals(m_Bsim) ~ m_Bsim$frame$year)
abline(h = 0, lty = 2, col = "red")
# The model performs fractionally better for the 1998 level, but the difference
# does not seem to be substantial.
# 3) Check the residuals for autocorrelation.
par(mfrow = c(1, 2))
acf(residuals(m_Bsim)); pacf(residuals(m_Bsim)) 
par(mfrow = c(1, 1))
# There are some significant lags in the plots. The model is not valid.

# The model is not valid. The other possible formats, including a GLMM with cell
# ID as a random effect, all had independence issues and it seems like GLMMs 
# won't be able to solve the independence issues. Therefore, we will use 
# Wilcoxon signed-rank test, which is a non-parametric test, to compare the
# dissimilarities in the two matrices.

# Recalculate the beta diversity components so that we have the untransformed 
# values. 
beta.pair(study_comm_1998, index.family = "sorensen") -> pairwise_1998
beta.pair(study_comm_2018, index.family = "sorensen") -> pairwise_2018
wilcox.test(x = as.vector(pairwise_1998$beta.sim),  
            y = as.vector(pairwise_2018$beta.sim),
            alternative = "two.sided")
# 2018's communities have significantly higher Bsim values (p < 2.2e-16).
wilcox.test(x = as.vector(pairwise_1998$beta.sne),  
            y = as.vector(pairwise_2018$beta.sne),
            alternative = "two.sided")
# 1998's communities have significantly higher Bsne values (p = 5.075e-08).
wilcox.test(x = as.vector(pairwise_1998$beta.sor),  
            y = as.vector(pairwise_2018$beta.sor),
            alternative = "two.sided")
# 2018's communities have significantly higher Bsor values (p < 2.2e-16).

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
 ggsave(filename = "taxonomic diversity/figs/p_taxonomic_beta.tiff", dpi = 400, 
        units = "mm", width = 120, height = 75)

# ggsave(p_functional_beta, filename = "functional diversity/figs/func_beta.tiff",
#       units = "mm", dpi = 400, width = 120, height = 75)


## 3.3: PAIRWISE COMPARISONS FOR THE SAME SQUARES ####


# Aim of the section is to compare each square, by using taxonomic beta 
# diversity measures, to see where the biggest changes in species composition 
# took place, and to relate this to grid squares' functional diversity. 

# Firstly, calculate pairwise taxonomic beta diversity for the same squares. We
# are actually comparing each square to its past/future-self.
rownames(study_comm_1998) == rownames(study_comm_2018)
# Define a function that will calculate pairwise beta diversity for each square.
func_pairwise_beta <- function(x){
  study_comm_1998[x, ] -> a
  study_comm_2018[x, ] -> b
  rbind.data.frame(a, b, make.row.names = T) -> c
  beta.pair(c)
}
# Run the function.
as.data.frame(t(sapply(X = 1:34, FUN = func_pairwise_beta))) -> pw_beta_self
# Save the output to be used later in the functional diversity script.
write.csv2(x = as.matrix(pw_beta_self), 
           file = "taxonomic diversity/output/self_beta.csv")




