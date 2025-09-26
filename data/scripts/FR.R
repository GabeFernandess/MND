
# Libraries
library(readxl)
library(tidyverse)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(sjstats)
library(car)
library(pwr)
library(viridis)
library(paletteer)
library(tidyr)
library(patchwork)
library(ggbeeswarm)
library(MuMIn)
library(nlme)
library(knitr)   

#mean FR - using time in weeks from baseline as variable 

#------------------

# Set the path to your Excel file
mu_data <- read_excel("C:/Users/gfernand/OneDrive/Research projects/MVIC/paper/github/MU.xlsx")
  mutate(
    Visit        = factor(Visit, levels = c("v1","v2","v3","v4","v5")),
    Side         = factor(Side,  levels = c("L","R")),
    MRC          = factor(MRC,   levels = c("5","4","3","2","1")),
    Participant  = factor(Participant),
    MU           = factor(MU),
    sex          = factor(sex,        levels = c("male","female")),
    mu_id        = interaction(Participant, Contraction, Side, MU, drop = TRUE),
    Time         = as.numeric(Time)
  )


mu_data <- mu_data %>%
  filter(!is.na(MRC))

colnames(mu_data)



#DATA MANIPULATION ----

names(mu_data)[names(mu_data) ==  "Time"] <- "time"
names(mu_data)[names(mu_data) ==  "MeanIFR_250ms"] <- "Mean_FR_250ms"

mu_data$limb <- ifelse(mu_data$MRC ==5, "Pre-symptomatic", "Symptomatic")


mu_data$limb <- factor(mu_data$limb,
                       levels = c("Pre-symptomatic", "Symptomatic"))



levels(as.factor(mu_data$Visit))
levels(as.factor(mu_data$limb))

#####-----
#####-----Demographics


library(tidyr)

# 0) Filter to baseline
mu0 <- mu_data %>% filter(time == 0)

# 1) Total MU rows at baseline
total_mu_rows_baseline <- nrow(mu0)

total_mu_rows_baseline

# 2) Distinct participants at baseline + who has L / R / both
#    First, which sides each participant has at baseline:
pt_side <- mu0 %>%
  distinct(Participant, Side) %>%
  group_by(Participant) %>%
  summarise(
    has_L   = any(Side == "L"),
    has_R   = any(Side == "R"),
    n_sides = n_distinct(Side),
    .groups = "drop"
  )

participants_summary_baseline <- tibble(
  n_participants = n_distinct(mu0$Participant),
  n_with_L_only  = sum(pt_side$has_L & !pt_side$has_R),
  n_with_R_only  = sum(pt_side$has_R & !pt_side$has_L),
  n_with_both    = sum(pt_side$n_sides == 2),
  n_distinct_limbs = nrow(distinct(mu0, Participant, Side))
)

#####-----
#####-----

#--- 1. CODE : "best" contraction per Participant?Visit?Side ----------------
contract_means <- mu_data %>%
  group_by(Participant, Visit, Side, Contraction) %>%
  summarise(
    mean_FR250 = mean(Mean_FR_250ms, na.rm=TRUE),
    ALSFRS     = first(ALSFRS),
    limb       = first(limb),
    .groups    = "drop"
  )

highest_FR <- contract_means %>%
  group_by(Participant, Visit, Side) %>%
  slice_max(mean_FR250, n=1, with_ties=FALSE) %>%
  ungroup()

# now filter mu_data to only those rows
mu_data <- mu_data %>%
  inner_join(highest_FR,
             by = c("Participant","Visit","Side","Contraction","ALSFRS","limb"))

#time diff between visits
mu_data %>%
  select(Participant, Visit, time) %>%
  distinct() %>%
  arrange(Participant, time) %>%
  group_by(Participant) %>%
  mutate(
    time_diff_weeks = time - lag(time)
  ) %>%
  filter(!is.na(time_diff_weeks)) %>%
  ungroup() %>%
  summarise(
    mean_diff_weeks = mean(time_diff_weeks),
    SD   = sd(time_diff_weeks),
  )


mu_data <- mu_data %>%
  mutate(
    time = time /4.345
  )

# 2. Collapse to MU level then Participant level
mu_collapsed <- mu_data %>%
  group_by(Participant, Side, Visit, time, limb, MU) %>%
  summarise(mean_FR_MU = mean(Mean_FR_250ms, na.rm = TRUE), .groups = "drop")


plot_data <- mu_collapsed %>%
  group_by(Participant, Side, Visit, limb, time) %>%
  summarise(mean_FR = mean(mean_FR_MU, na.rm = TRUE), .groups = "drop")




# 2) Also allow different residual variances in Strong vs. Weak limbs
model_lme_full <- lme(
  fixed       = Mean_FR_250ms ~ time * limb + Side,
  random      = ~ 1 | Participant,
  data        = mu_data,
  correlation = corCompSymm(form = ~time | Participant),
  weights     = varIdent(form = ~1 | limb),
  na.action   = na.exclude
)
summary(model_lme_full)
anova(model_lme_full) %>% 
  kable(digits = 3)

model <- model_lme_full

#-----
#-------------

# 4. Summaries

summary(model)
anova(model) %>% 
  kable(digits = 3)
r.squaredGLMM(model)

# 5. Time trend estimation (main effect)
time_trend <- emtrends(model, specs = ~ 1, var = "time")
summary(time_trend)
test(time_trend)
# Extract slope for time
# Extract slope correctly
tt <- summary(emtrends(model, ~ 1, var = "time"))
slope <- tt$time.trend[1]   # This is the slope estimate

# Residual SD
resid_sd <- sigma(model)

# Cohen's d
cohens_d <- slope / resid_sd
cohens_d

tt   # check the slope value
resid_sd  # check the residual SD

sigma(model)

# Get raw residuals and compute SD
resid_sd_raw <- sd(residuals(model), na.rm = TRUE)

# Now compute d
cohens_d <- slope / resid_sd_raw
cohens_d


# values from your model
slope     <- -0.3160   # 
SE_slope  <- 0.0623    # 
resid_sd  <- 10.59107 #<---sigma(model)

# Cohen's d
cohens_d  <- slope / resid_sd

# CI for the slope first
alpha <- 0.05
z_crit <- qnorm(1 - alpha/2)
lower_slope <- slope - z_crit * SE_slope
upper_slope <- slope + z_crit * SE_slope

# CI for Cohen's d
lower_d <- lower_slope / resid_sd
upper_d <- upper_slope / resid_sd

tibble(
  cohens_d = cohens_d,
  lower_95 = lower_d,
  upper_95 = upper_d
)


# 6. Colour palette for participants
# 1. Get all unique participants across v1–v5
id_fac <- interaction(plot_data$Participant, plot_data$Side, drop = TRUE)
n_ids  <- nlevels(id_fac)

# 2. Generate a palette for all participants
base_palette   <- viridis(256, option = "viridis")
subset_palette <- base_palette[1:round(256 * 0.85)]
palette_all    <- colorRampPalette(subset_palette)(n_ids)

plot_data <- plot_data %>%
  mutate(colour = palette_all[as.integer(id_fac)])

mu_collapsed <- mu_collapsed %>%
  mutate(id_fac = interaction(Participant, Side, drop = TRUE),
         colour = palette_all[as.integer(id_fac)])

# 7. Predict fixed effects across time (averaged across Side & Contraction)
newd <- expand.grid(
  time = seq(min(mu_data$time, na.rm = TRUE),
             max(mu_data$time, na.rm = TRUE),
             length = 100),
  Side = levels(mu_data$Side),
  limb = levels (mu_data$limb),
  Contraction = unique(mu_data$Contraction)
)

# for lme: level = 0 means “no random effects” (equivalent to re.form = NA in lmer)
newd$pred <- predict(model, newdata = newd, level = 0)

newd_avg <- newd %>%
  group_by(time) %>%
  summarise(pred = mean(pred), .groups = "drop")

# 8. Extract modelled emmeans + CI
time_vals <- seq(min(mu_data$time,   na.rm=TRUE),
                 max(mu_data$time,   na.rm=TRUE),
                 length.out = 100)

newd <- expand.grid(
  time        = time_vals,
  limb        = levels(mu_data$limb)
)


emm_df <- as.data.frame(
  emmeans(model, ~ time, at = list(time = time_vals))
)



## 1. Points + participant trajectories + marginal means (no ribbon/line)
mvic_plot_points <- ggplot(plot_data, aes(x = time, y = mean_FR, colour = colour)) +
  # raw MU means
  #geom_point(data = mu_collapsed,
  #     aes(x = time, y = mean_FR_MU, colour = colour),
  #   shape = 18, size = 5.5, alpha = 0.05,
  #  position = position_nudge(x = -0.7)) +
  scale_colour_identity() +
  # individual trajectories
  geom_line(aes(group = interaction(Participant, Side)),
            linewidth = 0.9, alpha = 0.2) +
  # marginal means
  geom_point(aes(fill = colour, colour = colour),
             shape = 23, size = 3, stroke = 1, alpha = 0.6) +
  scale_fill_identity() +
  theme_light(base_size = 22) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(size = 20),
    axis.title.x     = element_text(size = 20),
    axis.title.y     = element_text(size = 20, vjust = 0.5, margin = margin(r = 20)),
    strip.text.x     = element_text(size = 20),
    plot.margin      = margin(t = 10, r = 10, b = 10, l = 20),
    legend.position  = "none"
  ) +
  labs(
    x = "Time (months from baseline)",
    y = "Maximal firing rates (Hz)"
  ) +
  scale_x_continuous(
    breaks = pretty(plot_data$time, n = 5)
  )

print(mvic_plot_points)



# 2. Model‐predicted emmeans + CI ribbon only

# 1. Compute R² (marginal & conditional)
r2 <- r.squaredGLMM(model)
r2m <- round(r2[1], 3)
r2c <- round(r2[2], 3)

r2m *100
r2c *100


# 1) Slopes (weeks) + 95% CI, then convert to /month
time_trend <- emtrends(model,specs =~1, var = "time")
tt_df <- as.data.frame(summary(time_trend, infer = TRUE))

# rename & convert
tt_df <- tt_df %>%
  rename(slope_month = time.trend,
         lwr_month   = lower.CL,
         upr_month  = upper.CL) 

names(tt_df)

# 2) Annotation labels (slope ± CI)

ann_df <- emtrends(model, specs= ~1, var= "time") %>%
  as.data.frame() %>%
  rename(
    slope      = time.trend,
    lower.CL   = lower.CL,
    upper.CL   = upper.CL
  ) %>%
  mutate(   slope_month   = slope, 
            lwr_month     = lower.CL, 
            upr_month     = upper.CL,
            label = sprintf(
              "%.2f (%.2f,%.2f) Hz/month",
              slope_month,
              lwr_month,
              upr_month
            ),
            # pick x/y positions for your panel (adjust as needed)
            x = max(emm_trim$time) - 1,
            y = max(emm_trim$emmean) + 1)

summary(time_trend)

# 1. Ranges per limb
rng <- mu_data %>%
  summarise(tmin = min(time, na.rm = TRUE),
            tmax = max(time, na.rm = TRUE),
            .groups = "drop")

emm_trim <- map_dfr(seq_len(nrow(rng)), function(i){
  tt  <- seq(rng$tmin[i], rng$tmax[i], length.out = 100)
  
  tmp <- emmeans(
    model,
    ~ time,  # FIXED: include interaction
    at = list(time = tt,
              Side = levels(mu_data$Side),  # OK if these are model covariates
              Contraction = unique(mu_data$Contraction)),  # include only if used
    cov.reduce = FALSE
  ) %>%
    as.data.frame()
})

#------------------


mvic_plot_points2 <- mvic_plot_points +
  geom_ribbon(
    data        = emm_trim,
    inherit.aes = FALSE,
    aes(x = time, ymin = lower.CL, ymax = upper.CL),
    fill        = "#F4A000", alpha = 0.2
  ) +
  geom_line(
    data        = emm_trim,
    inherit.aes = FALSE,
    aes(x = time, y = emmean),
    colour      = "#F4A000", linewidth = 1.2
    
  )

mvic_plot_points2



  
mvic_plot_model_slopes3 <- mvic_plot_points2 +
  annotate(
    "text",
    x      = Inf,            
    y      = Inf,            
    label  = ann_df$label,   
    hjust  = 1.1,            
    vjust  = 4,          
    size   = 5
  ) 

print(mvic_plot_model_slopes3) 

#-----------------

#ggsave(file = "MUFR_mvic.tiff", units="in", width = 14, height= 9, dpi = 600, compression = "lzw")

