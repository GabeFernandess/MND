#--- LIBRARIES -----------------------------------------------------------
library(readxl)
library(tidyverse)
library(dplyr)
library(nlme)   
library(emmeans)
library(MuMIn)
library(ggbeeswarm)
library(viridis)
library(knitr)
library(lme4)
library(lmerTest)
library(sjstats)
library(car)
library(pwr)
library(robustlm)
library(paletteer)
library(tidyr)
library(patchwork)

#data frame
rms_data  <- read_excel("RMS.xlsx")

colnames(rms_data)

rms_data <- rms_data%>% 
mutate(
  Visit         = factor(Visit),
    Side        = factor(Side),
    MRC         = factor(MRC,   levels = c("5","4","3","2","1")),
    Participant = factor(Participant),
    ALSFRS      = factor(ALSFRS),
    sex         = factor(sex,   levels = c("male","female")),
    Time        = as.numeric(Time)
  ) %>% 
  filter(!is.na(MRC)) %>%
  mutate(
  
      limb = case_when(
        MRC == "5"        ~ "Pre-symptomatic",
        !is.na(MRC)       ~ "Symptomatic",     
        TRUE              ~ NA_character_     
      ),
      limb = factor(limb, levels = c("Pre-symptomatic", "Symptomatic")),
      time = Time/4.345 #time in months 
    )

#"best" contraction per Participant?Visit?Side ----------------
contract_means <- rms_data %>%
  group_by(Participant, Visit, Side) %>%
  summarise(
    mean_RMS = mean(meanRMS250ms_DD , na.rm=TRUE),
    limb       = first(limb),
    .groups    = "drop"
  )



rms_collapsed <- rms_data %>%
  group_by(Participant, Side, limb, time) %>%
  summarise(
    mean_RMS = mean(meanRMS250ms_DD, na.rm=TRUE),
            .groups = "drop"
    )

plot_data <- rms_collapsed 

rms_data %>%
  group_by(Participant) %>%
  summarise(n_time = n_distinct(time), .groups="drop") %>%
  count(n_time)


#-----------------

#--- MODEL-----------------------------------------------------

model <- lmer(meanRMS250ms_DD ~ time * limb + Side + (1 | Participant), data=rms_data)



anova(model)

#--- 5. SUMMARIES & EFFECT SIZES -----------------------------------------
summary(model)
anova(model) %>% kable(digits=3)

overall_trend <- emtrends(model, ~ 1, var = "time")
summary(overall_trend, infer = TRUE)  


summary(overall_trend, infer = TRUE) %>%
  as.data.frame() %>%
  dplyr::select(time.trend, SE, df, lower.CL, upper.CL, t.ratio, p.value) %>%
  knitr::kable(digits = 3)

t_seq <- seq(min(rms_data$time, na.rm = TRUE),
             max(rms_data$time, na.rm = TRUE),
             length.out = 100)
emm_df <- emmeans(
  model,specs = ~ 1 | ~ time,
  at = list(time = t_seq),
  weights = "proportional"
) %>% as.data.frame()


id_fac    <- interaction(plot_data$Participant, plot_data$Side, drop=TRUE)
n_ids     <- nlevels(id_fac)
base_pal  <- viridis(256, option="turbo")
subset_pal <- base_pal[1:round(256*0.85)]
palette_all <- colorRampPalette(subset_pal)(n_ids)

plot_data <- plot_data %>%
  mutate(colour = palette_all[as.integer(id_fac)])
rms_collapsed <- rms_collapsed %>%
  mutate(colour = palette_all[as.integer(interaction(Participant,Side,drop=TRUE))])


#plot _RMS

rms_plot_points <- ggplot() +
   geom_line(data = plot_data,
            aes(x=time, y=mean_RMS,
                group=interaction(Participant,Side),
                colour = colour),
            size=0.9, alpha=0.3) +
  geom_point(data = plot_data,
             aes(x=time, y=mean_RMS, fill = colour, colour = colour),
             shape=21, size=3, stroke=0.8, alpha=0.65) +
  
  scale_colour_identity() +
  scale_fill_identity() +
  
  geom_ribbon(data=emm_df,
              aes(x=time, ymin=lower.CL, ymax=upper.CL),
              fill="black", alpha=0.15) +
  geom_line(data=emm_df,
            aes(x=time, y=emmean),
            colour="black", size=1.3) +
  
  theme_light(base_size=22) +
  theme(
    panel.grid    = element_blank(),
    axis.text.x   = element_text( hjust=1),
    legend.position="none") +
  labs(
    x = "Time (months from baseline)",
    y = "RMS-EMG (mV)"
  )

p1 <- rms_plot_points

p1

# --- TIME TREND -------------------------------------

time_trend <- emtrends(model, ~ 1, var = "time", weights = "proportional")
tt_df <- summary(time_trend, infer = TRUE) %>%
  as.data.frame() %>%
  dplyr::rename(
    slope_month = time.trend,
    lwr_month   = lower.CL,
    upr_month   = upper.CL
  )

emm0 <- emmeans(model, ~ 1, at = list(time = 0), weights = "proportional") %>%
  as.data.frame()
b0 <- emm0$emmean[1]

t_seq <- seq(min(rms_data$time, na.rm = TRUE),
             max(rms_data$time, na.rm = TRUE),
             length.out = 200)

trend_line_df <- tibble::tibble(
  time    = t_seq,
  fit     = b0 + tt_df$slope_month[1] * t_seq,
  lower   = b0 + tt_df$lwr_month[1]   * t_seq,
  upper   = b0 + tt_df$upr_month[1]   * t_seq
)

ann_df <- tibble::tibble(
  label = sprintf("Slope = %.2f (%.2f, %.2f) mV/month",
                  tt_df$slope_month[1], tt_df$lwr_month[1], tt_df$upr_month[1]),
  x = Inf, y = Inf
)


#plot _RMS_slope

rms_plot_model_slopes <- rms_plot_points +
  ggplot2::geom_ribbon(
    data = trend_line_df,
    ggplot2::aes(x = time, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    fill = "black", alpha = 0.1
  ) +
  ggplot2::geom_line(
    data = trend_line_df,
    ggplot2::aes(x = time, y = fit),
    inherit.aes = FALSE,
    colour = "black", linewidth = 1.5
  ) +
  ggplot2::geom_text(
    data = ann_df,
    ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.12, vjust = 3, size = 5
  )

p1.1 <- print(rms_plot_model_slopes)


#ggsave(file = "rms_mvic_months.tiff", units="in", width = 14, height= 9, dpi = 600, compression = "lzw")

