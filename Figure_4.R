library(ggplot2)
library(scales)
# show_col(hue_pal()(5))
library(gridExtra)
library(gridExtra)
library(tidyverse)

rm(list = ls())

#Scenarios considered:

#---Delayed effects

tau = 4
medianE1 = 6
medianE2 = 9
medianC = 6
time = seq(0,by=0.001,10)
hr_delayed_effects = numeric()
for(i in 1:length(time)){
  hr_delayed_effects[i] = ifelse(time[i] <= tau, medianC/medianE1, medianC/medianE2)
}

#---Crossing hazards

tau = 4
medianE1 = 9
medianE2 = 4
medianC = 6
time = seq(0,by=0.001,10)
hr_crossing = numeric()
for(i in 1:length(time)){
  hr_crossing[i] = ifelse(time[i] <= tau, medianC/medianE1, medianC/medianE2)
}

plot(time,hr_crossing,type="l")

#---Decreasing effects

tau = 4
medianC1 = 6
medianC2 = 9
medianE = 9
time = seq(0,by=0.001,10)
hr_decreasing = numeric()
for(i in 1:length(time)){
  hr_decreasing[i] = ifelse(time[i] <= tau, medianC1/medianE, medianC2/medianE)
}


df_delayed = data.frame(time = time, hr = hr_delayed_effects)

df_crossing = data.frame(time = time, hr = hr_crossing)

df_decreasing = data.frame(time = time, hr = hr_decreasing)

#Plots

p1 <- ggplot(aes(x = time, y = hr), data = df_delayed) +
  geom_hline(yintercept=1, color = "grey", linetype = "dashed", size = 1) +
  ylim(0.6,1.1) +
  geom_line(size = 1) +
  scale_x_continuous("Time (months)") +
  labs(y = "Hazard ratio", color = "",title = "A)") +
  scale_color_manual(values=c("black", "black")) +
  theme_light() +
  theme(text = element_text(size = 16), legend.position = "none")
p1

p2 <- ggplot(aes(x = time, y = hr), data = df_crossing) +
  geom_hline(yintercept=1, color = "grey", linetype = "dashed", size = 1) +
  geom_line(size = 1) +
  scale_x_continuous("Time (months)") +
  labs(y = "Hazard ratio", color = "",title = "B)") +
  scale_color_manual(values=c("black", "black")) +
  theme_light() +
  theme(text = element_text(size = 16), legend.position = "none") +
  ylim(0.5,1.6)
p2

p3 <- ggplot(aes(x = time, y = hr), data = df_decreasing) +
  geom_hline(yintercept=1, color = "grey", linetype = "dashed", size = 1) +
  geom_line(size = 1) +
  scale_x_continuous("Time (months)") +
  labs(y = "Hazard ratio", color = "",title = "C)") +
  scale_color_manual(values=c("black", "black")) +
  theme_light() +
  theme(text = element_text(size = 16), legend.position = "none") +
  ylim(0.6,1.1)
p3


scenarios = grid.arrange(p1, p2, p3,
                         nrow = 1, ncol = 3)
class(scenarios) <- c("scenarios", class(scenarios))

ggsave("/Users/jose/Dropbox/Research/2019/RMST as an alternative to HR in clinical trials with delayed effects/Journal of Applied Statistics - 1st Review/Figure_4.pdf",scenarios,
       width = 10, height = 4, dpi = 300)
