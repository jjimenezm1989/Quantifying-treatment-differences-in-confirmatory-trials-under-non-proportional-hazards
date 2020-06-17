library(ggplot2)
library(scales)
# show_col(hue_pal()(5))
library(gridExtra)
require("survival")
library(survminer)

rm(list = ls())


delayed_effect = function(delay_c = 4, delay_e = 4){
  
  n_c = 100000
  n_e = 100000 
  rec_period = 12 
  rec_power = 1
  rate_c_1 = log(2) / 6
  rate_c_2 = log(2) / 9
  rate_e_1 = log(2) / 6 
  rate_e_2 = log(2) / 9
  max_cal_t  = 36
  n_events = 258
  
  # simulate recruitment times from power model:
  
  rec_c = rec_period * runif(n_c) ^ (1 / rec_power)
  rec_e = rec_period * runif(n_e) ^ (1 / rec_power)
  
  # control event times are exponentially distributed:
  
  t_1_c = rexp(n_c, rate = rate_c_1)
  t_2_c = rexp(n_c, rate = rate_c_2)
  t_c = ifelse(t_1_c < delay_c, t_1_c, delay_c + t_2_c)
  
  # experimental event times come from 2-piece exponential distribution:
  
  t_1_e = rexp(n_e, rate = rate_e_1)
  t_2_e = rexp(n_e, rate = rate_e_2)
  t_e = ifelse(t_1_e < delay_e, t_1_e, delay_e + t_2_e)
  
  # (calendar) event times, relative to start of trial.
  
  cal_t_c = rec_c + t_c
  cal_t_e = rec_e + t_e
  
  if (is.null(max_cal_t)){
    max_cal_t <- sort(c(cal_t_c, cal_t_e))[n_events]
  }
  
  # does the patient have an event before the data cut-off:
  
  event_c = cal_t_c <= max_cal_t
  event_e = cal_t_e <= max_cal_t
  
  # if patient's event time is censored, calculate their follow-up time:
  
  obs_t_c = ifelse(event_c, t_c, max_cal_t - rec_c)
  obs_t_e = ifelse(event_e, t_e, max_cal_t - rec_e)
  
  # store in data frame with group label:
  
  df = data.frame(time = c(obs_t_c, obs_t_e),
                  event = c(event_c, event_e),
                  group = rep(c("control", "experimental"), c(n_c, n_e)))
  
  # round time to 2 dp
  df$time = round(df$time, 2)
  
  df
}


medianE1 = 6
medianE2 = 9
medianC1 = 6
medianC2 = 9
time = seq(0,by=0.01,10)

#null deltaC = deltaE
tauC = 4
tauE = 4
hC =  ifelse(time < tauC, log(2)/medianC1, log(2)/medianC2)
hE =  ifelse(time < tauE, log(2)/medianE1, log(2)/medianE2)
HR = hE/hC

df_all_equal = data.frame(hazard = c(hE, hC),
                          type = c(rep("Experimental",length(time)),rep("Control",length(time))),
                          time = rep(time,2))

df_all_equal$type = factor(df_all_equal$type,levels = c("Experimental","Control"), labels = c("Experimental","Control"),ordered = TRUE)

p1 <- ggplot(aes(x = time, y = hazard, color = type), data = df_all_equal) +
  geom_line(size = 1.5) +
  labs(y = "Hazard", x = "Time (months)", color = "", title = "A)") +
  scale_color_manual(values=c("#F8766D", "#619CFF")) +
  theme_light() +
  theme(text = element_text(size = 16),
        legend.position = "none")

sim_data_all_equal = delayed_effect(delay_c = 4, delay_e = 4)

fit_all_equal <- survfit(Surv(time, event) ~ group, data = sim_data_all_equal)

s1 = ggsurvplot(fit_all_equal, data = sim_data_all_equal,
           censor = FALSE,
           legend.labs = c("Control", "Experimental"),
           legend.title = "",
           palette = c("#619CFF","#F8766D"),
           legend = "none",
           title = "A)")

survival_plot1 = s1$plot + geom_vline(xintercept = 4, linetype='dashed', color = "#619CFF") + geom_vline(xintercept = 4, linetype='dashed', color = "#F8766D")


#null deltaE >= deltaC
tauC = 4
tauE = 6
hC =  ifelse(time < tauC, log(2)/medianC1, log(2)/medianC2)
hE =  ifelse(time < tauE, log(2)/medianE1, log(2)/medianE2)
HR = hE/hC

df_EgeqC = data.frame(hazard = c(hE, hC),
                          type = c(rep("Experimental",length(time)),rep("Control",length(time))),
                          time = rep(time,2))

df_EgeqC$type = factor(df_EgeqC$type,levels = c("Experimental","Control"), labels = c("Experimental","Control"),ordered = TRUE)


p2 <- ggplot(aes(x = time, y = hazard, color = type), data = df_EgeqC) +
  geom_line(size = 1.5) +
  labs(y = "Hazard", x = "Time (months)", color = "", title = "B)") +
  scale_color_manual(values=c("#F8766D", "#619CFF")) +
  theme_light() +
  theme(text = element_text(size = 16),
        legend.position = c(0.8, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent"))


sim_data_EgeqC = delayed_effect(delay_c = 4, delay_e = 6)

fit_EgeqC <- survfit(Surv(time, event) ~ group, data = sim_data_EgeqC)

s2 = ggsurvplot(fit_EgeqC, data = sim_data_EgeqC,
           censor = FALSE,
           legend.labs = c("Control", "Experimental"),
           legend.title = "",
           palette = c("#619CFF","#F8766D"),
           legend = c(0.8,0.9),
           font.legend = c(13,"italic"),
           title = "B)")

survival_plot2 = s2$plot + geom_vline(xintercept = 4, linetype='dashed', color = "#619CFF") + geom_vline(xintercept = 6, linetype='dashed', color = "#F8766D")


#plots

#Hazard plots
scenarios = grid.arrange(p1, p2,
                         nrow = 1, ncol = 2)
class(scenarios) <- c("scenarios", class(scenarios))

ggsave("/Users/jose/Dropbox/Research/2019/RMST as an alternative to HR in clinical trials with delayed effects/Journal of Applied Statistics - 1st Review/Figure_9.pdf",scenarios,
       width = 10, height = 4, dpi = 300)

#Survival plots
scenarios = grid.arrange(survival_plot1, survival_plot2,
                         nrow = 1, ncol = 2)
ggsave("/Users/jose/Dropbox/Research/2019/RMST as an alternative to HR in clinical trials with delayed effects/Journal of Applied Statistics - 1st Review/Figure_10.pdf",scenarios,
       width = 10, height = 4, dpi = 300)