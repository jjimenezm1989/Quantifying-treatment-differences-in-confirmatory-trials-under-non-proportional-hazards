library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(gridExtra)
library(survival)
library(survRM2)
library(survminer)
# show_col(hue_pal()(4))

rm(list = ls())

delayed_effect_sim = function(n_c = 165, n_e = 165, rec_period = 12, rec_power = 1, med_c = 6, rate_e_1 = log(2) / 6, rate_e_2 = log(2) / 9, delay = 2, max_cal_t  = 36, n_events = NULL, model = NULL){
  
  if (is.null(max_cal_t) && is.null(n_events)) stop("either max_cal_t or n_events must be specified.")
  if ((!is.null(max_cal_t)) && (!is.null(n_events))) stop("one of max_cal_t and n_events must be NULL.")
  if (is.null(max_cal_t) && (n_events > n_c + n_e)) stop("number of events not reached.")
  
  if( !is.null(model))
  {
    n_c <- model$n_c
    n_e <- model$n_e
    rec_period = model$rec_period
    rec_power = model$rec_power
    med_c = model$med_c
    rate_e_1 = model$rate_e_1
    rate_e_2 = model$rate_e_2 
    delay = model$delay 
    max_cal_t  = model$max_cal_t 
    n_events = model$n_events
  }
  
  # simulate recruitment times from power model:
  
  rec_c = rec_period * runif(n_c) ^ (1 / rec_power)
  rec_e = rec_period * runif(n_e) ^ (1 / rec_power)
  
  # control event times are exponentially distributed:
  
  t_c = rexp(n_c, rate = log(2) / med_c)
  
  # experimental event times come from 2-piece exponential distribution:
  
  t_1_e = rexp(n_e, rate = rate_e_1)
  t_2_e = rexp(n_e, rate = rate_e_2)
  t_e = ifelse(t_1_e < delay, t_1_e, delay + t_2_e)
  
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

crossing_hazards_sim = function(n_c = 165, n_e = 165, 
                                rec_period = 12, 
                                rec_power = 1, 
                                med_c = 6, 
                                rate_e_1 = log(2) / 6, 
                                rate_e_2 = log(2) / 9, 
                                threshold = 2, 
                                max_cal_t  = 36, 
                                n_events = NULL, 
                                model = NULL){
  
  if (is.null(max_cal_t) && is.null(n_events)) stop("either max_cal_t or n_events must be specified.")
  if ((!is.null(max_cal_t)) && (!is.null(n_events))) stop("one of max_cal_t and n_events must be NULL.")
  if (is.null(max_cal_t) && (n_events > n_c + n_e)) stop("number of events not reached.")
  
  if( !is.null(model))
  {
    n_c <- model$n_c
    n_e <- model$n_e
    rec_period = model$rec_period
    rec_power = model$rec_power
    med_c = model$med_c
    rate_e_1 = model$rate_e_1
    rate_e_2 = model$rate_e_2 
    threshold = model$threshold 
    max_cal_t  = model$max_cal_t 
    n_events = model$n_events
  }
  
  # simulate recruitment times from power model:
  
  rec_c = rec_period * runif(n_c) ^ (1 / rec_power)
  rec_e = rec_period * runif(n_e) ^ (1 / rec_power)
  
  # control event times are exponentially distributed:
  
  t_c = rexp(n_c, rate = log(2) / med_c)
  
  # experimental event times come from 2-piece exponential distribution:
  
  t_1_e = rexp(n_e, rate = rate_e_1)
  t_2_e = rexp(n_e, rate = rate_e_2)
  t_e = ifelse(t_1_e < threshold, t_1_e, threshold + t_2_e)
  
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

decreasing_hazards_sim = function(n_c = 165, n_e = 165, 
                                  rec_period = 12, 
                                  rec_power = 1, 
                                  med_e = 9, 
                                  rate_c_1 = log(2) / 6, 
                                  rate_c_2 = log(2) / 9, 
                                  threshold = 2, 
                                  max_cal_t  = 36, 
                                  n_events = NULL, 
                                  model = NULL){
  
  if (is.null(max_cal_t) && is.null(n_events)) stop("either max_cal_t or n_events must be specified.")
  if ((!is.null(max_cal_t)) && (!is.null(n_events))) stop("one of max_cal_t and n_events must be NULL.")
  if (is.null(max_cal_t) && (n_events > n_c + n_e)) stop("number of events not reached.")
  
  if( !is.null(model))
  {
    n_c <- model$n_c
    n_e <- model$n_e
    rec_period = model$rec_period
    rec_power = model$rec_power
    med_e = model$med_e
    rate_c_1 = model$rate_c_1
    rate_c_2 = model$rate_c_2 
    threshold = model$threshold 
    max_cal_t  = model$max_cal_t 
    n_events = model$n_events
  }
  
  # simulate recruitment times from power model:
  
  rec_c = rec_period * runif(n_c) ^ (1 / rec_power)
  rec_e = rec_period * runif(n_e) ^ (1 / rec_power)
  
  # control event times are exponentially distributed:
  
  t_e = rexp(n_e, rate = log(2) / med_e)
  
  # experimental event times come from 2-piece exponential distribution:
  
  t_1_c = rexp(n_c, rate = rate_c_1)
  t_2_c = rexp(n_c, rate = rate_c_2)
  t_c = ifelse(t_1_c < threshold, t_1_c, threshold + t_2_c)
  
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

rmst_zs = function(data,t_cut){
  
  df = as.data.frame(data)
  
  # #Calculate t_cut
  # data_minimax2 = df %>% group_by(group) %>% summarize(maxTime= max(time)) 
  # data_minimax3 = data_minimax2 %>% summarize(minimax = min(maxTime))
  # t_cut = data_minimax3$minimax
  
  if (length(t_cut) != 1) stop("t_cut must be length 1")
  
  ## use survival::survfit to do KM estimation by group:
  fit <- survival::survfit(Surv(time, event) ~ group, data = df)
  
  if(any(min(fit$time[fit$n.risk == 1]) < t_cut)) {
    
    return(NA)
    
    
  }else {
    ## use survRM2:: to do RMST estimation by group:
    results <- rmst2(time = df$time, 
                     status = as.numeric(df$event), 
                     arm = ifelse(df$group == "control", 0, 1),
                     tau = t_cut)
  }
  
  mean_1 <- results$RMST.arm1$rmst["Est."]
  se_1 <- results$RMST.arm1$rmst["se"]
  
  mean_0 <- results$RMST.arm0$rmst["Est."]
  se_0 <- results$RMST.arm0$rmst["se"]
  
  mean_diff <- mean_1 - mean_0
  se_diff <- sqrt(se_0 ^ 2 + se_1 ^2)
  
  z_diff <- mean_diff / se_diff
  
  z_diff
  
  
}
get_risk_table <- function(dt){
  
  # arrange the data set in increasing order of survival time:
  
  dt <- dt[order(dt$time),]
  
  # split into 2 data sets: one for control; one for experimental:
  
  dt_c <- dt[dt$group == "control",]
  dt_e <- dt[dt$group == "experimental",]
  
  # number of patients on each arm
  
  n_c <- length(dt_c$time)
  n_e <- length(dt_e$time)
  
  # the number of patients at risk will decrease by 1 after each event/censored observation.
  
  at_risk_c <- n_c - 1:n_c + 1
  at_risk_e <- n_e - 1:n_e + 1
  
  # create a risk table just for the control arm data...
  
  risk_table_c <- data.frame(t = dt_c$time,
                             n_c = at_risk_c,
                             d_c = as.numeric(dt_c$event))
  
  # ...where there no patients/events on the experimental arm:
  
  risk_table_c$d_e <- 0
  risk_table_c$n_e <- NA
  
  # create a risk table just for the experimental arm data...
  
  risk_table_e <- data.frame(t = dt_e$time,
                             n_e = at_risk_e,
                             d_e = as.numeric(dt_e$event))
  
  # ...where there are no patients/events on the control arm:
  
  risk_table_e$d_c <- 0
  risk_table_e$n_c <- NA
  
  # put the risk tables on top of each other...
  
  risk_table <- rbind(risk_table_c, risk_table_e)
  
  # ...and reorder by event/censoring times (across both arms):
  
  risk_table <- risk_table[order(risk_table$t),]
  
  # whenever is.na(n_e) == TRUE, this means that the event/censored observation on this
  # row was from the control arm. To fill in the number at risk on the experimental
  # arm we look at the subsequent row, repeating if necessary, until we find a row
  # where is.na(n_e) == FALSE.
  # similarly for n_c when is.na(n_c) == TRUE.
  
  risk_table <- risk_table %>% tidyr::fill(n_e, n_c, .direction = "up")
  
  # at the bottom of the risk table, it's still possible that is.na(n_e) == TRUE if
  # all subsequent events/censorings are from the control arm. In this case the
  # number at risk is zero. Similarly for the control arm.
  
  risk_table$n_c[is.na(risk_table$n_c)] <- 0
  risk_table$n_e[is.na(risk_table$n_e)] <- 0
  
  # now we deal with ties. We group together the data that have the same value of "t",
  # work out how many patients were at risk just prior to "t", and how many events
  # happened at "t":
  
  risk_table <- risk_table %>%
    group_by(t) %>%
    summarize(n_e = max(n_e),
              d_e = sum(d_e),
              n_c = max(n_c),
              d_c = sum(d_c)) %>%
    as.data.table()
  
  # we only keep the "t" where there was at least one event:
  
  risk_table <- risk_table[risk_table$d_e > 0 | risk_table$d_c > 0,]
  
  # calculate number of events, number at risk across arms.
  
  risk_table$n <- risk_table$n_e + risk_table$n_c
  risk_table$d <- risk_table$d_e + risk_table$d_c
  
  # calculate the number censored between consecutive event times:
  
  risk_table$l <- risk_table$n - risk_table$d - c(risk_table$n[-1], 0)
  risk_table$l_c <- risk_table$n_c - risk_table$d_c - c(risk_table$n_c[-1], 0)
  risk_table$l_e <- risk_table$n_e - risk_table$d_e - c(risk_table$n_e[-1], 0)
  
  # return the completed risk table:
  
  risk_table
  
}

# significance level (one-sided)
alpha <- 0.025

# number of simulations
M <- 10000

t_cut = seq(1,by=1,20)

# model parameters Delayed effects
model_DE <- list(n_c = 165,
                 n_e = 165,
                 rec_period = 12, 
                 rec_power = 1, 
                 med_c = 6, 
                 rate_e_1 = log(2) / 6, 
                 rate_e_2 = log(2) / 9, 
                 delay = 4, 
                 max_cal_t  = 36, 
                 n_events = 258)

# model parameters Crossing hazards
model_CH <- list(n_c = 165,
                 n_e = 165,
                 rec_period = 12, 
                 rec_power = 1, 
                 med_c = 6, 
                 rate_e_1 = log(2) / 9, 
                 rate_e_2 = log(2) / 4, 
                 threshold = 9, 
                 max_cal_t  = 36, 
                 n_events = 258)

# model parameters Decreasing hazards
model_DH <- list(n_c = 165, 
                 n_e = 165, 
                 rec_period = 12, 
                 rec_power = 1, 
                 med_e = 9, 
                 rate_c_1 = log(2) / 6, 
                 rate_c_2 = log(2) / 9, 
                 threshold = 9, 
                 max_cal_t  = 36, 
                 n_events = 258)

# model parameters PH
model_PH <- list(n_c = 165,
                 n_e = 165,
                 rec_period = 12, 
                 rec_power = 1, 
                 med_c = 6, 
                 rate_e_1 = log(2) / 9, 
                 rate_e_2 = log(2) / 9, 
                 delay = 0, 
                 max_cal_t  = 36, 
                 n_events = 258)

# simulate data

sim_data_DH <- replicate(M, decreasing_hazards_sim(model = model_DH), simplify = FALSE)
sim_data_CH <- replicate(M, crossing_hazards_sim(model = model_CH), simplify = FALSE)
sim_data_DE <- replicate(M, delayed_effect_sim(model = model_DE), simplify = FALSE) 
sim_data_PH <- replicate(M, delayed_effect_sim(model = model_PH), simplify = FALSE) 

result <- data.table(t_cut = t_cut,
                     power_rmstDH = rep(0, length(t_cut)),
                     power_rmstCH = rep(0, length(t_cut)),
                     power_rmstDE = rep(0, length(t_cut)),
                     power_rmstPH = rep(0, length(t_cut)))

for (i in 1:length(t_cut)) {
  
  
  print(paste0("Simulation progress: ", i, " out of ", length(t_cut)))
  
  # calculate standard logrank
  rmst_Z_DH <- lapply(sim_data_DH, rmst_zs, t_cut=t_cut[i])
  rmst_Z_CH <- lapply(sim_data_CH, rmst_zs, t_cut=t_cut[i])
  rmst_Z_DE <- lapply(sim_data_DE, rmst_zs, t_cut=t_cut[i])
  rmst_Z_PH <- lapply(sim_data_PH, rmst_zs, t_cut=t_cut[i])
  
  
  # calculate power
  result$power_rmstDH[i] <- mean(rmst_Z_DH > qnorm(1-alpha)) 
  result$power_rmstCH[i] <- mean(rmst_Z_CH > qnorm(1-alpha)) 
  result$power_rmstDE[i] <- mean(rmst_Z_DE > qnorm(1-alpha)) 
  result$power_rmstPH[i] <- mean(rmst_Z_PH > qnorm(1-alpha)) 
  
  
}


#Figure3A

# simulate data

sim_data <- delayed_effect_sim(model = model_DE)
fit<- survfit(Surv(time, event) ~ group, data = sim_data)

#Calculate t_cut
data_minimax2 = sim_data %>% group_by(group) %>% summarize(maxTime= max(time)) 
data_minimax3 = data_minimax2 %>% summarize(minimax = min(maxTime))
t_cut_val = data_minimax3$minimax

#Calculate window WLR
rt = last(get_risk_table(sim_data)$t)

# Basic survival curves
p = ggsurvplot(fit, data = sim_data,
               size = 0.3,
               title = "A)",
               censor.size = 4,
               xlim = c(0,40),
               legend.title = " ",
               legend = c(0.8,0.9),
               legend.labs = c("Control", "Experimental")) 
p1 = p$plot + geom_segment(aes(x = rt, y = 0, xend = rt, yend = 0.5),color = "grey", size = 0.3) +
  geom_segment(aes(x = rt, y = 0.5, xend = 0, yend = 0.5), color = "grey", size = 0.3) +
  geom_segment(aes(x = t_cut_val, y = 0, xend = t_cut_val, yend = 0.7), color = "grey", size = 0.3) +
  geom_segment(aes(x = t_cut_val, y = 0.7, xend = 0, yend = 0.7), color = "grey", size = 0.3) +
  annotate("text", x = 15, y = 0.55, label = "Window log-rank") +
  annotate("text", x = 15, y = 0.75, label = "Window RMST")

p1

#Figure 3B

df = data.frame(time = rep(t_cut,4),
                scenario = c(rep("Proportional Hazards",length(t_cut)),rep("Delayed Effects",length(t_cut)),rep("Crossing Hazards",length(t_cut)),rep("Decreasing Effects",length(t_cut))),
                power = c(result$power_rmstPH,result$power_rmstDE,result$power_rmstCH,result$power_rmstDH))

df$scenario = factor(df$scenario,levels = c("Proportional Hazards","Delayed Effects","Decreasing Effects","Crossing Hazards"), labels = c("Proportional Hazards","Delayed Effects","Decreasing Effects","Crossing Hazards"),ordered = TRUE)

p2 <- ggplot(aes(x = time, y = power, color = scenario), data = df) +
  geom_line(size = 1.2) +
  scale_x_continuous(expression(paste("Cutting time ","(",t^"*",")"))) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"), limits = c(0,1)) +
  labs(y = "Power", color = "",title = "B)") +
  scale_color_manual(values=c("#F8766D", "#B79F00", "#00BA38", "#00BFC4")) +
  theme_light() +
  theme(text = element_text(size = 14),
        legend.position = c(0.17,0.85),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic", size = 8),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent"))

scenarios = grid.arrange(p1, p2,
                         nrow = 1, ncol = 2)
class(scenarios) <- c("scenarios", class(scenarios))

ggsave("/Users/jose/Dropbox/Research/2019/RMST as an alternative to HR in clinical trials with delayed effects/Journal of Applied Statistics - 1st Review/Figure_3.pdf",scenarios,
       width = 10, height = 4, dpi = 300)
