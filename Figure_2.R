library(data.table)
library(dplyr)
library(ggplot2)

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
calculate_zs <- function(risk_table){
  
  #if (class(risk_table) == "list") risk_table = risk_table$risk_table
  
  n_e <- max(risk_table$n_e)
  n_c <- max(risk_table$n_c)
  
  
  # formulas for U and V[U] in Leton and Zuluaga pg. 596.
  
  u <- with(risk_table, sum(w * (d_c - d * n_c / n)))
  v_u <- with(risk_table, sum(w^2 * n_c * n_e * d * (n - d) / n / n / (n - 1), na.rm = TRUE))
  
  z_u <- u / sqrt(v_u)
  
  z_u
  
  
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
calculate_weights <- function(risk_table, method = "logrank", hr_fun = NULL, rho = 0, gamma = 0){
  
  if(method == "logrank")
  {
    risk_table$w <- 1
  }
  
  else if(method == "theta")
  {
    risk_table$w <- -log(hr_fun(risk_table$t))
  }
  
  else if(method == "hr_weight")
  {
    #risk_table$w <- (hr_fun(risk_table$t)^(-1) -hr_fun(max(risk_table$t))^(-1)) / (hr_fun(0)^(-1) -hr_fun(max(risk_table$t))^(-1))
    risk_table$w <- (hr_fun(risk_table$t)^(-1) - 1) / (hr_fun(0)^(-1) - 1)
    
  }
  
  else if(method == "fh"){
    
    # Kaplan-Meier estimate of survival in the lumped data:
    risk_table$s <- exp(cumsum(log(1 - risk_table$d / risk_table$n)))
    # Fleming-Harrington (rho, gamma) weights:
    risk_table$w <- risk_table$s ^ rho * (1 - risk_table$s) ^ gamma
    
  }
  
  
  
  risk_table
}

#################################################################################################################
n_c = 330/2
n_e = 330/2
rec_period = 12
rec_power = 1
med_c = 6
rate_e_1 = log(2)/6
rate_e_2 = log(2)/9
delay = 4
max_cal_t  = NULL
end_event = 200

# significance level (one-sided)
alpha <- 0.025

# number of simulations
M <- 10000

rho_vector = seq(0, by = 0.1, 1)
gamma_vector = seq(0, by = 0.1, 1)

result <- mat.or.vec(length(rho_vector),length(gamma_vector))

for(ii in 1:length(rho_vector)){
  for(jj in 1:length(gamma_vector)){
    
    cat("Rho = ", rho_vector[ii], "Gamma = ", gamma_vector[jj], "\n")
  
  # model parameters TS
  model_DE <- list(n_c = n_c,
                   n_e = n_e,
                   rec_period = rec_period, 
                   rec_power = rec_power, 
                   med_c = med_c, 
                   rate_e_1 = rate_e_1, 
                   rate_e_2 = rate_e_2, 
                   delay = delay, 
                   max_cal_t  = max_cal_t, 
                   n_events = end_event)
  
  # simulate data
  sim_data_DE <- replicate(M, delayed_effect_sim(model = model_DE), simplify = FALSE) 
  
  # create risk table
  risk_table_DE <- lapply(sim_data_DE, get_risk_table)
  
  # calculate standard logrank
  #logrank_Z_DE <- lapply(lapply(risk_table_DE, calculate_weights, method = "logrank"), calculate_zs)
  fh_Z_DE <- lapply(lapply(risk_table_DE, calculate_weights, method = "fh",rho=rho_vector[ii],gamma=gamma_vector[jj]), calculate_zs)
  
  
  # calculate power
  #result$power_logrank_DE[i] <- mean(logrank_Z_DE > qnorm(1-alpha))
  result[ii,jj] <- mean(fh_Z_DE > qnorm(1-alpha))
  
  }
}

save.image("~/LR_vs_FH_DE_10000.RData")

df_0 = mat.or.vec(length(rho_vector)*length(gamma_vector),3)

count = 1
for(i in 1:length(rho_vector)){
  for(j in 1:length(gamma_vector)){
    df_0[count,1] = rho_vector[i]
    df_0[count,2] = gamma_vector[j]
    df_0[count,3] = result[i,j]
    count = count + 1
  }
}

df = data.frame(rho = df_0[,1],
                gamma = df_0[,2],
                Power = df_0[,3])


p = ggplot(df, aes(rho, gamma, fill= Power)) + 
  geom_tile() +
  scale_fill_gradient(low="#F8766D", high="#619CFF") +
  theme_minimal() + 
  labs(y = expression(gamma), x = expression(rho), color = "Power") +
  theme(text = element_text(size = 12),
        legend.text = element_text(face = "italic"))

ggsave("Figure_2.pdf",p,
       width = 6, height = 3, dpi = 300)

