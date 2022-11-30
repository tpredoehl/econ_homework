###############################################################################
###
### set-up
###
###############################################################################
pkgs.installed <- installed.packages()
pkgs.required <- c(
  "car", 
  "AER", 
  "stargazer", 
  "tictoc", 
  "quantmod", 
  "xts", 
  "readr", 
  "latex2exp", 
  "gridExtra", 
  "summarytools", 
  "qwraps2", 
  "nortest", 
  "moments", 
  "xtable", 
  "sm", 
  "astsa", 
  "tseries", 
  "ggplot2", 
  "ggplotgui", 
  "shiny", 
  "tidyverse", 
  "gridExtra"
  )
pkgs.missing <- pkgs.required[which(!pkgs.required %in% pkgs.installed)]
lapply(pkgs.missing, install.packages, character.only = TRUE)
lapply(pkgs.required, require, character.only = TRUE)

rm(list=ls()) # clear all variables from enviroment/workspace

independent_AR_simulation<- function(phi_vec, c_vec, sigma_u_vec, T){
  # Simulate two independent RW and run spurious regression
  
  # function to simulate K independent AR process with :
  # T       = n.f simulated dates (simulated sample size)
  # phi_vec = [\phi_1, ..., phi_N]' = autoregressive coefficients for each AR process
  # c_vec   = [c_1, ..., c_N] = constant in AR process
  # number of AR(1) processes to simulate
  N <- length(phi_vec) 
  # T = T[6]
  # create empty matrix T X N that will contain all simulated AR processes
  AR_sim_mat <- NA*matrix(0,T+1,N) # empty vector ( (T+1) X 1 )
  u_sim_mat <- NA*matrix(0,T+1,N) # empty vector ( (T+1) X 1 )
  
  # loop over each AR process ind1
  for (ind1 in seq(1,N)) {
    # ind1 <- 1
    c       <-       c_vec[ind1]
    phi     <-     phi_vec[ind1]
    sigma_u <- sigma_u_vec[ind1]
    
    # simulate vector of innovations zero mean, and sd given as input
    # fix "seed" in order to produce same random numbers every time next line of code is called!
    # reference: https://stackoverflow.com/questions/13605271/reasons-for-using-the-set-seed-function
    as.numeric(Sys.time())-> t; set.seed((t - floor(t)) * 1e8 -> seed) 
    
    u                 <- rnorm(n = T+1, mean = 0, sd = sigma_u) #WN
    u_sim_mat[,ind1]  <- u
    
    # generate starting observation for the process: 
    # unconditional mean (stationary processes) or 0 (for non-stat processes only!, as phi=1 and c/(1-phi)=NA)
    if (abs(phi)<1){
      AR_sim_mat[1,ind1] <- c/(1-phi) ; #start from unconditional value of the process
    } else{
      AR_sim_mat[1,ind1] <- 0 ; #start from unconditional value of the process
    }
    
    # simulate all other observations t=2 to T+1
    for (t in seq(2,T+1)) AR_sim_mat[t,ind1] <- phi*AR_sim_mat[t-1,ind1] + u[t]
  }
  # cut the first row from the data
  # length(AR_sim_mat); nrow(AR_sim_mat)
  # length(u_sim_mat); nrow(u_sim_mat)
  AR_sim_mat <- AR_sim_mat[-1,]
  u_sim_mat <- u_sim_mat[-1,]
  
  # Output of function
  # 1 = simulated AR process
  # 2 = simulated innovations
  return(list(AR_sim_mat = AR_sim_mat, u_sim_mat = u_sim_mat))
}

###############################################################################
###
### run MC
###
###############################################################################
# N. of MC simulations
Nsim <- 100000
# Nsim <- 200

# DGP parameters : as above #############
phi_vec     = c(1, 1, 0.8, 0.8)
c_vec       = c(0,0, 0, 0)
sigma_u_vec = c(1, 1, 1, 1) 
T           = c(6, 12, 60, 120, 240, 360, 480)

# empty matrix to store intercept, t-stat(intercept), beta, t-stat (beta), and R^2
MC_mat <- NA*matrix(data = 0, nrow = length(T)*Nsim, ncol = 21)
colnames(MC_mat) <- c(
  "T", 
  "spur_RW_intercept", "spur_RW_t_intercept", "spur_RW_beta", "spur_RW_t_beta", "spur_RW_R2",
  "valid_RW_intercept", "valid_RW_t_intercept", "valid_RW_beta", "valid_RW_t_beta", "valid_RW_R2",
  "spur_AR_intercept", "spur_AR_t_intercept", "spur_AR_beta", "spur_AR_t_beta", "spur_AR_R2",
  "valid_AR_intercept", "valid_AR_t_intercept", "valid_AR_beta", "valid_AR_t_beta", "valid_AR_R2"
)
tic()
for (t in T) {
  #t = 6
  t_idx <- which(t == T)
  for(sim_idx in seq(1,Nsim)) {
    # sim_idx = 1
    output_temp <- independent_AR_simulation(phi_vec, c_vec, sigma_u_vec, t)
    AR_sim <- output_temp$AR_sim_mat
    reg_sim_RW_SPURS  <- lm(AR_sim[,2] ~ AR_sim[,1])
    reg_sim_RW_VALID  <- lm(AR_sim[,2] ~ AR_sim[,1] + lag(AR_sim[,2]))
    reg_sim_AR_SPURS  <- lm(AR_sim[,4] ~ AR_sim[,3])
    reg_sim_AR_VALID  <- lm(AR_sim[,4] ~ AR_sim[,3] + lag(AR_sim[,4]))
    mat_idx = (t_idx-1) * Nsim + sim_idx
    MC_mat[mat_idx, 1] <- t  # T
    MC_mat[mat_idx, 2] <- summary(reg_sim_RW_SPURS)$coefficients[1,1]  # intercept : value
    MC_mat[mat_idx, 3] <- summary(reg_sim_RW_SPURS)$coefficients[1,3]  # intercept : t-stat
    MC_mat[mat_idx, 4] <- summary(reg_sim_RW_SPURS)$coefficients[2,1]  # beta : value
    MC_mat[mat_idx, 5] <- summary(reg_sim_RW_SPURS)$coefficients[2,3]  # beta : t-stat
    MC_mat[mat_idx, 6] <- summary(reg_sim_RW_SPURS)$r.squared          # R^2
    
    MC_mat[mat_idx, 7] <- summary(reg_sim_RW_VALID)$coefficients[1,1]  # intercept : value
    MC_mat[mat_idx, 8] <- summary(reg_sim_RW_VALID)$coefficients[1,3]  # intercept : t-stat
    MC_mat[mat_idx, 9] <- summary(reg_sim_RW_VALID)$coefficients[2,1]  # beta : value
    MC_mat[mat_idx, 10] <- summary(reg_sim_RW_VALID)$coefficients[2,3]  # beta : t-stat
    MC_mat[mat_idx, 11] <- summary(reg_sim_RW_VALID)$r.squared          # R^2
    
    MC_mat[mat_idx, 12] <- summary(reg_sim_AR_SPURS)$coefficients[1,1]  # intercept : value
    MC_mat[mat_idx, 13] <- summary(reg_sim_AR_SPURS)$coefficients[1,3]  # intercept : t-stat
    MC_mat[mat_idx, 14] <- summary(reg_sim_AR_SPURS)$coefficients[2,1]  # beta : value
    MC_mat[mat_idx, 15] <- summary(reg_sim_AR_SPURS)$coefficients[2,3]  # beta : t-stat
    MC_mat[mat_idx, 16] <- summary(reg_sim_AR_SPURS)$r.squared          # R^2
    
    MC_mat[mat_idx, 17] <- summary(reg_sim_AR_VALID)$coefficients[1,1]  # intercept : value
    MC_mat[mat_idx, 18] <- summary(reg_sim_AR_VALID)$coefficients[1,3]  # intercept : t-stat
    MC_mat[mat_idx, 19] <- summary(reg_sim_AR_VALID)$coefficients[2,1]  # beta : value
    MC_mat[mat_idx, 20] <- summary(reg_sim_AR_VALID)$coefficients[2,3]  # beta : t-stat
    MC_mat[mat_idx, 21] <- summary(reg_sim_AR_VALID)$r.squared          # R^2
  }  
}
toc()

###############################################################################
###
### graph a.i
###
###############################################################################
MC_R2 <- MC_mat %>% 
  as.data.frame() %>% 
  as_tibble() %>%
  select(T, spur_RW_R2) %>% 
  mutate(T = as.factor(T))

graph1ai <- ggplot(MC_R2, aes(x = spur_RW_R2)) +
  geom_histogram(color = "#787878", fill = "#0099F8", position = 'stack', binwidth = 0.05, aes(y = after_stat(ndensity))) + 
  geom_density(aes(y = after_stat(ndensity)), lwd = .2,
               linetype = 1,
               colour = 2) +
  facet_grid( . ~ T, space = "free_x" ) +
  theme_bw() +
  ggtitle(label = TeX(r"(Distribution of $R^2$)")) + 
  theme(text = element_text()) +
  scale_y_continuous(
    breaks = c(.5, .10, .25, .5, .75, .90, .95),
    labels = scales::percent
  ) +
  scale_x_continuous(
    breaks = c(0, .5, 1),
    labels = scales::percent
  ) +
  ylab(label = "Quantile") + 
  xlab(label = TeX(r"($R^2$)"))

###############################################################################
###
### graph a.ii
###
###############################################################################
MC_beta_dist <- MC_mat %>% 
  as.data.frame() %>% 
  as_tibble() %>%
  select(T, spur_RW_t_beta) %>% 
  mutate(T = as.factor(T))

graph1aii <- ggplot(MC_beta_dist, aes(x = spur_RW_t_beta)) +
  geom_histogram(bins = 35, color = "#787878", fill = "#0099F8", position = 'stack', aes(y = after_stat(ndensity))) + 
  geom_density(aes(y = after_stat(ndensity)), lwd = .2,
               linetype = 1,
               colour = 2) +
  facet_grid( . ~ T, space = "free_x" ) +
  theme_bw() +
  ggtitle(label = TeX(r"(Distribution of $t_\beta$)")) + 
  theme(text = element_text()) +
  scale_y_continuous(
    breaks = c(.5, .10, .25, .5, .75, .90, .95),
    labels = scales::percent
  ) +
  scale_x_continuous(
    breaks = c(-50, 0, 50),
    labels = scales::number
  ) +
  ylab(label = "Quantile") + 
  xlab(label = TeX(r"($t_\beta$)"))

###############################################################################
###
### graph a.iii
###
###############################################################################
alpha <- 0.05;
crit_value <- qnorm(1-(alpha/2))

MC_beta_reject <- MC_mat %>% 
  as.data.frame() %>% 
  as_tibble() %>%
  select(T, spur_RW_t_beta, valid_RW_t_beta, spur_AR_t_beta, valid_AR_t_beta) %>% 
  mutate(
    t_crit = crit_value,
    reject_spur_RW = ifelse(abs(spur_RW_t_beta) > t_crit,1,0),
    reject_valid_RW = ifelse(abs(valid_RW_t_beta) > t_crit,1,0),
    reject_spur_AR = ifelse(abs(spur_AR_t_beta) > t_crit,1,0),
    reject_valid_AR = ifelse(abs(valid_AR_t_beta) > t_crit,1,0)
  ) %>% 
  group_by(T) %>% 
  summarise(
    pct_reject_spur_RW = sum(reject_spur_RW)/Nsim,
    pct_reject_valid_RW = sum(reject_valid_RW)/Nsim,
    pct_reject_spur_AR = sum(reject_spur_AR)/Nsim,
    pct_reject_valid_AR = sum(reject_valid_AR)/Nsim,
  )

graph1aiii <- ggplot(MC_beta_reject, aes(x = T)) +
  geom_line(linetype = "dotted", aes(y = pct_reject_spur_RW)) +
  geom_line(linetype = "dotted", aes(y = pct_reject_valid_RW)) +
  geom_line(linetype = "solid", aes(y = pct_reject_spur_AR)) +
  geom_line(linetype = "solid", aes(y = pct_reject_valid_AR)) +
  labs(
    title = TeX(r"(% of regressions which reject $H_0: \beta = 0$)"), 
    subtitle = TeX(r"(with $t = \frac{\beta - 0}{\sigma_{\beta}}>1.96$)")
  ) +
  geom_text(aes(x = 400, y = 0.80, label = "LRM1: spurious regression RW")) +
  geom_text(aes(x = 400, y = 0.30, label = "LRM2: spurious regression AR(1)")) +
  geom_text(aes(x = 400, y = 0.12, label = "LRM3: valid regression RW")) +
  geom_text(aes(x = 400, y = 0.03, label = "LRM4: valid regression AR(1)")) +
  ylab(label = "") +
  xlab(label = "T") +
  theme_bw() +
  theme(text = element_text(size = 8)) +
  scale_y_continuous(
    breaks = c(.05, .25, .5, .75, .95),
    labels = scales::percent
  )

###############################################################################
###
### save graphs
###
###############################################################################
ggsave(filename = "q1/1ai_chart.png", plot = graph1ai)
ggsave(filename = "q1/1aii_chart.png", plot = graph1aii)
ggsave(filename = "q1/1aiii_chart.png", plot = graph1aiii)
