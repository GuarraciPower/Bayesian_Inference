if (!require(pacman)) install.packages("pacman")
p_load(rjags, coda, nimble, R2OpenBUGS, ggplot2, here, dplyr, ggpubr, tidyr)

# Load data
projdata <- as.data.frame(read.csv(here("data/projectdata.txt")))|>
  mutate(less_133pc = ifelse(Poverty == "<133% FPL", 1,0),
         btn133_400_pc = ifelse(Poverty == "133% to <400% FPL", 1,0),
         great_400pc = ifelse(Poverty == ">400% FPL", 1,0))

# Data prep for bugs
model_data <- list(
  Y = projdata$Vaccinated,
  N = projdata$Sample.Size,
  btn133_400_pc = projdata$btn133_400_pc,
  great_400pc = projdata$great_400pc,
  J = nrow(projdata)
)

# Initial values
model_inits <- list(
  list(beta0 = 0, beta1 = 0, beta2 = 0)
)

# Parameters to monitor
parameters <- c("beta0", "beta1", "beta2")

#Define the model
model1 <- function(){
  for (i in 1:J){
    Y[i] ~ dbin(p[i], N[i])
    logit(p[i]) <- beta0 + beta1*btn133_400_pc[i] + beta2*great_400pc[i]
  }
  #priors
  beta0 ~ dnorm(0, 0.001)
  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0, 0.001)
}

write.model(model1, here("models/model1code.txt"))
file.show(here("models/model1code.txt"))

model.out <- bugs(model_data, model_inits, 
                  parameters = parameters, model.file = here("models/model1code.txt"),
                  n.chains = 1, n.iter = 10000, n.burnin = 5000, codaPkg = TRUE, debug = TRUE)

out <- read.bugs(model.out)
summary(out)
plot(out)
mcmc_samples <- as.mcmc(out)
mcmc_df <- as.data.frame(mcmc_samples)
mcmc_df$iteration <- 1:nrow(mcmc_df)

HPDinterval(as.mcmc(as.matrix(out)))

# Create density plots using ggplot2

points_data <- data.frame(x = mcmc_df$beta0, x1 = mcmc_df$beta1,x2 = mcmc_df$beta2,
                          x3 = mcmc_df$deviance,y = rep(0, nrow(mcmc_df)))

densities<- ggarrange(p_beta0 <- ggplot(mcmc_df, aes(x = beta0)) +
  geom_density(fill = "blue", alpha = 0.1) +
  geom_point(data = points_data, aes(x = x, y = y)) +
  labs(title = "Posterior Distribution of beta0", x = "beta0", y = "Density")+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),

p_beta1 <- ggplot(mcmc_df, aes(x = beta1)) +
  geom_density(fill = "green", alpha = 0.1) +
  geom_point(data = points_data, aes(x = x1, y = y)) +
  labs(title = "Posterior Distribution of beta1", x = "beta1", y = "Density") +
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),

p_beta2 <- ggplot(mcmc_df, aes(x = beta2)) +
  geom_density(fill = "red", alpha = 0.1) +
  geom_point(data = points_data, aes(x = x2, y = y)) +
  labs(title = "Posterior Distribution of beta2", x = "beta2", y = "Density")+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),

p_dev <- ggplot(mcmc_df, aes(x = deviance)) +
  geom_density(fill = "yellow", alpha = 0.1) +
  geom_point(data = points_data, aes(x = x3, y = y)) +
  labs(title = "Posterior Distribution of deviance", x = "Deviance", y = "Density")+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), nrow = 2, ncol = 2)

ggsave(here("outputs/densityplots.png"), plot = densities, width = 12, height = 10)

## History plot $ posterior distributions

# Create trace plots

mcmc_long <- pivot_longer(mcmc_df, cols = -iteration, names_to = "Parameter", values_to = "Value")

trace_plot <- ggplot(mcmc_long, aes(x = iteration, y = Value, color = Parameter)) +
  geom_line() +
  scale_color_manual(values = c("blue", "green", "red", "yellow")) +
  facet_wrap(~ Parameter, scales = "free_y") +
  labs(title = "Trace Plots of MCMC Samples", x = "Iteration", y = "Parameter Value") +
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

# Save the trace plot
ggsave(here("outputs/traceplots.png"), plot = trace_plot, width = 12, height = 10)

par(mfrow=c(2,2))
traceplot(out)
densplot(out)

## Correlation and autocorrelation plots
crosscorr.plot(out)
autocorr.plot(out)
