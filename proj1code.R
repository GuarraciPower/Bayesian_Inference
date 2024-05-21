if (!require(pacman)) install.packages("pacman")
p_load(rjags, coda, nimble, R2OpenBUGS, ggplot2, here)

# Load data
projdata <- as.data.frame(read.csv(here("data/projectdata.txt")))
