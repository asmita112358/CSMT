library(ggplot2)
library(patchwork)
library(ggpubr)
library(tidyverse)
library(data.table)
library(ggpattern)
setwd("/Users/asmitaroy/Downloads/CSMT/sim_outputs")
##Size and power plots for the relationship between n and K
p1 <- read.csv("power150.csv")[,-1]
p2 <- read.csv("power250.csv")[,-1]
p3 <- read.csv("power350.csv")[,-1]
p4 <- read.csv("power450.csv")[,-1]
p5 <- read.csv("power550.csv")[,-1]

# Define a lookup table for sample_size and corresponding K values

custom_labeller <- function(df) {
  paste0("n = ", df$n, ", K = ", df$K)
}


n = rep(seq(150,550, by = 100), each = 5)
K = rep(c(6, 7, 9, 10, 11),each = 5)
sample_size = custom_labeller(data.frame(n,K))
pow_data <- data.frame(sample_size, rbind(p1, p2, p3, p4, p5))
colnames(pow_data)[2] = "Signal strength"
colnames(pow_data)[7:10] = c("Subbel_se", "ABtest_se", "MaxP_se", "Sobel_se")


# Reshape data to long format for easier plotting
pow_data_long <- pow_data %>%
  pivot_longer(cols = starts_with("Subbel"):starts_with("Sobel"),
               names_to = "Method",
               values_to = "Value") %>%
  mutate(Method = factor(Method, levels = c("Subbel", "ABtest", "MaxP", "Sobel")))
ylims = range(pow_data_long$Value)
power = ggplot(pow_data_long, aes(x = `Signal strength`, y = Value, color = Method)) +
  geom_line() +
  facet_wrap(~ sample_size,nrow = 1, ncol = 5, scales = "free_y") + 
  ylim(ylims)+# One grid for each sample size
  labs(x = "Signal Strength", y = "Power",
       title = "Empirical power by Signal strength for different sample sizes") +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "lightgrey", color = "black"),
        panel.border = element_rect(color = "black", fill = NA))


s = list()
d = list()

s[[1]] <- read.csv("sparse150.csv")[,-1]
s[[2]] <- read.csv("sparse250.csv")[,-1]
s[[3]] <- read.csv("sparse350.csv")[,-1]
s[[4]] <- read.csv("sparse450.csv")[,-1]
s[[5]] <- read.csv("sparse550.csv")[,-1]

d[[1]] <- read.csv("dense150.csv")[,-1]
d[[2]] <- read.csv("dense250.csv")[,-1]
d[[3]] <- read.csv("dense350.csv")[,-1]
d[[4]] <- read.csv("dense450.csv")[,-1]
d[[5]] <- read.csv("dense550.csv")[,-1]



N.rep = 500
size_sparse = do.call(rbind,lapply(s, function(x)colMeans(x <= 0.05)))
se_size_sparse = sqrt(size_sparse*(1-size_sparse)/N.rep)
colnames(size_sparse) <- colnames(se_size_sparse) <-  c("Subbel", "ABtest", "MaxP", "Sobel")


size_dense = do.call(rbind,lapply(d, function(x)colMeans(x <= 0.05)))
se_size_dense = sqrt(size_dense*(1-size_dense)/N.rep)
colnames(size_dense) <- colnames(se_size_dense) <-  c("Subbel", "ABtest", "MaxP", "Sobel")
Null_type <- rep(c("Sparse null", "Dense null"), each = 5)
n = seq(150,550, by = 100)
K = c(6, 7, 9, 10, 11)
sample_size = rep(custom_labeller(data.frame(n,K)), times = 2)
df_size = data.frame(sample_size, Null_type, rbind(size_sparse, size_dense))
df_size_se = data.frame(sample_size, Null_type, rbind(se_size_sparse, se_size_dense))


df_size_long <- df_size %>%
  pivot_longer(cols = starts_with("Subbel"):starts_with("Sobel"),
               names_to = "Method",
               values_to = "Size") %>%
  mutate(Method = factor(Method, levels = c("Subbel", "ABtest", "MaxP", "Sobel")))

df_se_size_long <- df_size_se %>% 
  pivot_longer(cols = starts_with("Subbel"):starts_with("Sobel"),
               names_to = "Method",
               values_to = "se_size") %>%
  mutate(Method = factor(Method, levels = c("Subbel", "ABtest", "MaxP", "Sobel")))

df_size_long$se = df_se_size_long$se_size

size = ggplot(df_size_long, aes(x = Method, y = Size, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + 
  scale_pattern_manual(
    "Method",
    values = c("none", "stripe", "none", "none")
  )+
  geom_col_pattern(
    aes(pattern = Method),
    position = "dodge",
    pattern_angle = 45,
    pattern_density = .1,
    pattern_spacing = .04,
    pattern_fill = "black",
    color = "black"
  ) +
  geom_errorbar(aes(ymin = Size - 1.96 * se, ymax = Size + 1.96 * se), 
                width = 0.2, position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0.05)+
  facet_grid(rows = vars(Null_type), cols = vars(sample_size)) +
  theme_minimal() +
  labs(x = "Method", y = "Size", title = "Empirical size by method for different sample sizes") +
  theme(legend.position = "bottom", 
        strip.text.x = element_text(size = 12, face = "bold"),  # Style for Sample_size facet labels
        strip.text.y = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "lightgrey", color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5) )
allp <- ggarrange(size, power, nrow = 2, heights = c(1.3,1))
ggsave("var_K2.png", allp, device = png, width = 9, height = 7.5, units = "in", dpi = 1200, bg = "white")




##All power plots

p2 <- read.csv("power100.csv")[,-1]
p1 <- read.csv("power200.csv")[,-1]
p3 <- read.csv("power300.csv")[,-1]

colnames(p1) <- colnames(p2) <- colnames(p3) <- c("Signal strength", "Subbel", "ABtest", "MaxP", "Sobel")
n = rep(seq(100,300, by = 100), each = 9)
K = rep(c(5,7,8),each = 9)
sample_size = custom_labeller(data.frame(n,K))
pow_data1 <- data.frame(sample_size, rbind(p1, p2, p3))

pv1<- read.csv("power100_var.csv")[,-1]
pv2 <- read.csv("power200_var.csv")[,-1]
pv3 <- read.csv("power300_var.csv")[,-1]
colnames(pv1) <- colnames(pv2) <- colnames(pv3) <- c("Signal ratio", "Subbel", "ABtest", "MaxP", "Sobel")
n = rep(seq(100,300, by = 100), each = 7)
K = rep(c(5,7,8),each = 7)
sample_size = custom_labeller(data.frame(n,K))

pow_data2 <- data.frame(sample_size, rbind(pv1, pv2, pv3))

# Reshape data to long format for easier plotting
pow_data_long <- pow_data1 %>%
  pivot_longer(cols = starts_with("Subbel"):starts_with("Sobel"),
               names_to = "Method",
               values_to = "Value") %>%
  mutate(Method = factor(Method, levels = c("Subbel", "ABtest", "MaxP", "Sobel")))
ylims = range(pow_data_long$Value)
power1 = ggplot(pow_data_long, aes(x = `Signal.strength`, y = Value, color = Method)) +
  geom_line() +
  facet_wrap(~ sample_size,nrow = 1, ncol = 5, scales = "free_y") + 
  ylim(ylims)+# One grid for each sample size
  labs(x = expression(paste("Signal Strength = ", alpha, " = ",beta)), y = "Power",
       title = "Empirical power by Signal strength for different sample sizes") +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "lightgrey", color = "black"),
        panel.border = element_rect(color = "black", fill = NA))

pow_data_long <- pow_data2 %>%
  pivot_longer(cols = starts_with("Subbel"):starts_with("Sobel"),
               names_to = "Method",
               values_to = "Value") %>%
  mutate(Method = factor(Method, levels = c("Subbel", "ABtest", "MaxP", "Sobel")))
ylims = range(pow_data_long$Value)
power2 = ggplot(pow_data_long, aes(x = `Signal.ratio`, y = Value, color = Method)) +
  geom_line() +
  facet_wrap(~ sample_size,nrow = 1, ncol = 5, scales = "free_y") + 
  ylim(ylims)+# One grid for each sample size
  labs(x = expression(paste("Signal Ratio = ",alpha,"/",beta)), y = "Power",
       title = "Empirical power by Signal ratio for different sample sizes") +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "lightgrey", color = "black"),
        panel.border = element_rect(color = "black", fill = NA))

allpow <- ggarrange(power1,power2, nrow = 2, common.legend = TRUE)
ggsave("all_power.png", allpow, device = png, width = 8, height = 6, units = "in", dpi = 1200, bg = "white")

