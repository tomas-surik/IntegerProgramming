"""
VISUALISE SENSITIVITY ANALYSIS

This code visualises the results of optimisation and sensitivity analysis.
"""
library(tidyverse)
library(gridExtra)

result_AB <- read_csv("C:/Data/result_AB.csv")
result_AB <- subset(result_AB, result_AB$A_cost >= 300)
result_AB <- subset(result_AB, result_AB$A_cost <= 500)
result_AB <- subset(result_AB, result_AB$B_cost >= 500)
result_AB <- subset(result_AB, result_AB$B_cost <= 700)

result_AB$solution <- as.factor(result_AB$options_selected)
result_AB$optimal <- ifelse(result_AB$options_selected == "[[1, 7, 8]]", "Yes","No")
levels(result_AB$solution) <-  c("[1, 7, 8]","[1, 8, 10]","[2, 8, 10]")

ggplot(aes(x = A_cost, y = B_cost, color = solution), data = result_AB) + geom_point()
p <- ggplot(aes(x = A_cost, y = B_cost, color = solution), data = result_AB) + geom_point(shape = 19) + scale_color_brewer(palette="Dark2")
p <- p + ylim(c(500, 700)) + xlim(300, 500) + geom_hline(yintercept = 600, color = "red") + geom_vline(xintercept = 400, color = "red")
p <- p + theme_bw() + xlab("A Cost ($)") + ylab("B Cost ($)") + labs(color = "Options Selected")
p <- p + scale_color_manual(values = c("orange", "blue", "forestgreen"))
p <- p + labs(title="Sensitivity to A Cost and B Cost")
p


result_C <- read_csv("C:/Data/result_C.csv")
result_C <- subset(result_C, result_C$C_cost >= 2000)
result_C <- subset(result_C, result_C$C_cost <= 4000)
result_C <- subset(result_C, result_C$C_base <= 4)

result_C$solution <- as.factor(result_C$options_selected)
levels(result_C$solution) <- c("[1, 7, 8]","[1, 8, 10]","[7, 8, 10]") 
result_C$C_base <- as.factor(result_C$C_base)

q <- ggplot(aes(x = C_base, y = C_cost, color = solution), data = result_C) + geom_point(shape = 15)
q <- q + theme_bw() + xlab("C Base") + ylab("C Cost ($)")
q <- q + scale_color_manual(values = c("orange", "blue", "lightblue"))
q <- q +  geom_hline(yintercept = 3000, color = "red") + geom_vline(xintercept = 2, color = "red")
q <- q + labs(color = "Options Selected") + labs(title="Sensitivity to C Base and C Cost") 
q

grid.arrange(p, q, ncol=2)
