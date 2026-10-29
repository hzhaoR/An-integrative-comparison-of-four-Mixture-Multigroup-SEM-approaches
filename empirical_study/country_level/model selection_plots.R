
library(ggplot2)
library(ggpubr)

# Table
BICtable <- cbind(c(1:9), BIC_G[1:9], LL_S[1:9])
colnames(BICtable) <- c("Clusters", "BIC", "LogL")
BICtable <- as.data.frame(BICtable)
BICtable$Clusters <- as.integer(BICtable$Clusters)

# Table
AICtable <- cbind(c(1:9), AIC[1:9], LL_S[1:9])
colnames(AICtable) <- c("Clusters", "AIC", "LogL")
AICtable <- as.data.frame(AICtable)
AICtable$Clusters <- as.integer(AICtable$Clusters)

# xtable(BICtable)

# Table_CHull
CHull_table <- CHull(loglik = LL_S[1:9], nrpar = nr.par_obs[1:9], nsclust = c(1, 9))
colnames(CHull_table) <- c("Clusters", "LogL", "nr.par", "CHull")
CHull_table <- as.data.frame(CHull_table)
CHull_table$nr.par <- as.integer(CHull_table$nr.par)
# CHull_table <- nice_table(CHull_table)
# flextable::save_as_docx(CHull_table, path = "CHull_table.docx")


# Plot_BIC
j <- ggplot(data = as.data.frame(BICtable), aes(x = Clusters, y = BIC)) + geom_line() + 
  theme_classic() +
  geom_point() + scale_x_continuous(breaks = seq(1, 9)) + ggtitle("BIC per model") + 
  theme(text = element_text(size=10), legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 50, vjust = 0.5)) 
# +  scale_y_continuous(breaks = BICtable$BIC)
# ggsave("Fig6left.pdf", plot = last_plot(), width = 6.5, height = 7, dpi = 300, units = "in", device = cairo_pdf)
# dev.off()

# Plot_BIC
k <- ggplot(data = as.data.frame(AICtable), aes(x = Clusters, y = AIC)) + geom_line() + 
  theme_classic() +
  geom_point() + scale_x_continuous(breaks = seq(1, 9)) + ggtitle("AIC per model") + 
  theme(text = element_text(size=10), legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 50, vjust = 0.5)) 


# Plot_CHull
l <- ggplot(data = CHull_table, aes(x = nr.par, y = LogL)) + geom_line() + 
  theme_classic() +
  geom_point() + scale_x_continuous(breaks = nr.par_obs) + ggtitle("LogL per model") + 
  theme(text = element_text(size=10), legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
theme(axis.text.y = element_text(angle = 50, vjust = 0.5)) 

p1 <- ggarrange(plotlist = list(a, b, c), ncol = 3)
p2 <- ggarrange(plotlist = list(d, e, f), ncol = 3)
p3 <- ggarrange(plotlist = list(g, h, i), ncol = 3)
p4 <- ggarrange(plotlist = list(j, k, l), ncol = 3)

final_plot <- ggarrange(
  p1, p2, p3, p4,
  ncol = 1,
  nrow = 4
)

ggsave("Fig_MixSEM_final.pdf", plot = last_plot(), width = 9, height = 12, dpi = 300, units = "in", device = cairo_pdf)

