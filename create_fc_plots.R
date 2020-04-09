library("ggplot2")
library("ggpubr")

# load limma data
leflu_data <- read.csv("limma_results_leflunomide.csv", as.is=TRUE)
flucon_data <- read.csv("limma_results_fluconazole.csv", as.is=TRUE)
ifos_data <- read.csv("limma_results_ifosfamide.csv", as.is=TRUE)

# filter for adj.P.Val < 0.05
leflu_significant <- leflu_data[leflu_data$adj.P.Val < 0.05,]
flucon_significant <- flucon_data[flucon_data$adj.P.Val < 0.05,]
ifos_significant <- ifos_data[ifos_data$adj.P.Val < 0.05,]

# create histograms
leflu_hist <- ggplot(leflu_significant, aes(x=logFC)) + geom_histogram(binwidth = 0.2, color="black", fill="#66e0ff") + labs(title="Leflunomide", x="Log Fold Change", y="count") + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=8), axis.title.y=element_text(size=8))
flucon_hist <- ggplot(flucon_significant, aes(x=logFC)) + geom_histogram(binwidth = 0.15, color="black", fill="#80ffbf") + labs(title="Fluconazole", x="Log Fold Change", y="count") + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=8), axis.title.y=element_text(size=8))
ifos_hist <- ggplot(ifos_significant, aes(x=logFC)) + geom_histogram(binwidth = 0.1, color="black", fill="#bb99ff") + labs(title="Ifosfamide", x="Log Fold Change", y="count") + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=8), axis.title.y=element_text(size=8))

# arrange histograms
all_hist <- ggarrange(leflu_hist, flucon_hist, ifos_hist, 
                      labels = c("A", "B", "C"), 
                      ncol = 3, nrow = 1)

# save to a png
ggsave("logfc_hist.png", plot = all_hist, width=6.5, height=2.5, units="in")


## create volcano plots using the unfiltered data

# create data for the probes i want to label
leflu_label <- leflu_data[(leflu_data$logFC > 1.5 | leflu_data$logFC < -1.5) & (leflu_data$P.Value < 0.05),]
flucon_label <- flucon_data[(flucon_data$logFC > 1.5 | flucon_data$logFC < -1.5) & (flucon_data$P.Value < 0.05),]
ifos_label <- ifos_data[(ifos_data$logFC > 1.5 | ifos_data$logFC < -1.5) & (ifos_data$P.Value < 0.05),]

# get top 3 ids to highlight
leflu_highlight <- leflu_data[which(leflu_data$X %in% c("1370269_at", "1392946_at", "1395403_at")),]
flucon_highlight <- flucon_data[which(flucon_data$X %in% c("1371076_at", "1387316_at", "1395403_at")),]
ifos_highlight <- ifos_data[which(ifos_data$X %in% c("1370355_at", "1369864_a_at", "1387874_at")),]

# plot settings for all 3 subplots
v_lines <- list(geom_vline(xintercept = 1.5, col = "red", linetype = "dotted", size = 0.5), geom_vline(xintercept = -1.5, col = "red", linetype = "dotted", size = 0.5), geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dotted", size = 0.5))
volcano_labels <- list(labs(x="Log Fold Change", y="-log10(Nominal P-Value)"))
volcano_theme <- list(theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=8), axis.title.y=element_text(size=8)))

# create plots for each chemical
leflu_scat <- ggplot(leflu_data, aes(x=logFC, y=-log10(P.Value))) + 
  geom_point(color="#66e0ff", size=0.5) + 
  geom_point(data=leflu_highlight, aes(x=logFC,y=-log10(P.Value)), color='red', size=0.6) +
  labs(title="Leflunomide") +
  geom_text(data = leflu_label, aes(x=logFC, y=-log10(P.Value), label=X), hjust = 0, nudge_x=0.25, size=1.2) + 
  expand_limits(x = c(-3.5, 9.2)) + 
  v_lines + volcano_labels + volcano_theme
flucon_scat <- ggplot(flucon_data, aes(x=logFC, y=-log10(P.Value))) + 
  geom_point(color="#80ffbf", size=0.5) + 
  geom_point(data=flucon_highlight, aes(x=logFC,y=-log10(P.Value)), color='red', size=0.6) +
  labs(title="Fluconazole") +
  geom_text(data = flucon_label, aes(x=logFC, y=-log10(P.Value), label=X), hjust = 0, nudge_x=0.25, size=1.2) + 
  expand_limits(x = c(-5, 5)) + 
  v_lines + volcano_labels + volcano_theme
ifos_scat <- ggplot(ifos_data, aes(x=logFC, y=-log10(P.Value))) + 
  geom_point(color="#bb99ff", size=0.5) + 
  geom_point(data=ifos_highlight, aes(x=logFC,y=-log10(P.Value)), color='red', size=0.6) +
  labs(title="Ifosfamide") +
  geom_text(data = ifos_label, aes(x=logFC, y=-log10(P.Value), label=X), hjust = 0, nudge_x=0.1, size=1.2) + 
  expand_limits(x = c(-3, 3)) + 
  v_lines + volcano_labels + volcano_theme

# arrange the plots on one graph
all_scat <- ggarrange(leflu_scat, flucon_scat, ifos_scat, 
                      labels = c("A", "B", "C"), 
                      ncol = 3, nrow = 1)

# save to a png
ggsave("logfc_scat_volcano_notsig.png", plot = all_scat, width=6.5, height=2.5, units="in")