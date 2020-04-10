library(dplyr)
library(ggplot2)
library(ggpubr)

# Read in the data files
# limma
leflu_limma <- read.csv('/projectnb/bf528/users/group4/project3/analysis/limma_results_leflunomide_new.csv', as.is = TRUE)
flucon_limma <- read.csv('/projectnb/bf528/users/group4/project3/analysis/limma_results_fluconazole_new.csv', as.is = TRUE)
ifos_limma <- read.csv('/projectnb/bf528/users/group4/project3/analysis/limma_results_ifosfamide_new.csv', as.is = TRUE)

# RNA-seq
#example_deseq <- read.csv("/project/bf528/project_3/results/example_deseq_results.csv", as.is = TRUE)
leflu_deseq <- read.csv('/projectnb/bf528/users/group4/project3/samples/deseq2/deseq_results1.csv', as.is = TRUE)
flucon_deseq <- read.csv('/projectnb/bf528/users/group4/project3/samples/deseq2/deseq_results2.csv', as.is = TRUE)
ifos_deseq <- read.csv('/projectnb/bf528/users/group4/project3/samples/deseq2/deseq_results3.csv', as.is = TRUE)

# map refseq to the row number
leflu_deseq$REFSEQID <- affy_map$REFSEQ[leflu_deseq$X]
flucon_deseq$REFSEQID <- affy_map$REFSEQ[flucon_deseq$X]
ifos_deseq$REFSEQID <- affy_map$REFSEQ[ifos_deseq$X]

# Load in the affy map
affy_map <- read.csv('/project/bf528/project_3/refseq_affy_map.csv', as.is=TRUE)

# Get the probe id from the limma results
#example_probe_ids <- example_limma$X
leflu_probe_ids <- leflu_limma$X
flucon_probe_ids <- flucon_limma$X
ifos_probe_ids <- ifos_limma$X

# Find probe id in affy map
# for (row in 1:length(example_probe_ids)){
#   refseq <- list(affy_map$REFSEQ[which(affy_map$PROBEID == example_probe_ids[row])])
#   example_limma$REFSEQID[row] <- refseq
# }
for (row in 1:length(leflu_probe_ids)){
  refseq <- list(affy_map$REFSEQ[which(affy_map$PROBEID == leflu_probe_ids[row])])
  leflu_limma$REFSEQID[row] <- refseq
}
for (row in 1:length(flucon_probe_ids)){
  refseq <- list(affy_map$REFSEQ[which(affy_map$PROBEID == flucon_probe_ids[row])])
  flucon_limma$REFSEQID[row] <- refseq
}
for (row in 1:length(ifos_probe_ids)){
  refseq <- list(affy_map$REFSEQ[which(affy_map$PROBEID == ifos_probe_ids[row])])
  ifos_limma$REFSEQID[row] <- refseq
}

# Expand rows that map to multiple refseqids and delete rows that map to no refseqids
# Subset to only probeid, refseqid, logFC, and aveexpr
# example_limma_mapped <- example_limma[FALSE, names(example_limma) %in% c("X", "REFSEQID", "logFC", "AveExpr")]
# as.data.frame(example_limma_mapped)
leflu_limma_mapped <- leflu_limma[FALSE, names(leflu_limma) %in% c("X", "REFSEQID", "logFC", "AveExpr")]
as.data.frame(leflu_limma_mapped)
flucon_limma_mapped <- flucon_limma[FALSE, names(flucon_limma) %in% c("X", "REFSEQID", "logFC", "AveExpr")]
as.data.frame(flucon_limma_mapped)
ifos_limma_mapped <- ifos_limma[FALSE, names(ifos_limma) %in% c("X", "REFSEQID", "logFC", "AveExpr")]
as.data.frame(ifos_limma_mapped)

# loop through
# for (row in 1:nrow(example_limma)) {
#   # filter for significant probes
#   if ((example_limma[row, "P.Value"] < 0.05) & (example_limma[row, "logFC"] < -1.5 | example_limma[row, "logFC"] > 1.5)) {
#     probe<- example_limma[row,"X"] 
#     refseq <- example_limma[row, "REFSEQID"] 
#     logfc <- example_limma[row, "logFC"]
#     aveexpr <- example_limma[row, "AveExpr"]
#     # if there is a match
#     if (length(refseq[[1]])>0) {
#       # taking into account multiple refseqids
#       for (id in refseq[[1]]) {
#         id_row <- data.frame(probe,id, logfc, aveexpr) 
#         names(id_row)<-c("X","REFSEQID", "logFC", "AveExpr") 
#         # add a new row to the dataframe
#         example_limma_mapped <- rbind(example_limma_mapped, id_row)
#       }
#     }
#   }
# }
for (row in 1:nrow(leflu_limma)) {
  # filter for significant probes
  if (leflu_limma[row, "P.Value"] < 0.05 & (leflu_limma[row, "logFC"] < -1.5 | leflu_limma[row, "logFC"] > 1.5)) {
    probe<- leflu_limma[row,"X"] 
    refseq <- leflu_limma[row, "REFSEQID"] 
    logfc <- leflu_limma[row, "logFC"]
    aveexpr <- leflu_limma[row, "AveExpr"]
    # if there is a match
    if (length(refseq[[1]])>0) {
      # taking into account multiple refseqids
      for (id in refseq[[1]]) {
        id_row <- data.frame(probe,id, logfc, aveexpr) 
        names(id_row)<-c("X","REFSEQID", "logFC", "AveExpr") 
        # add a new row to the dataframe
        leflu_limma_mapped <- rbind(leflu_limma_mapped, id_row)
      }
    }
  }
}
for (row in 1:nrow(flucon_limma)) {
  # filter for significant probes
  if (flucon_limma[row, "P.Value"] < 0.05 & (flucon_limma[row, "logFC"] < -1.5 | flucon_limma[row, "logFC"] > 1.5)) {
    probe<- flucon_limma[row,"X"] 
    refseq <- flucon_limma[row, "REFSEQID"] 
    logfc <- flucon_limma[row, "logFC"]
    aveexpr <- flucon_limma[row, "AveExpr"]
    # if there is a match
    if (length(refseq[[1]])>0) {
      # taking into account multiple refseqids
      for (id in refseq[[1]]) {
        id_row <- data.frame(probe,id, logfc, aveexpr) 
        names(id_row)<-c("X","REFSEQID", "logFC", "AveExpr") 
        # add a new row to the dataframe
        flucon_limma_mapped <- rbind(flucon_limma_mapped, id_row)
      }
    }
  }
}
for (row in 1:nrow(ifos_limma)) {
  # filter for significant probes
  if (ifos_limma[row, "P.Value"] < 0.05 & (ifos_limma[row, "logFC"] < -1.5 | ifos_limma[row, "logFC"] > 1.5)) {
    probe<- ifos_limma[row,"X"] 
    refseq <- ifos_limma[row, "REFSEQID"] 
    logfc <- ifos_limma[row, "logFC"]
    aveexpr <- ifos_limma[row, "AveExpr"]
    # if there is a match
    if (length(refseq[[1]])>0) {
      # taking into account multiple refseqids
      for (id in refseq[[1]]) {
        id_row <- data.frame(probe,id, logfc, aveexpr) 
        names(id_row)<-c("X","REFSEQID", "logFC", "AveExpr") 
        # add a new row to the dataframe
        ifos_limma_mapped <- rbind(ifos_limma_mapped, id_row)
      }
    }
  }
}


# collapsed the mapped list based on refseqid, taking the median of the logFC
#example_limma_mapped.filtered <- example_limma_mapped %>% group_by(REFSEQID) %>% summarise(logFC=median(logFC), AveExpr=median(AveExpr))
leflu_limma_mapped.filtered <- leflu_limma_mapped %>% group_by(REFSEQID) %>% summarise(logFC=median(logFC), AveExpr=median(AveExpr))
flucon_limma_mapped.filtered <- flucon_limma_mapped %>% group_by(REFSEQID) %>% summarise(logFC=median(logFC), AveExpr=median(AveExpr))
ifos_limma_mapped.filtered <- ifos_limma_mapped %>% group_by(REFSEQID) %>% summarise(logFC=median(logFC), AveExpr=median(AveExpr))


# Get only the significant DEGs from RNAseq
#example_deseq <- subset(example_deseq, (pvalue < 0.05) & (log2FoldChange < -1.5 | log2FoldChange > 1.5))
leflu_deseq <- subset(leflu_deseq, pvalue < 0.05 & (log2FoldChange < -1.5 | log2FoldChange > 1.5))
flucon_deseq <- subset(flucon_deseq, pvalue < 0.05 & (log2FoldChange < -1.5 | log2FoldChange > 1.5))
ifos_deseq <- subset(ifos_deseq, pvalue < 0.05 & (log2FoldChange < -1.5 | log2FoldChange > 1.5))

# Assuming there are no replicates...
# See how many of the filtered map are in the deseq results
# count <- 0
# for (row in 1:nrow(example_deseq)) {
#   id_deseq <- example_deseq[row, "X"]
#   logfc <- example_deseq[row, "log2FoldChange"]
#   limma_row <- example_limma_mapped.filtered[which(example_limma_mapped.filtered$REFSEQID == id_deseq), ]
#   if (nrow(limma_row) != 0) {
#     if(sign(logfc) == sign(limma_row[["logFC"]])) {
#       count <- count+1
#     }
#   }
# }
# example_concordance <- (2 * count)/(nrow(example_deseq) + nrow(example_limma_mapped.filtered))
count <- 0
for (row in 1:nrow(leflu_deseq)) {
  id_deseq <- leflu_deseq[row, "X"]
  logfc <- leflu_deseq[row, "log2FoldChange"]
  limma_row <- leflu_limma_mapped.filtered[which(leflu_limma_mapped.filtered$REFSEQID == id_deseq), ]
  if (nrow(limma_row) != 0) {
    if(sign(logfc) == sign(limma_row[["logFC"]])) {
      count <- count+1
    }
  }
}
leflu_concordance <- (2 * count)/(nrow(leflu_deseq) + nrow(leflu_limma_mapped.filtered))
count <- 0
for (row in 1:nrow(flucon_deseq)) {
  id_deseq <- flucon_deseq[row, "X"]
  logfc <- flucon_deseq[row, "log2FoldChange"]
  limma_row <- flucon_limma_mapped.filtered[which(flucon_limma_mapped.filtered$REFSEQID == id_deseq), ]
  if (nrow(limma_row) != 0) {
    if(sign(logfc) == sign(limma_row[["logFC"]])) {
      count <- count+1
    }
  }
}
flucon_concordance <- (2 * count)/(nrow(flucon_deseq) + nrow(flucon_limma_mapped.filtered))
count <- 0
for (row in 1:nrow(ifos_deseq)) {
  id_deseq <- ifos_deseq[row, "X"]
  logfc <- ifos_deseq[row, "log2FoldChange"]
  limma_row <- ifos_limma_mapped.filtered[which(ifos_limma_mapped.filtered$REFSEQID == id_deseq), ]
  if (nrow(limma_row) != 0) {
    if(sign(logfc) == sign(limma_row[["logFC"]])) {
      count <- count+1
    }
  }
}
ifos_concordance <- (2 * count)/(nrow(ifos_deseq) + nrow(ifos_limma_mapped.filtered))


# Make scatter plot with treatment on x-axis and concordance of y-axis (like in figure 2a)
# Make a table of the number of DEGs and concordance
example_scatter_deseq <- data.frame("Concordance" = example_concordance, "Effect"=nrow(example_deseq))
ggplot(example_scatter_deseq, aes(x=Effect, y=Concordance)) +
  geom_point()

scatter_lables <- c("LEF", "FLU", "IFO")
scatter_deseq <- data.frame("Concordance"=c(leflu_concordance, flucon_concordance, ifos_concordance), "Effect"=c(nrow(leflu_deseq), nrow(flucon_deseq), nrow(ifos_deseq)))
scatter_deseq_plot <- ggplot(scatter_deseq, aes(x=Effect, y=Concordance)) +
  geom_point() +
  labs(x="RNA-Seq Effect Size", y="Concordance") + 
  geom_text(data=scatter_labels, aes(x=logFC, y=-log10(P.Value), label=X), hjust = 0, nudge_x=0.1, size=1.2)

# Make a table for microarrays
example_scatter_limma <- data.frame("Concordance"=example_concordance, "Effect"=nrow(example_limma_mapped.filtered))
ggplot(example_scatter_limma, aes(x=Effect, y=Concordance)) +
  geom_point() +

scatter_limma <- data.frame("Concordance"=c(leflu_concordance, flucon_concordance, ifos_concordance), "Effect"=c(nrow(leflu_limma), nrow(flucon_limma), nrow(ifos_limma)))
scatter_limma_plot <- ggplot(scatter_limma, aes(x=Effect, y=Concordance)) +
  geom_point() +
  labs(x="Microarray Effect Size", y="Concordance") + 
  geom_text(data=scatter_labels, aes(x=logFC, y=-log10(P.Value), label=X), hjust = 0, nudge_x=0.1, size=1.2)

scatter <- ggarrange(scatter_limma_plot, scatter_deseq_plot, 
                      labels = c("A", "B"),
                      ncol = 2, nrow = 1)

# Find median of the base mean from deseq results
example_deseq_med <- median(example_deseq$baseMean)
example_deseq_above <- subset(example_deseq, baseMean > example_deseq_med) # genes above median
example_deseq_below <- subset(example_deseq, baseMean < example_deseq_med) # genes below median

leflu_deseq_med <- median(leflu_deseq$baseMean)
leflu_deseq_above <- subset(leflu_deseq, baseMean > leflu_deseq_med) # genes above median
leflu_deseq_below <- subset(leflu_deseq, baseMean < leflu_deseq_med) # genes below median

flucon_deseq_med <- median(flucon_deseq$baseMean)
flucon_deseq_above <- subset(flucon_deseq, baseMean > flucon_deseq_med) # genes above median
flucon_deseq_below <- subset(flucon_deseq, baseMean < flucon_deseq_med) # genes below median

ifos_deseq_med <- median(ifos_deseq$baseMean)
ifos_deseq_above <- subset(ifos_deseq, baseMean > ifos_deseq_med) # genes above median
ifos_deseq_below <- subset(ifos_deseq, baseMean < ifos_deseq_med) # genes below median

# find median of the average expression from limma results
example_limma_med <- median(example_limma_mapped.filtered$AveExpr)
example_limma_above <- subset(example_limma_mapped.filtered, AveExpr > example_limma_med)
example_limma_below <- subset(example_limma_mapped.filtered, AveExpr < example_limma_med)

leflu_limma_med <- median(leflu_limma_mapped.filtered$AveExpr)
leflu_limma_above <- subset(leflu_limma_mapped.filtered, AveExpr > leflu_limma_med)
leflu_limma_below <- subset(leflu_limma_mapped.filtered, AveExpr < leflu_limma_med)

flucon_limma_med <- median(flucon_limma_mapped.filtered$AveExpr)
flucon_limma_above <- subset(flucon_limma_mapped.filtered, AveExpr > flucon_limma_med)
flucon_limma_below <- subset(flucon_limma_mapped.filtered, AveExpr < flucon_limma_med)

ifos_limma_med <- median(ifos_limma_mapped.filtered$AveExpr)
ifos_limma_above <- subset(ifos_limma_mapped.filtered, AveExpr > ifos_limma_med)
ifos_limma_below <- subset(ifos_limma_mapped.filtered, AveExpr < ifos_limma_med)

# Recompute concordance for each of the above and below groups
# Concordance for below groups
count <- 0
for (row in 1:nrow(example_deseq_below)) {
  id_deseq <- example_deseq_below[row, "X"]
  logfc <- example_deseq_below[row, "log2FoldChange"]
  below_row <- example_limma_below[which(example_limma_below$REFSEQID == id_deseq), ]
  if (nrow(below_row) != 0) {
    if(sign(logfc) == sign(below_row[["logFC"]])) {
      count <- count+1
    }
  }
}
example_below_concord <- (2*count)/(nrow(example_deseq_below) + nrow(example_limma_below))
count <- 0
for (row in 1:nrow(leflu_deseq_below)) {
  id_deseq <- leflu_deseq_below[row, "X"]
  logfc <- leflu_deseq_below[row, "log2FoldChange"]
  below_row <- leflu_limma_below[which(leflu_limma_below$REFSEQID == id_deseq), ]
  if (nrow(below_row) != 0) {
    if(sign(logfc) == sign(below_row[["logFC"]])) {
      count <- count+1
    }
  }
}
leflu_below_concord <- (2*count)/(nrow(leflu_deseq_below) + nrow(leflu_limma_below))
count <- 0
for (row in 1:nrow(flucon_deseq_below)) {
  id_deseq <- flucon_deseq_below[row, "X"]
  logfc <- flucon_deseq_below[row, "log2FoldChange"]
  below_row <- flucon_limma_below[which(flucon_limma_below$REFSEQID == id_deseq), ]
  if (nrow(below_row) != 0) {
    if(sign(logfc) == sign(below_row[["logFC"]])) {
      count <- count+1
    }
  }
}
flucon_below_concord <- (2*count)/(nrow(flucon_deseq_below) + nrow(flucon_limma_below))
count <- 0
for (row in 1:nrow(ifos_deseq_below)) {
  id_deseq <- ifos_deseq_below[row, "X"]
  logfc <- ifos_deseq_below[row, "log2FoldChange"]
  below_row <- ifos_limma_below[which(ifos_limma_below$REFSEQID == id_deseq), ]
  if (nrow(below_row) != 0) {
    if(sign(logfc) == sign(below_row[["logFC"]])) {
      count <- count+1
    }
  }
}
ifos_below_concord <- (2*count)/(nrow(ifos_deseq_below) + nrow(ifos_limma_below))

# Concordance for above groups
count <- 0
for (row in 1:nrow(leflu_deseq_above)) {
  id_deseq <- leflu_deseq_above[row, "X"]
  logfc <- leflu_deseq_above[row, "log2FoldChange"]
  above_row <- leflu_limma_above[which(leflu_limma_above$REFSEQID == id_deseq), ]
  if (nrow(above_row) != 0) {
    if(sign(logfc) == sign(above_row[["logFC"]])) {
      count <- count+1
    }
  }
}
leflu_above_concord <- (2*count)/(nrow(leflu_deseq_above) + nrow(leflu_limma_above))
count <- 0
for (row in 1:nrow(flucon_deseq_above)) {
  id_deseq <- flucon_deseq_above[row, "X"]
  logfc <- flucon_deseq_above[row, "log2FoldChange"]
  above_row <- flucon_limma_above[which(flucon_limma_above$REFSEQID == id_deseq), ]
  if (nrow(above_row) != 0) {
    if(sign(logfc) == sign(above_row[["logFC"]])) {
      count <- count+1
    }
  }
}
flucon_above_concord <- (2*count)/(nrow(flucon_deseq_above) + nrow(flucon_limma_above))
count <- 0
for (row in 1:nrow(ifos_deseq_above)) {
  id_deseq <- ifos_deseq_above[row, "X"]
  logfc <- ifos_deseq_above[row, "log2FoldChange"]
  above_row <- ifos_limma_above[which(ifos_limma_above$REFSEQID == id_deseq), ]
  if (nrow(above_row) != 0) {
    if(sign(logfc) == sign(above_row[["logFC"]])) {
      count <- count+1
    }
  }
}
ifos_above_concord <- (2*count)/(nrow(ifos_deseq_above) + nrow(ifos_limma_above))

# Make a barplot of all of these values
# make into a table
bar <- data.frame("Concordance"=c(leflu_above_concord, leflu_below_concord, leflu_concordance,
                                  flucon_above_concord, flucon_below_concord, flucon_concordance,
                                  ifos_above_concord, ifos_below_concord, ifos_concordance), 
                  "Analysis"=c("Above", "Below", "Overall", "Above", "Below", "Overall", "Above", "Below", "Overall"),
                  "Chemical"=c("Leflunomide","Leflunomide","Leflunomide",
                               "Fluconazole","Fluconazole","Fluconazole",
                               "Ifosfamide","Ifosfamide","Ifosfamide"))
# Plot data with a side-by-side bar chart
ggplot(bar, aes(fill=Analysis, x=Chemical, y=Concordance)) + 
  geom_bar(stat="identity", position="dodge") +
  labs(title="Above and Below Median Concordance")





