read_norm1 <- read.csv("/projectnb/bf528/users/group4/project3/samples/deseq2/deseq_norm_counts1.csv")
read_norm2 <- read.csv("/projectnb/bf528/users/group4/project3/samples/deseq2/deseq_norm_counts2.csv")
read_norm3 <- read.csv("/projectnb/bf528/users/group4/project3/samples/deseq2/deseq_norm_counts3.csv")

read_de1 <-read.csv("/projectnb/bf528/users/group4/project3/samples/deseq2/deseq_results1.csv")
read_de2 <-read.csv("/projectnb/bf528/users/group4/project3/samples/deseq2/deseq_results2.csv")
read_de3 <-read.csv("/projectnb/bf528/users/group4/project3/samples/deseq2/deseq_results3.csv")

de1 <- read_de1[order(read_de1$log2FoldChange, decreasing = TRUE),]
de2 <- read_de3[order(read_de2$log2FoldChange, decreasing = TRUE),]
de3 <- read_de2[order(read_de3$log2FoldChange, decreasing = TRUE),]

top_de1 = de1[1:50,]
top_de2 = de2[1:50,]
top_de3 = de3[1:50,]
gene_list = rbind(top_de1,top_de2,top_de3)

joined1 <- merge(read_norm1,gene_list,by="X")
joined2 <- merge(read_norm2,gene_list,by="X")
joined3 <- merge(read_norm3,gene_list,by="X")

norm_only1 = joined1[,1:3]
norm_only2 = joined2[,1:4]
norm_only3 = joined3[,1:4]

merge1 <- merge(norm_only1,norm_only2,by="X")
merge2 <- merge(merge1,norm_only3,by="X")

final <- merge2[,-1]
rownames(final) <- merge2[,1]
final <- as.matrix(final, rownames = 1)
heatmap(final)
