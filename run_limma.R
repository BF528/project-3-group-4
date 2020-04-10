# try re-doing limma for just the matching controls

library("limma")

# sample info dataframe with array_id and chemical columns
samples <- read.csv('/project/bf528/project_3/groups/group_3_mic_info.csv',as.is=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt',
                  sep='\t',
                  as.is=TRUE,
                  header=TRUE,
                  row.names=1,
)

## LEFLUNOMIDE
# subset the full expression matrix to just those in this comparison
rma.subset.leflu <- rma[paste0('X',samples$array_id[samples$chemical=='LEFLUNOMIDE' | (samples$chemical=='Control' & samples$vehicle=='CORN_OIL_100_%')])]

# construct a design matrix modeling treatment vs control for use by limma
design <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle=='CORN_OIL_100_%'],
    levels=c('Control','LEFLUNOMIDE')
  )
)
colnames(design) <- c('Intercept','LEFLUNOMIDE')

# run limma
fit <- lmFit(rma.subset.leflu, design)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset.leflu), adjust='BH')

# filter for unadjusted p-value < 0.05
t.filtered <- subset(t, adj.P.Val < 0.05)
t.filtered <- subset(t, P.Value < 0.05)
# write out the results to file
write.csv(t,'limma_results_leflunomide_new.csv')



## FLUCONAZOLE
# subset the full expression matrix to just those in this comparison
rma.subset <- rma[paste0('X',samples$array_id[samples$chemical=='FLUCONAZOLE' | (samples$chemical=='Control' & samples$vehicle=='CORN_OIL_100_%')])]

# construct a design matrix modeling treatment vs control for use by limma
design <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle=='CORN_OIL_100_%'],
    levels=c('Control','FLUCONAZOLE')
  )
)
colnames(design) <- c('Intercept','FLUCONAZOLE')

# run limma
fit <- lmFit(rma.subset, design)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')

# filter for unadjusted p-value < 0.05
t.filtered <- subset(t, adj.P.Val < 0.05)
t.filtered <- subset(t, P.Value < 0.05)
# write out the results to file
write.csv(t,'limma_results_fluconazole_new.csv')


## IFOSFAMIDE
# subset the full expression matrix to just those in this comparison
rma.subset <- rma[paste0('X',samples$array_id[samples$chemical=='IFOSFAMIDE' | (samples$chemical=='Control' & samples$vehicle=='SALINE_100_%')])]

# construct a design matrix modeling treatment vs control for use by limma
design <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle=='SALINE_100_%'],
    levels=c('Control','IFOSFAMIDE')
  )
)
colnames(design) <- c('Intercept','IFOSFAMIDE')

# run limma
fit <- lmFit(rma.subset, design)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')

# filter for unadjusted p-value < 0.05
t.filtered <- subset(t, adj.P.Val < 0.05)
t.filtered <- subset(t, P.Value < 0.05)
# write out the results to file
write.csv(t,'limma_results_ifosfamide_new.csv')
