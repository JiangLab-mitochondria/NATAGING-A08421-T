library(edgeR)
set.seed(1234)

data <- read.csv("TE-3-samples-HighMut-vs-LowMut-gene-counts-featurecounts-results.txt", sep = "\t")
group <- factor(c("Mut", "Mut", "Mut", 
                  "Ctr", "Ctr", "Ctr"))
group <- relevel(group, ref = "Ctr")
y <- DGEList(counts = data, group = group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
df <- as.data.frame(topTags(qlf, n=200000))
df$gene.name <- rownames(df)
df <- df %>%
  mutate(gene.type = case_when(logFC >= log2(1.5) & PValue <= 0.05 ~ "up",
                               logFC <= -log2(1.5) & PValue <= 0.05 ~ "down",
                               TRUE ~ "ns")) 