###Transcriptomic data obtained from the GDC portal is used. The rows are gene names and the columns are samples with cells being TPM values.### 
dat1 <- read.csv("tumor_tpm.csv")
rownames(dat1) <- dat1[,1]
dat1 <- dat1[,-1]
oxp <- readLines("oxp.txt")
ppp <- readLines("ppp.txt")
sel_gn <- c(oxp, ppp)
sum(sel_gn %in% rownames(dat1))
sel_index <- c()
na_index <- c()
for (i in 1:length(sel_gn)) {
  if (sel_gn[i] %in% rownames(dat1)) {sel_index[i] = which(rownames(dat1) == sel_gn[i])
  } else {na_index[i] = paste("index", i, "rm")}
}
sel_index <- sel_index[-which(is.na(sel_index))]
selexprs <- dat1[sel_index,] 
selexp_scaled <- t(apply(selexprs, 1, scale))
colnames(selexp_scaled) <- colnames(selexprs)
boxplot(t(selexp_scaled))

#######

oxp_index <- c()
na_oxp_index <- c()
for (i in 1:length(oxp)) {
  if (oxp[i] %in% rownames(dat1)) {oxp_index[i] = which(rownames(dat1) == oxp[i])
  } else {na_oxp_index[i] = paste("index", i, "rm")}
}
oxp_index <- oxp_index[-which(is.na(oxp_index))]
oxp_exprs <- dat1[oxp_index,] 
oxp_sd <- c()
for (i in 1:nrow(oxp_exprs)) {
  oxp_sd[i] <- sd(as.numeric(oxp_exprs[i,]))
}

ppp_index <- c()
na_ppp_index <- c()
for (i in 1:length(ppp)) {
  if (ppp[i] %in% rownames(dat1)) {ppp_index[i] = which(rownames(dat1) == ppp[i])
  } else {na_ppp_index[i] = paste("index", i, "rm")}
}
#ppp_index <- ppp_index[-which(is.na(ppp_index))]
ppp_exprs <- dat1[ppp_index,] 
ppp_sd <- c()
for (i in 1:nrow(ppp_exprs)) {
  ppp_sd[i] <- sd(as.numeric(ppp_exprs[i,]))
}

oxp_sd <- as.data.frame(oxp_sd, nrow = length(oxp_sd), ncol = 1)
oxp_sd[,2] <- rownames(oxp_exprs)
colnames(oxp_sd) <- c("sd", "gene_name")
plot(oxp_sd[,1])

ppp_sd <- as.data.frame(ppp_sd, nrow = length(ppp_sd), ncol = 1)
ppp_sd[,2] <- rownames(ppp_exprs)
colnames(ppp_sd) <- c("sd", "gene_name")
plot(ppp_sd[,1])

library(ggplot2)
ggplot(oxp_sd, aes(x = "", y = sd)) +
  geom_violin(trim = FALSE, fill = "skyblue", color = "black") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal(base_size = 15) +
  theme(text = element_text(family = "Times New Roman")) +
  labs(x = "OXP_Genes", y = "Standard Deviation")

library(ggplot2)
ggplot(ppp_sd, aes(x = "", y = sd)) +
  geom_violin(trim = FALSE, fill = "skyblue", color = "black") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal(base_size = 15) +
  theme(text = element_text(family = "Times New Roman")) +
  labs(x = "PPP_Genes", y = "Standard Deviation")

quantile(oxp_sd[,1])
quantile(ppp_sd[,1])

##Rules to remove genes is removing genes having less than 50% of the median sd##
oxp_sd_cutoff <- subset(oxp_sd, subset = (oxp_sd[,1] > 14))
oxp_sd_cutoff_gn <- as.vector(oxp_sd_cutoff[,2])
oxp_ppp_cutoff <- subset(ppp_sd, subset = (ppp_sd[,1] > 12))
oxp_ppp_cutoff_gn <- as.vector(oxp_ppp_cutoff[,2])
sel_u_gn2 <- c(oxp_sd_cutoff_gn, oxp_ppp_cutoff_gn)
sum(sel_u_gn2 %in% rownames(dat1))
sel_u_index2 <- c()
na_u_index2 <- c()
for (i in 1:length(sel_u_gn2)) {
  if (sel_u_gn2[i] %in% rownames(dat1)) {sel_u_index2[i] = which(rownames(dat1) == sel_u_gn2[i])
  } else {na_u_index2[i] = paste("index", i, "rm")}
}
u_selexprs2 <- dat1[sel_u_index2,] 
selexp_u_scaled2 <- t(apply(u_selexprs2, 1, scale))
colnames(selexp_u_scaled2) <- colnames(u_selexprs2)
boxplot(t(selexp_u_scaled2))
library(ConsensusClusterPlus)
na_rm_i <- which(is.na(selexp_u_scaled2[,1]))
#selexp_scaled <- selexp_scaled[-na_rm_i,]
results_1_u2 <- ConsensusClusterPlus(t(selexp_u_scaled2), maxK = 10,
                                     reps = 100, pFeature = 0.8, pItem = 1, 
                                     clusterAlg = "hc", distance = "euclidean", innerLinkage = "ward.D2",
                                     seed = 12345678, plot = "pdf")
ksel <- results_1_u2[[5]]
ksel_matrix <- ksel[["consensusMatrix"]]
colnames(ksel_matrix) <- rownames(selexp_u_scaled2)
rownames(ksel_matrix) <- rownames(selexp_u_scaled2)
cluster_assignments_u2 <- ksel$consensusClass
annot_col <- data.frame(Cluster = factor(cluster_assignments_u2))
annot_row <- as.data.frame(c(rep("OXPHOS", times = 81),  rep("PPP", times = 20)))
rownames(annot_row) <- rownames(selexp_u_scaled2)                    
colnames(annot_row) <- c("Pathway")
library(pheatmap)
annot_colour <- list(
  Pathway = c(OXPHOS = "pink", PPP = "red")
)
pheatmap(ksel_matrix,
         annotation_row = annot_row, annotation_col = annot_col, color = colorRampPalette(c("white", "blue"))(100),
         clustering_method = "ward.D2", fontfamily = "Times", border_color = NA, 
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         fontsize_col = 7.2, annotation_colors = annot_colour, fontsize = 12, cellwidth = 5.4,
         show_rownames = F, treeheight_row = 20, treeheight_col = 25, seed = 12345678)

##Extracting the genes of the clusters selected##
gn_cluster <- as.data.frame(cluster_assignments_u2)
gn_cluster[,2] <- rownames(gn_cluster)
colnames(gn_cluster) <- c("cluster_id", "gene_name")
c_3_i <- c(which(gn_cluster[,1] == 3))
c_3_gn <- gn_cluster[c_3_i,2]
c_5_i <- c(which(gn_cluster[,1] == 5))
c_5_gn <- gn_cluster[c_5_i,2]

write(c_3_gn, file = "c_3_gn.txt")
write(c_5_gn, file = "c_5_gn.txt")

sub_c3 <- selexp_u_scaled2[c_3_gn,]
sub_c5 <- selexp_u_scaled2[c_5_gn,]
c3_median <- c()
for (i in 1:ncol(sub_c3)) {
  c3_median[i] <- median(sub_c3[,i]) 
}
c5_median <- c()
for (i in 1:ncol(sub_c5)) {
  c5_median[i] <- median(sub_c5[,i]) 
}
x <- c3_median
y <- c5_median
group <- ifelse(x > 0 & y > 0, "A",   
                ifelse(x < 0 & y > 0, "B",    
                       ifelse(x < 0 & y < 0, "C",   
                              ifelse(x > 0 & y < 0, "D", "Origin")))) 
group_colors <- c("A" = "red", "B" = "blue", "C" = "green", "D" = "purple", "Origin" = "black")
par(family = "Times New Roman")
plot(x, y, col = group_colors[group], pch = 19, cex = 0.5, xlim = c(-1.2, 2.5), ylim = c(-1, 3),
     xlab = "Median OXPHOS gene expression (z-score)",
     ylab = "Median PPP gene expression (z-score)", cex.lab = 1.4, family="Times")
abline(h = 0, v = 0, lwd = 2, col = "black", lty = 2)
legend("topright", 
       legend = c("Mixed", "PPP Leaning", "Quiescent", "OXPHOS Leaning"), 
       col = group_colors[c("A", "B", "C", "D")], 
       pch = 19, cex = 1.5,
       title = "Metabolic Group")

gA_i0 <- which(x > 0 & y > 0)
gB_i0 <- which(x < 0 & y > 0)
gC_i0 <- which(x < 0 & y < 0)
gD_i0 <- which(x > 0 & y < 0)
met_group <- c()
met_group[gA_i0] <- c("A")
met_group[gB_i0] <- c("B")
met_group[gC_i0] <- c("C")
met_group[gD_i0] <- c("D")
met_group <- as.data.frame(met_group)
rownames(met_group) <- colnames(selexp_u_scaled2)
met_group[,2] <- rownames(met_group)
colnames(met_group)[2] <- c("patient_id")

p_order <- subset.data.frame(met_group, select = 1)
p_order <- sort_by.data.frame(p_order, met_group)
patient_list <- rownames(p_order)
p_order <- as.data.frame(gsub("A", "Mixed", p_order[,1]))
p_order <- as.data.frame(gsub("B", "PPP Leaning", p_order[,1]))
p_order <- as.data.frame(gsub("C", "Quiescent", p_order[,1]))
p_order <- as.data.frame(gsub("D", "OXPHOS Leaning", p_order[,1]))
colnames(p_order) <- c("Metabolic Group")
rownames(p_order) <- patient_list
sel_median_gene <- rbind(sub_c3, sub_c5)
sel_median_gene_ordered <- sel_median_gene[,rownames(p_order)]
negbreak <- seq(-2,0, by = 0.05)
posbreak <- seq(0,2, by = 0.05)
cbreak <- c(negbreak, posbreak[-1])
row_annot <- as.data.frame(c(rep("Cluster 3, OXPHOS", times = 26), rep("Cluster 5, PPP", times = 8)))
rownames(row_annot) <- c(c_3_gn, c_5_gn)
colnames(row_annot) <- c("Pathway")
annot_colour <- list(
  Pathway = c("Cluster 3, OXPHOS" = "#21897e", "Cluster 5, PPP" = "#e9c46a")
)
pheatmap(sel_median_gene_ordered, color = colorRampPalette(c("#d7f0f9","#76bce8", "#1c3f9a", "#000000", "#7a0e0e", "#e04b3a", "#fca283"))(81),
         breaks = cbreak, cluster_rows = F, fontsize_row = 10.5, fontsize_col = 1.5, fontfamily = "Times",
         gaps_row = 26, fontsize = 10, cluster_cols = F, annotation_col = p_order, annotation_row = row_annot, annotation_colors = annot_colour,
         gaps_col = c(max(grep("Mixed", p_order[,1])), 
                      max(grep("PPP Leaning", p_order[,1])),
                      max(grep("Quiescent", p_order[,1]))
         )
)

##survival_analysis##
survival_dat <- read.csv("tumor_vitals.csv")
survival_dat <- survival_dat[,c(12,14,15)]
dead_i <- which(survival_dat$demographic.vital_status == "Dead")
alive_i <- which(survival_dat$demographic.vital_status == "Alive")
survival_dat$demographic.vital_status[dead_i] <- c("1")
survival_dat$demographic.vital_status[alive_i] <- c("0")
p_order[,2] <- rownames(p_order)
colnames(p_order)[2] <- c("patient_id")
library(dplyr)
surv_data <- inner_join(survival_dat, p_order, by = "patient_id")
rownames(surv_data) <- surv_data[,1]
surv_data <- surv_data[,-1]
p_rm_i <- which(surv_data$diagnoses.days_to_last_follow_up == "'--")
survival_data_cleaned <- surv_data[-p_rm_i,]
colnames(survival_data_cleaned) <- c("vitals", "days", "met_group")
survival_data_cleaned[,1] <- as.numeric(survival_data_cleaned[,1])
survival_data_cleaned[,2] <- as.numeric(survival_data_cleaned[,2])
survival_data_cleaned <- as.data.frame(survival_data_cleaned)
library(survival)
library(survminer)
library(dplyr)
sfit <- survfit(Surv(days, vitals)~met_group, data = survival_data_cleaned)
fit <- coxph(Surv(days, vitals)~met_group, data = survival_data_cleaned)
summary(sfit)
summary(fit)
fit_km <- survfit(Surv(days, vitals) ~ met_group, data = survival_data_cleaned)
ggsurvplot(fit_km, data = survival_data_cleaned, pval = TRUE,
           risk.table = TRUE,  legend.title = "Metabolic Group",
           legend.labs = c("Mixed", "OXPHOS Leaning", "PPP Leaning", "Quiescent"),
           ggtheme = theme_classic2(base_family = "Times New Roman"), 
           xlab = "Days",
           ylab = "Survival probability",
           risk.table.y.text.col = TRUE,
           risk.table.y.text = TRUE,
           tables.height = 0.15,
           surv.plot.height = 1.2, 
           risk.table.height = 0.2,
           
           font.x = c(13, "black"), 
           font.tickslab = 9, 
           risk.table.fontsize = 4, 
           break.time.by = 500, 
           palette = c(
             "Mixed" = "#F3A6EC",       # light pink
             "OXPHOS Leaning" = "#F08070", # coral red
             "PPP Leaning" = "#94C9F5", # sky blue
             "Quiescent" = "#21D343"     # vivid green
           ))

##Stacked bar plots####Factor_analysis##Type
ln_type <- read.csv("tumor_status_ln_type.csv")
ln_type <- subset(ln_type, select = c(12,14,15))
met_annot <- inner_join(ln_type, p_order, by = "patient_id")
rownames(met_annot) <- met_annot[,1]
met_annot <- met_annot[,-1]
gA_i <- grep("Mixed", met_annot$`Metabolic Group`)
gB_i <- grep("PPP Leaning", met_annot$`Metabolic Group`)
gC_i <- grep("Quiscent", met_annot$`Metabolic Group`)
gD_i <- grep("OXPHOS Leaning", met_annot$`Metabolic Group`)
ga <- met_annot[gA_i,]
gb <- met_annot[gB_i,]
gc <- met_annot[gC_i,]
gd <- met_annot[gD_i,]
groups <- list(ga, gb, gc, gd)
type <- matrix(NA, nrow = 4, ncol = 5)
for (i in 1:4) {
  x <- groups[[i]]
  np <- length(grep("primary", x$diagnoses.classification_of_tumor))
  nm <- length(grep("metastasis", x$diagnoses.classification_of_tumor))
  nr <- length(grep("recurrence", x$diagnoses.classification_of_tumor))
  nsp <- length(grep("Subsequent Primary", x$diagnoses.classification_of_tumor))
  nsyp <- length(grep("Synchronous primary", x$diagnoses.classification_of_tumor))
  type[i,] <- c(np,nm,nr,nsp,nsyp)
}
colnames(type) <- c("Primary", "Metastasis", "Recurrence", "Subsequent Primary", "Synchronous primary")
rownames(type) <- c("Mixed","PPP Leanig", "Quiscent", "OXPHOS Leaning")
type_prop <- prop.table(type, margin = 2)
par(mar = c(5, 4, 4, 8))  # enlarge right margin (default is c(5,4,4,2))
barplot((type_prop), beside = FALSE, col = terrain.colors(5),
        legend.text = TRUE,
        args.legend = list(x = "topright", cex = 0.8, inset = c(-0.15, 0)),
        names.arg = colnames(type_prop),
        main = "Tumor Types",
        ylab = "Proportion")

##for ln status##
ln <- matrix(NA, nrow = 4, ncol = 4)
for (i in 1:4) {
  x <- groups[[i]]
  n0 <- length(grep("N0", x$diagnoses.ajcc_pathologic_n))
  n1 <- length(grep("N1", x$diagnoses.ajcc_pathologic_n))
  n2 <- length(grep("N2", x$diagnoses.ajcc_pathologic_n))
  n3 <- length(grep("N3", x$diagnoses.ajcc_pathologic_n))
  ln[i,] <- c(n0,n1,n2,n3)
}
colnames(ln) <- c("N0","N1", "N2", "N3")
rownames(ln) <- c("Mixed","PPP Leaning", "Quiscent", "OXPHOS Leaning")
ln_prop <- prop.table(ln, margin = 2)
par(mar = c(5, 4, 4, 8))
barplot(ln_prop, beside = FALSE, col = terrain.colors(5),
        legend.text = TRUE,
        args.legend = list(x = "topright", cex = 0.8, inset = c(-0.15, 0)),
        names.arg = colnames(ln_prop),
        main = "Tumor LN status",
        ylab = "Proportion")

##Factor_analysis##Grade
factor_dat <- read.csv("tumor_status_factor_inputs.csv")
factor_dat <- subset(factor_dat, select = c(12,14:19))
met_annot <- inner_join(factor_dat, p_order, by = "patient_id")
rownames(met_annot) <- met_annot[,1]
met_annot <- met_annot[,-1]
gA_i <- grep("Mixed", met_annot$`Metabolic Group`)
gB_i <- grep("PPP Leaning", met_annot$`Metabolic Group`)
gC_i <- grep("Quiscent", met_annot$`Metabolic Group`)
gD_i <- grep("OXPHOS Leaning", met_annot$`Metabolic Group`)
ga <- met_annot[gA_i,]
gb <- met_annot[gB_i,]
gc <- met_annot[gC_i,]
gd <- met_annot[gD_i,]
groups <- list(ga, gb, gc, gd)
type <- matrix(NA, nrow = 4, ncol = 4)
for (i in 1:4) {
  x <- groups[[i]]
  ng1 <- length(grep("G1", x$diagnoses.tumor_grade))
  ng2 <- length(grep("G2", x$diagnoses.tumor_grade))
  ng3 <- length(grep("G3", x$diagnoses.tumor_grade))
  ng4 <- length(grep("G4", x$diagnoses.tumor_grade))
  type[i,] <- c(ng1,ng2,ng3,ng4)
}
colnames(type) <- c("G1", "G2", "G3", "G4")
rownames(type) <- c("Mixed","PPP Leanig", "Quiscent", "OXPHOS Leaning")
type_prop <- prop.table(type, margin = 2)
par(mar = c(5, 4, 4, 8))  # enlarge right margin (default is c(5,4,4,2))
barplot((type_prop), beside = FALSE, col = terrain.colors(5),
        legend.text = TRUE,
        args.legend = list(x = "topright", cex = 0.8, inset = c(-0.15, 0)),
        names.arg = colnames(type_prop),
        main = "Tumor Grade",
        ylab = "Proportion")

##Stage##
stage <- matrix(NA, nrow = 4, ncol = 7)
for (i in 1:4) {
  x <- groups[[i]]
  ns0 <- length(grep("Stage 0", x$diagnoses.ajcc_pathologic_stage))
  ns1 <- length(grep("Stage I", x$diagnoses.ajcc_pathologic_stage))
  ns2 <- length(grep("Stage II", x$diagnoses.ajcc_pathologic_stage))
  ns3 <- length(grep("Stage III", x$diagnoses.ajcc_pathologic_stage))
  ns4a <- length(grep("Stage IVA", x$diagnoses.ajcc_pathologic_stage))
  ns4b <- length(grep("Stage IVB", x$diagnoses.ajcc_pathologic_stage))
  ns4c <- length(grep("Stage IVC", x$diagnoses.ajcc_pathologic_stage))
  stage[i,] <- c(ns0, ns1, ns2, ns3, ns4a, ns4b, ns4c)
}
colnames(stage) <- c("Stage 0", "Stage I", "Stage II", "Stage III", "Stage IVA", "Stage IVB", "Stage IVC")
rownames(stage) <- c("Mixed","PPP Leanig", "Quiscent", "OXPHOS Leaning")
stage_prop <- prop.table(stage, margin = 2)
par(mar = c(5, 4, 4, 8))  # enlarge right margin (default is c(5,4,4,2))
barplot((stage_prop), beside = FALSE, col = terrain.colors(5),
        legend.text = TRUE,
        args.legend = list(x = "topright", cex = 0.8, inset = c(-0.15, 0)),
        names.arg = colnames(stage_prop),
        main = "Tumor Stage",
        ylab = "Proportion")

##T_Satge##
t_stage <- matrix(NA, nrow = 4, ncol = 4)
for (i in 1:4) {
  x <- groups[[i]]
  nt1 <- length(grep("T1", x$diagnoses.ajcc_pathologic_t))
  nt2 <- length(grep("T2", x$diagnoses.ajcc_pathologic_t))
  nt3 <- length(grep("T3", x$diagnoses.ajcc_pathologic_t))
  nt4 <- length(grep("T4", x$diagnoses.ajcc_pathologic_t))
  t_stage[i,] <- c(nt1, nt2, nt3, nt4)
}
colnames(t_stage) <- c("T1", "T2", "T3", "T4")
rownames(t_stage) <- c("Mixed","PPP Leanig", "Quiscent", "OXPHOS Leaning")
tstage_prop <- prop.table(t_stage, margin = 2)
par(mar = c(5, 4, 4, 8))  # enlarge right margin (default is c(5,4,4,2))
barplot((tstage_prop), beside = FALSE, col = terrain.colors(5),
        legend.text = TRUE,
        args.legend = list(x = "topright", cex = 0.8, inset = c(-0.15, 0)),
        names.arg = colnames(tstage_prop),
        main = "Tumor T Stage",
        ylab = "Proportion")

##Proliferation_Score_correlation##
proliferation <- read.csv("proliferation_score.csv")
sample_map <- read.csv("sample_map.csv")
colnames(proliferation)[1] <- c("case_id")
colnames(sample_map)[6] <- c("case_id")
library(dplyr)
proli_joint <- inner_join(proliferation, sample_map, by = "case_id")
proli_joint <- proli_joint[,c(2,13)]
c3_median <- as.data.frame(c3_median)
c3_median$patient_id <- colnames(sub_c3)
c5_median <- as.data.frame(c5_median)
c5_median$patient_id <- colnames(sub_c5)
library(dplyr)
proli_median <- inner_join(proli_joint, c3_median, by = "patient_id")
proli_median <- inner_join(proli_median, c5_median, by = "patient_id")
proli_median <- na.omit(proli_median)
quantile(proli_median$c3_median)
proli_c3median <- subset(proli_median, subset = (proli_median$c3_median < 1 & proli_median$c3_median > -1))
library(ggplot2)
library(ggExtra)
c3_proli <- ggplot(proli_c3median, aes(x = c3_median, y = Proliferation)) +
  geom_point(inherit.aes = TRUE,color = "grey40", size = 1.5) +
  labs(x = "Cluster 3 OXPHOS Median Expression", y = "Proliferation Score") +
  geom_smooth(method = "lm", color = "blue", fill = "lightblue") +
  theme_minimal(base_size = 18, base_family = "Times New Roman") +
  theme(
    text = element_text(family = "Times New Roman")
  )
ggMarginal(
  c3_proli,
  type = "density",
  margins = "both",
  size = 5,
  colour = "black",
  fill = c("#21897e")  # top, right
)
fit <- lm(Proliferation ~ c3_median, data = proli_c3median)
summary(fit)
cor(proli_c3median$Proliferation, proli_c3median$c3_median)
quantile(proli_median$c5_median)
proli_c5median <- subset(proli_median, subset = (proli_median$c5_median < 1 & proli_median$c5_median > -1))
library(ggplot2)
c5_proli <- ggplot(proli_c5median, aes(x = c5_median, y = Proliferation)) +
  geom_point(inherit.aes = TRUE,color = "grey40", size = 1.5) +
  labs(x = "Cluster 5 PPP Median Expression", y = "Proliferation Score") +
  xlim(-1,1) +
  geom_smooth(method = "lm", color = "blue", fill = "lightblue") +
  theme_minimal(base_size = 18, base_family = "Times New Roman") +
  theme(
    text = element_text(family = "Times New Roman")
  )
ggMarginal(
  c5_proli,
  type = "density",
  margins = "both",
  size = 5,
  colour = "black",
  fill = c("#e9c46a")
)
fit <- lm(Proliferation ~ c5_median, data = proli_c5median)
summary(fit)

##Boxplot of proliferation score##
xc3 <- proli_median$c3_median
yc5 <- proli_median$c5_median
gmix_i <- which(xc3 > 0 & yc5 > 0)
gppp_i <- which(xc3 < 0 & yc5 > 0)
gquis_i <- which(xc3 < 0 & yc5 < 0)
goxp_i <- which(xc3 > 0 & yc5 < 0)
met_group <- c()
met_group[gmix_i] <- c("Mixed")
met_group[goxp_i] <- c("OXPHOS Leaning")
met_group[gppp_i] <- c("PPP Leaning")
met_group[gquis_i]<- c("Quiscent")
proli_median$met_group <- met_group
ggboxplot(proli_median, x = "met_group", y = "Proliferation", fill = "met_group") +
  scale_fill_manual(values = c("Mixed" = "pink", 
                               "OXPHOS Leaning" = "orange", 
                               "PPP Leaning" = "blue",
                               "Quiscent" = "green")) +
  annotate("text", x = 1.1, y = 2,
           label = paste0("Kruskal-Wallis, p = ", signif(kruskal.test(Proliferation ~ met_group, data = proli_median)$p.value, 3)),
           hjust = 0, size = 5, family = "Times New Roman") +
  theme_minimal(base_size = 14, base_family = "Times New Roman") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    text = element_text(family = "Times New Roman"),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold")
  ) +
  labs(x = "", y = "Proliferation", fill = "Metabolic Group")

##Performing GSEA of the met groups against the others##
dat2 <- read.csv("tumor_counts.csv")
rownames(dat2) <- dat2[,1]
dat2 <- dat2[,-1]
met_order_p_index <- data.frame(
  mixed <- rep("Others", time = 472),
  oxp_leaning <- rep("Others", time = 472),
  ppp_leaning <- rep("Others", time = 472),
  quiscent <- rep("Others", time = 472)
)
colnames(met_order_p_index) <- c("Mixed", "OXP_leaning", "PPP_leaning", "Quiscent")
met_order_p_index$Mixed[gA_i0] <- c("Mixed")
met_order_p_index$OXP_leaning[gD_i0] <- c("OXP_leaning")
met_order_p_index$PPP_leaning[gB_i0] <- c("PPP_leaning")
met_order_p_index$Quiscent[gC_i0] <- c("Quiscent")
rownames(met_order_p_index) <- colnames(dat2)
met_order_p_index$Mixed <- as.factor(met_order_p_index$Mixed)
met_order_p_index$Mixed <- relevel(met_order_p_index$Mixed, ref = "Others")
met_order_p_index$OXP_leaning <- as.factor(met_order_p_index$OXP_leaning)
met_order_p_index$OXP_leaning <- relevel(met_order_p_index$OXP_leaning, ref = "Others")
met_order_p_index$PPP_leaning <- as.factor(met_order_p_index$PPP_leaning)
met_order_p_index$PPP_leaning <- relevel(met_order_p_index$PPP_leaning, ref = "Others")
met_order_p_index$Quiscent <- as.factor(met_order_p_index$Quiscent)
met_order_p_index$Quiscent <- relevel(met_order_p_index$Quiscent, ref = "Others")

##Mixed##
library(DESeq2)
dds_mixed <- DESeqDataSetFromMatrix(
  countData = dat2,
  colData = met_order_p_index,
  design = ~ Mixed
)
keep_mixed <- rowSums(counts(dds_mixed)) > 4720
dds_mixed <- dds_mixed[keep_mixed,]
ddsM <- DESeq(dds_mixed)
resM <- results(ddsM)
resM_df <- as.data.frame(resM)
resM_filtered <- subset(resM_df, subset = (resM_df$padj < 0.05))
resM_filtered <- sort_by(resM_filtered, resM_filtered$log2FoldChange, decreasing = T)
resM_filtered$gene_symbol <- rownames(resM_filtered)

############################################################################################################
library(org.Hs.eg.db)
gn_mixed <- rownames(resM_filtered)
EG_IDs <- mget(gn_mixed, revmap(org.Hs.egSYMBOL),ifnotfound=NA)
EG_IDs <- as.data.frame(unlist(EG_IDs))
colnames(EG_IDs) <- c("Entrez_ID")
EG_IDs$gene_symbol <- rownames(EG_IDs)
library(dplyr)
gsea_input <- inner_join(resM_filtered, EG_IDs, by = "gene_symbol")
gsea_input <- gsea_input[,c(2,8)]
un_entrez_gn <- unique(gsea_input[,2])
un_entrez_gn_i <- match( un_entrez_gn, gsea_input[,2])
gsea_input <- gsea_input[un_entrez_gn_i,]
gsea_input <- na.omit(gsea_input)
gsea_gn_fc <- c(gsea_input[,1])
names(gsea_gn_fc) <- gsea_input[,2]
library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
gsea <- gseKEGG(geneList = gsea_gn_fc, keyType = "ncbi-geneid", by = "fgsea")
#saveRDS(gsea, file = "gsea_result.rds")
gsea_res <- as.data.frame(gsea)

#write.csv(gsea_res, file = "gsea_mixed.csv")
gseaplot2(gsea, geneSetID = c("hsa00480",
"hsa00190",
"hsa05208",
"hsa00040",
"hsa01230",
"hsa01200",
"hsa00030",
"hsa04514"), base_size = 15,  rel_heights = c(1, 0.4, 0.5))

library(patchwork)
gsea_plot <- gseaplot2(
  gsea,
  geneSetID = c(
    "hsa00480", "hsa00190", "hsa05208", "hsa00040", "hsa01230",
    "hsa01200", "hsa00030", "hsa04514"
  ),
  base_size = 15,
  rel_heights = c(1, 0.3, 0.3)
)

# Apply Times New Roman font to each plot
gsea_plot_fonts <- lapply(gsea_plot, function(p) {
  p + theme(text = element_text(family = "Times New Roman"),
            legend.text = element_text(size = 20)
            )
})

# Combine all plots into one figure using patchwork
combined_plot <- wrap_plots(gsea_plot_fonts, , heights = c(1, 0.3, 0.3))

# Display combined plot
print(combined_plot)

#################################################################################################
library(ggplot2)
library(patchwork)
library(cowplot)
library(gtable)

# ===== Selected pathways =====
selected_pathways <- c(
  "hsa00480", "hsa00190", "hsa05208", "hsa00040",
  "hsa01230", "hsa01200", "hsa00030", "hsa04514"
)

# ===== Create GSEA plots =====
gsea_plot <- gseaplot2(
  gsea,
  geneSetID = selected_pathways,
  base_size = 18,
  rel_heights = c(1, 0.3, 0.3)
)

# ===== Apply Times New Roman font and borders =====
gsea_plot_fonts <- lapply(gsea_plot, function(p) {
  p + theme(
    text = element_text(family = "Times New Roman"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    panel.border = element_rect(color = "black", fill = NA) # border box
  )
})

# ===== Extract legend from first plot =====
get_legend <- function(my_plot) {
  tmp <- ggplotGrob(my_plot + theme(legend.position = "right"))
  gtable::gtable_filter(tmp, "guide-box")
}
legend_only_plot <- ggdraw(get_legend(gsea_plot_fonts[[1]]))

# ===== Remove legends from all subplots =====
plots_no_legends <- lapply(gsea_plot_fonts, function(p) {
  p + theme(legend.position = "none")
})

# ===== Combine all plots vertically =====
combined_plot <- wrap_plots(plots_no_legends, heights = c(1, 0.3, 0.3))

# ===== Show results =====
print(combined_plot)    # Plots without legends
print(legend_only_plot)
###########################################################################################

##OXP_Leaning##
library(DESeq2)
dds_oxp <- DESeqDataSetFromMatrix(
  countData = dat2,
  colData = met_order_p_index,
  design = ~ OXP_leaning
)
keep_oxp <- rowSums(counts(dds_oxp)) > 4720
dds_oxp <- dds_oxp[keep_oxp,]
ddsO <- DESeq(dds_oxp)
resO <- results(ddsO)
resO_df <- as.data.frame(resO)
resO_filtered <- subset(resO_df, subset = (resO_df$padj < 0.05))
resO_filtered <- sort_by(resO_filtered, resO_filtered$log2FoldChange, decreasing = T)
resO_filtered$gene_symbol <- rownames(resO_filtered)

library(org.Hs.eg.db)
gn_oxp <- rownames(resO_filtered)
EG_IDs <- mget(gn_oxp, revmap(org.Hs.egSYMBOL),ifnotfound=NA)
EG_IDs <- as.data.frame(unlist(EG_IDs))
colnames(EG_IDs) <- c("Entrez_ID")
EG_IDs$gene_symbol <- rownames(EG_IDs)
library(dplyr)
gsea_input <- inner_join(resO_filtered, EG_IDs, by = "gene_symbol")
gsea_input <- gsea_input[,c(2,8)]
un_entrez_gn <- unique(gsea_input[,2])
un_entrez_gn_i <- match( un_entrez_gn, gsea_input[,2])
gsea_input <- gsea_input[un_entrez_gn_i,]
gsea_input <- na.omit(gsea_input)
gsea_gn_fc <- c(gsea_input[,1])
names(gsea_gn_fc) <- gsea_input[,2]
library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
gsea_oxp <- gseKEGG(geneList = gsea_gn_fc, keyType = "ncbi-geneid", by = "fgsea")
gsea_res_oxp <- as.data.frame(gsea_oxp)
#write.csv(gsea_res_oxp, file = "gsea_oxp.csv")
gseaplot2(gsea_oxp, geneSetID = c("hsa00190",
"hsa05208",
"hsa00240",
"hsa04066",
"hsa04010",
"hsa04630",
"hsa00590",
"hsa05200",
"hsa04151",
"hsa00040"), base_size = 15,  rel_heights = c(1, 0.4, 0.5))

library(patchwork)
gsea_plot <- gseaplot2(
  gsea_oxp,
  geneSetID = c("hsa00190",
                "hsa05208",
                "hsa00240",
                "hsa04066",
                "hsa04010",
                "hsa04630",
                "hsa00590",
                "hsa05200",
                "hsa04151",
                "hsa00040"),
  base_size = 18,
  rel_heights = c(1, 0.4, 0.5)
)

# Apply Times New Roman font to each plot
gsea_plot_fonts <- lapply(gsea_plot, function(p) {
  p + theme(text = element_text(family = "Times New Roman"))
})

# Combine all plots into one figure using patchwork
combined_plot <- wrap_plots(gsea_plot_fonts, ncol = 1, heights = c(1, 0.4, 0.5))

# Display combined plot
print(combined_plot)

#################################################################################################
library(ggplot2)
library(patchwork)
library(cowplot)
library(gtable)

# ===== Selected pathways =====
selected_pathways <- c("hsa00190",
                       "hsa05208",
                       "hsa00240",
                       "hsa04066",
                       "hsa04010",
                       "hsa04630",
                       "hsa00590",
                       "hsa05200",
                       "hsa04151",
                       "hsa00040")

# ===== Create GSEA plots =====
gsea_plot <- gseaplot2(
  gsea_oxp,
  geneSetID = selected_pathways,
  base_size = 18,
  rel_heights = c(1, 0.3, 0.3)
)

# ===== Apply Times New Roman font and borders =====
gsea_plot_fonts <- lapply(gsea_plot, function(p) {
  p + theme(
    text = element_text(family = "Times New Roman"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    panel.border = element_rect(color = "black", fill = NA) # border box
  )
})

# ===== Extract legend from first plot =====
get_legend <- function(my_plot) {
  tmp <- ggplotGrob(my_plot + theme(legend.position = "right"))
  gtable::gtable_filter(tmp, "guide-box")
}
legend_only_plot <- ggdraw(get_legend(gsea_plot_fonts[[1]]))

# ===== Remove legends from all subplots =====
plots_no_legends <- lapply(gsea_plot_fonts, function(p) {
  p + theme(legend.position = "none")
})

# ===== Combine all plots vertically =====
combined_plot <- wrap_plots(plots_no_legends, heights = c(1, 0.3, 0.3))

# ===== Show results =====
print(combined_plot)    # Plots without legends
print(legend_only_plot)
###########################################################################################

##PPP_leaning##
library(DESeq2)
dds_ppp <- DESeqDataSetFromMatrix(
  countData = dat2,
  colData = met_order_p_index,
  design = ~ PPP_leaning
)
keep_ppp <- rowSums(counts(dds_ppp)) > 4720
dds_ppp <- dds_ppp[keep_ppp,]
ddsP <- DESeq(dds_ppp)
resP <- results(ddsP)
resP_df <- as.data.frame(resP)
resP_filtered <- subset(resP_df, subset = (resP_df$padj < 0.05))
resP_filtered <- sort_by(resP_filtered, resP_filtered$log2FoldChange, decreasing = T)
resP_filtered$gene_symbol <- rownames(resP_filtered)

library(org.Hs.eg.db)
gn_ppp <- rownames(resP_filtered)
EG_IDs <- mget(gn_ppp, revmap(org.Hs.egSYMBOL),ifnotfound=NA)
EG_IDs <- as.data.frame(unlist(EG_IDs))
colnames(EG_IDs) <- c("Entrez_ID")
EG_IDs$gene_symbol <- rownames(EG_IDs)
library(dplyr)
gsea_input <- inner_join(resP_filtered, EG_IDs, by = "gene_symbol")
gsea_input <- gsea_input[,c(2,8)]
un_entrez_gn <- unique(gsea_input[,2])
un_entrez_gn_i <- match( un_entrez_gn, gsea_input[,2])
gsea_input <- gsea_input[un_entrez_gn_i,]
gsea_input <- na.omit(gsea_input)
gsea_gn_fc <- c(gsea_input[,1])
names(gsea_gn_fc) <- gsea_input[,2]
library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
gsea_ppp <- gseKEGG(geneList = gsea_gn_fc, keyType = "ncbi-geneid", by = "fgsea")
gsea_res_ppp <- as.data.frame(gsea_ppp)
#write.csv(gsea_res_ppp, file = "gsea_ppp.csv")
#dotplot(gsea_ppp, showCategory=20, split=".sign",  font.size = 10) + facet_grid(.~.sign)
gseaplot2(gsea_ppp, geneSetID = c("hsa00480",
"hsa00040",
"hsa00010",
"hsa01200",
"hsa00590",
"hsa00030",
"hsa05200",
"hsa00190",
"hsa03320",
"hsa04630"),  base_size = 15,  rel_heights = c(1, 0.4, 0.5))

library(patchwork)
gsea_plot <- gseaplot2(
  gsea_ppp,
  geneSetID = c("hsa00480",
                "hsa00040",
                "hsa00010",
                "hsa01200",
                "hsa00590",
                "hsa00030",
                "hsa05200",
                "hsa00190",
                "hsa03320",
                "hsa04630"),
  base_size = 18,
  rel_heights = c(1, 0.4, 0.5)
)

# Apply Times New Roman font to each plot
gsea_plot_fonts <- lapply(gsea_plot, function(p) {
  p + theme(text = element_text(family = "Times New Roman"))
})

# Combine all plots into one figure using patchwork
combined_plot <- wrap_plots(gsea_plot_fonts, ncol = 1, heights = c(1, 0.4, 0.5))

# Display combined plot
print(combined_plot)

#################################################################################################
library(ggplot2)
library(patchwork)
library(cowplot)
library(gtable)

# ===== Selected pathways =====
selected_pathways <- c("hsa00480",
                       "hsa00040",
                       "hsa00010",
                       "hsa01200",
                       "hsa00590",
                       "hsa00030",
                       "hsa05200",
                       "hsa00190",
                       "hsa03320",
                       "hsa04630")

# ===== Create GSEA plots =====
gsea_plot <- gseaplot2(
  gsea_ppp,
  geneSetID = selected_pathways,
  base_size = 18,
  rel_heights = c(1, 0.3, 0.3)
)

# ===== Apply Times New Roman font and borders =====
gsea_plot_fonts <- lapply(gsea_plot, function(p) {
  p + theme(
    text = element_text(family = "Times New Roman"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    panel.border = element_rect(color = "black", fill = NA) # border box
  )
})

# ===== Extract legend from first plot =====
get_legend <- function(my_plot) {
  tmp <- ggplotGrob(my_plot + theme(legend.position = "right"))
  gtable::gtable_filter(tmp, "guide-box")
}
legend_only_plot <- ggdraw(get_legend(gsea_plot_fonts[[1]]))

# ===== Remove legends from all subplots =====
plots_no_legends <- lapply(gsea_plot_fonts, function(p) {
  p + theme(legend.position = "none")
})

# ===== Combine all plots vertically =====
combined_plot <- wrap_plots(plots_no_legends, heights = c(1, 0.3, 0.3))

# ===== Show results =====
print(combined_plot)    # Plots without legends
print(legend_only_plot)
###########################################################################################


##Quiscent##
library(DESeq2)
dds_q <- DESeqDataSetFromMatrix(
  countData = dat2,
  colData = met_order_p_index,
  design = ~ Quiscent
)
keep_q <- rowSums(counts(dds_q)) > 4720
dds_q <- dds_q[keep_q,]
ddsQ <- DESeq(dds_q)
resQ <- results(ddsQ)
resQ_df <- as.data.frame(resQ)
resQ_filtered <- subset(resQ_df, subset = (resQ_df$padj < 0.05))
resQ_filtered <- sort_by(resQ_filtered, resQ_filtered$log2FoldChange, decreasing = T)
resQ_filtered$gene_symbol <- rownames(resQ_filtered)

library(org.Hs.eg.db)
gn_q <- rownames(resQ_filtered)
EG_IDs <- mget(gn_q, revmap(org.Hs.egSYMBOL),ifnotfound=NA)
EG_IDs <- as.data.frame(unlist(EG_IDs))
colnames(EG_IDs) <- c("Entrez_ID")
EG_IDs$gene_symbol <- rownames(EG_IDs)
library(dplyr)
gsea_input <- inner_join(resQ_filtered, EG_IDs, by = "gene_symbol")
gsea_input <- gsea_input[,c(2,8)]
un_entrez_gn <- unique(gsea_input[,2])
un_entrez_gn_i <- match( un_entrez_gn, gsea_input[,2])
gsea_input <- gsea_input[un_entrez_gn_i,]
gsea_input <- na.omit(gsea_input)
gsea_gn_fc <- c(gsea_input[,1])
names(gsea_gn_fc) <- gsea_input[,2]
library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
gsea_q <- gseKEGG(geneList = gsea_gn_fc, keyType = "ncbi-geneid", by = "fgsea")
gsea_res_q <- as.data.frame(gsea_q)
#write.csv(gsea_res_q, file = "gsea_q.csv")
#dotplot(gsea_q, showCategory=20, split=".sign",  font.size = 10) + facet_grid(.~.sign)
gseaplot2(gsea_q, geneSetID = c("hsa04630",
"hsa04514",
"hsa04022",
"hsa05200",
"hsa00790",
"hsa00030",
"hsa01200",
"hsa00040",
"hsa00480",
"hsa00190",
"hsa05208"), base_size = 15,  rel_heights = c(1, 0.4, 0.5))

library(patchwork)
gsea_plot <- gseaplot2(
  gsea_q,
  geneSetID = c("hsa04630",
                "hsa04514",
                "hsa04022",
                "hsa05200",
                "hsa00790",
                "hsa00030",
                "hsa01200",
                "hsa00040",
                "hsa00480",
                "hsa00190",
                "hsa05208"),
  base_size = 18,
  rel_heights = c(1, 0.4, 0.5)
)

# Apply Times New Roman font to each plot
gsea_plot_fonts <- lapply(gsea_plot, function(p) {
  p + theme(text = element_text(family = "Times New Roman"))
})

# Combine all plots into one figure using patchwork
combined_plot <- wrap_plots(gsea_plot_fonts, ncol = 1, heights = c(1, 0.4, 0.5))

# Display combined plot
print(combined_plot)

#################################################################################################
library(ggplot2)
library(patchwork)
library(cowplot)
library(gtable)

# ===== Selected pathways =====
selected_pathways <- c("hsa04630",
                       "hsa04514",
                       "hsa04022",
                       "hsa05200",
                       "hsa00790",
                       "hsa00030",
                       "hsa01200",
                       "hsa00040",
                       "hsa00480",
                       "hsa00190",
                       "hsa05208")

# ===== Create GSEA plots =====
gsea_plot <- gseaplot2(
  gsea_q,
  geneSetID = selected_pathways,
  base_size = 18,
  rel_heights = c(1, 0.3, 0.3)
)

# ===== Apply Times New Roman font and borders =====
gsea_plot_fonts <- lapply(gsea_plot, function(p) {
  p + theme(
    text = element_text(family = "Times New Roman"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    panel.border = element_rect(color = "black", fill = NA) # border box
  )
})

# ===== Extract legend from first plot =====
get_legend <- function(my_plot) {
  tmp <- ggplotGrob(my_plot + theme(legend.position = "right"))
  gtable::gtable_filter(tmp, "guide-box")
}
legend_only_plot <- ggdraw(get_legend(gsea_plot_fonts[[1]]))

# ===== Remove legends from all subplots =====
plots_no_legends <- lapply(gsea_plot_fonts, function(p) {
  p + theme(legend.position = "none")
})

# ===== Combine all plots vertically =====
combined_plot <- wrap_plots(plots_no_legends, heights = c(1, 0.3, 0.3))

# ===== Show results =====
print(combined_plot)    # Plots without legends
print(legend_only_plot)
###########################################################################################
