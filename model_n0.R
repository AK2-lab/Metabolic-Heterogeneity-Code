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
  labs(x = "OXP_Genes", y = "Standard Deviation")

library(ggplot2)
ggplot(ppp_sd, aes(x = "", y = sd)) +
  geom_violin(trim = FALSE, fill = "skyblue", color = "black") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal(base_size = 15) +
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
#selexp_u_scaled2 <- t(apply(u_selexprs2, 1, scale))
#colnames(selexp_u_scaled2) <- colnames(u_selexprs2)
#boxplot(t(selexp_u_scaled2))

train_means <- rowMeans(u_selexprs2, na.rm = TRUE)
train_sds   <- apply(u_selexprs2, 1, sd, na.rm = TRUE)

selexp_u_scaled2 <- sweep(u_selexprs2, 1, train_means, "-")
selexp_u_scaled2 <- sweep(selexp_u_scaled2, 1, train_sds, "/")
boxplot(t(selexp_u_scaled2))

##survival_data##
##survival_analysis##
survival_dat <- read.csv("tumor_vitals.csv")
survival_dat <- survival_dat[,c(12,14,15)]
dead_i <- which(survival_dat$demographic.vital_status == "Dead")
alive_i <- which(survival_dat$demographic.vital_status == "Alive")
survival_dat$demographic.vital_status[dead_i] <- c("1")
survival_dat$demographic.vital_status[alive_i] <- c("0")
p_rm_i <- which(survival_dat$diagnoses.days_to_last_follow_up == "'--")
survival_data_cleaned <- survival_dat[-p_rm_i,]
colnames(survival_data_cleaned) <- c("patient_id", "status", "time")
rownames(survival_data_cleaned) <- survival_data_cleaned[,1]
survival_data_cleaned <- survival_data_cleaned[,-1]
survival_data_cleaned[,1] <- as.numeric(survival_data_cleaned[,1])
survival_data_cleaned[,2] <- as.numeric(survival_data_cleaned[,2])
survival_data_cleaned <- as.matrix(survival_data_cleaned)

##Fitting Lasso and GLMNET## ##genes having signif in survival are selected for risk model creation with HR as the coeffecients##
x <- t(selexp_u_scaled2)
y <- survival_data_cleaned
p_list <- rownames(y)
x <- x[p_list,]
library(glmnet)
set.seed(12345678)
cvfit <- cv.glmnet(x, y, family = "cox", type.measure = "C")
plot(cvfit, family = "Times", cex.lab = 1.3)
lam <- cvfit$lambda.min
fit <- glmnet(x, y, family = "cox")
plot(fit, family = "Times", cex.lab = 1.3)
selected_coef <- coef(fit, s = lam)
gene_index <- which(selected_coef != 0)
selected_gene <- rownames(selected_coef)[gene_index]
lassoexprs <- x[,selected_gene]
mucoxd <- cbind(lassoexprs, y)
library(survival)
mucoxd1 <- as.data.frame(mucoxd)
cox_formula <- as.formula(paste("Surv(time, status) ~", paste(selected_gene, collapse = "+")))
mucox_all <- coxph(cox_formula, data = mucoxd1)
summary(mucox_all)
mu_var_sel_gen <- c("ATP5F1D", "ATP5F1E", "ATP5PB", "ATP5PF", "NDUFA2", "NDUFA4L2", "NDUFB5", "TCIRG1", 
                    "UQCRH", "ALDOC", "G6PD", "GPI", "PFKM", "PRPS2", "TALDO1")
x1 <- x[,mu_var_sel_gen]
library(glmnet)
set.seed(12345678)
cvfit <- cv.glmnet(x1, y, family = "cox", type.measure = "C")
plot(cvfit, family = "Times", cex.lab = 1.3)
lam <- cvfit$lambda.min
fit <- glmnet(x1, y, family = "cox")
plot(fit,  family = "Times", cex.lab = 1.3)
selected_coef <- coef(fit, s = lam)
gene_index <- which(selected_coef != 0)
selected_gene <- rownames(selected_coef)[gene_index]
x1 <- x1[,selected_gene]
##Risk Score Calculation##
rs <- (-0.11376153*(x[,"ATP5F1D"])) + 
  (0.15883126*(x[,"ATP5PB"])) + 
  (0.30887920*(x[,"ATP5PF"])) +
  (-0.15810955*(x[,"NDUFA2"])) +
  (-0.06622706*(x[,"NDUFA4L2"])) +
  (0.20485381*(x[,"NDUFB5"])) +
  (0.11144366*(x[,"TCIRG1"])) +
  (-0.29045146*(x[,"UQCRH"])) +
  (0.14899363*(x[,"ALDOC"])) +
  (0.13617248*(x[,"G6PD"])) +
  (0.12304353*(x[,"GPI"])) +
  (0.11039748*(x[,"PFKM"])) +
  (-0.17829467*(x[,"PRPS2"]))
risk_score <- cbind(x1, rs)
range(rs)
hrisk_i <- which(rs > 0)
lrisk_i <- which(rs < 0)

##Survivability Prediction##
cat <- "NA"[1:373]
cat[hrisk_i] <- "HR"
cat[lrisk_i] <- "LR"
risk0 <- cbind(y,cat)
risk0[,3] <- as.factor(risk0[,3])
risk0 <- as.data.frame(risk0)
risk0[,1] <- as.numeric(risk0[,1])
risk0[,2] <- as.numeric(risk0[,2])
risk0[,3] <- as.factor(cat)
sfit <- survfit(Surv(time, status)~as.factor(cat), data=risk0)
summary(sfit)
plot(sfit)
library(survminer)
ggsurvplot(sfit, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("HR", "LR"), legend.title="Risk-Score",  
           palette=c("dodgerblue2", "orchid2"), 
           risk.table.height=.15, xlab = "Days",
           break.time.by = 500)

##Risk Score Distribution##
hrisk <- rs[hrisk_i]
lrisk <- rs[lrisk_i]
hrisk <- sort(hrisk)
lrisk <- sort(lrisk)
risk <- c(lrisk, hrisk)
risk <- as.data.frame(risk)
category <- c(rep("orchid2", times = length(lrisk)), rep("dodgerblue", times = length(hrisk)))
risk[,2] <- category
risk[,3] <- c(1:373)
colnames(risk) <- c("risk","colour","index")
par(family = "Times New Roman")
plot(risk[,1], col = risk$colour, ylab = "Risk Score", pch = 16, xlab = "", cex.lab = 1.4, family="Times")
abline(v = length(lrisk), lty = 2, col = "red", lwd = 2)
legend("bottomright", legend = c("LR", "HR"),              
       col = c("orchid2", "dodgerblue"),
       pch = 16, cex = 1.5,                              
       title = "Risk Group")

##################################################################Reviewer Comment######################
##Survival Info##
cancer <- c("HNSC")
i <- 1
path <- paste0("/media/cuser17/DATA/Pan_Cancer/", cancer[i], "/clinical.cart.2025-08-20")
setwd(path)
survival <- read.csv("clinical.tsv", sep = "\t")
survival <- survival[,c(10,164)]
na_index <- which(survival$diagnoses.year_of_diagnosis == "'--")
survival[na_index,2] <- NA
survival <- survival[complete.cases(survival),]
unique_case <- unique(survival$cases.submitter_id)
unique_case_index <- match(unique_case, survival$cases.submitter_id)
survival <- survival[unique_case_index,]
colnames(survival)[1] <- c("Case.ID")
setwd("/home/cuser17/Documents/TCGA_HNSC_476/MS_work/consensus_clustering/Genelist_090625")
sample_map <- read.csv("sample_map.csv")  
sample_map <- sample_map[,c(6, 12)]  
library(dplyr)
combined <- inner_join(survival, sample_map, by = "Case.ID")
combined <- combined[,-1] 
combined$diagnoses.year_of_diagnosis <- as.numeric(combined$diagnoses.year_of_diagnosis)

combined$year_bin <- cut(
  combined$diagnoses.year_of_diagnosis,
  breaks = 4,
  include.lowest = TRUE
)

year_list <- split(combined, combined$year_bin)
risk_score_list <- c()
for (i in 1:length(year_list)) {
  patient_id <- year_list[[i]][,2]
  patient_idx <- match(patient_id, rownames(risk_score))
  patient_idx <- patient_idx[complete.cases(patient_idx)]
  risk_score_list[[i]] <- risk_score[patient_idx,]
}
##############################################################################################################
j <- 4
total_patient <- 1:nrow(risk_score_list[[j]])

##training set##
random_sample <- sample(total_patient, 0.7*length(total_patient), replace = FALSE, set.seed(123457)) 
test_order <- (random_sample)
test_set <- risk_score_list[[j]][test_order,]
test_patient <- rownames(test_set)
test_data <- cbind(test_set, risk[test_patient,2])
lr_test_i <- which(test_data[,14] < 0)
hr_test_i <- which(test_data[,14] > 0)
lr_test <- as.data.frame(sort((test_data[,14])[lr_test_i], decreasing = T))
hr_test <- as.data.frame(sort((test_data[,14])[hr_test_i], decreasing = F))
colnames(lr_test) <- c("risk_score")
colnames(hr_test) <- c("risk_score")
order <- rbind(lr_test, hr_test)
order_index <- rownames(order)
test_data <- test_data[order_index,]
colnames(test_data)[c(14,15)] <- c("risk_score", "colour")
plot((test_data[,14]), col = test_data[,15], ylab = "Risk Score", pch = 16, xlab = "",
     main = "Risk score distribution in Internal Training set", cex.lab = 1.4, family="Times")
abline(v = length(lr_test_i), lty = 2, col = "red", lwd = 2)
legend("bottomright", legend = c("LR", "HR"),              
       col = c("orchid2", "dodgerblue"),
       pch = 16, cex = 1.4,                   
       title = "Risk Group")
library(pheatmap)
heat_test_data <- test_set[,-c(14)]
heat_test_data <- heat_test_data[rownames(test_data),]
negbreak <- seq(-2,0, by =0.05)
posbreak <- seq(0,2, by =0.05)
cbreak <- c(negbreak, posbreak[-1])
risk_annot <- c(rep("LR", times = length(lr_test_i)), rep("HR", times = length(hr_test_i)))
risk_annot <- as.data.frame(risk_annot)
rownames(risk_annot) <- rownames(heat_test_data)
colnames(risk_annot) <- c("Risk Group")
pheatmap(t(heat_test_data), color = colorRampPalette(c("#d7f0f9","#76bce8", "#1c3f9a", "#000000", "#7a0e0e", "#e04b3a", "#fca283"))(81),
         fontsize_col = 1, fontsize = 11, cluster_rows = F, fontfamily = "Times", breaks = cbreak, fontsize_row = 11, cluster_cols = F,
         annotation_col = risk_annot, gaps_col = length(lr_test_i))


##validation set##
vald_index <- which(!(total_patient) %in% random_sample)
val_set <- risk_score_list[[j]][vald_index,]
val_patient <- rownames(val_set)
val_data <- cbind(val_set, risk[val_patient,2])
lr_test_i_v <- which(val_data[,14] < 0)
hr_test_i_v <- which(val_data[,14] > 0)
lr_test_v <- as.data.frame(sort((val_data[,14])[lr_test_i_v], decreasing = T))
hr_test_v <- as.data.frame(sort((val_data[,14])[hr_test_i_v], decreasing = F))
colnames(lr_test_v) <- c("risk_score")
colnames(hr_test_v) <- c("risk_score")
order_v <- rbind(lr_test_v, hr_test_v)
order_index_v <- rownames(order_v)
val_data <- val_data[order_index_v,]
colnames(val_data)[c(14,15)] <- c("risk_score", "colour")
plot((val_data[,14]), col = val_data[,15], ylab = "Risk Score", pch = 16, cex.lab = 1.3,
     main = "Risk score distribution in Internal Validation set")
abline(v = length(lr_test_i_v), lty = 2, col = "red", lwd = 2)
legend("bottomright", legend = c("LR", "HR"),              
       col = c("orchid2", "dodgerblue"),
       pch = 16,cex = 1.4,                             
       title = "Risk Group")
library(pheatmap)
heat_val_data <- val_set[,-c(14)]
heat_val_data <- heat_val_data[rownames(val_data),]
negbreak <- seq(-2,0, by =0.05)
posbreak <- seq(0,2, by =0.05)
cbreak <- c(negbreak, posbreak[-1])
risk_annot <- c(rep("LR", times = length(lr_test_i_v)), rep("HR", times = length(hr_test_i_v)))
risk_annot <- as.data.frame(risk_annot)
rownames(risk_annot) <- rownames(heat_val_data)
colnames(risk_annot) <- c("Risk Group")
pheatmap(t(heat_val_data), color = colorRampPalette(c("#d7f0f9","#76bce8", "#1c3f9a", "#000000", "#7a0e0e", "#e04b3a", "#fca283"))(81),
         fontsize_col = 3, fontfamily = "Times", fontsize = 11, cluster_rows = F, breaks = cbreak, fontsize_row = 11, cluster_cols = F,
         annotation_col = risk_annot, border_color = NA, gaps_col = length(lr_test_i_v))

##tdROC curve in training##
train <- test_set[,c(14)][rownames(test_data)]
train_os <- y[rownames(test_data),]
training_data <- cbind(train, train_os)
training_data <- as.data.frame(training_data)
library(survivalROC)
roc <- survivalROC(Stime = c(training_data$time), status = c(training_data$status), marker = c(training_data$train),
                   predict.time = 365, method = "KM")
plot(roc$FP, roc$TP, type = "l", col = "blue", xlab = "False Positive Rate",
     ylab = "True Positive Rate", main = "Time-dependent ROC for survival", cex.lab = 1.4 )
abline(a = 0, b = 1, lty = 2, col = "red")
text(x = 0.85, y = 0.12, cex = 1.4, label = paste0("AUC = ", sprintf("%.4f",roc$AUC)))

##tdROC curve in validation##
validation <- val_set[,c(14)][rownames(val_data)]
val_os <- y[rownames(val_data),]
validation_data <- cbind(validation, val_os)
validation_data <- as.data.frame(validation_data)
roc_v <- survivalROC(Stime = c(validation_data$time), status = c(validation_data$status), marker = c(validation_data$validation),
                     predict.time = 365, method = "KM")
plot(roc_v$FP, roc_v$TP, type = "l", col = "blue", xlab = "False Positive Rate",
     ylab = "True Positive Rate", main = "Time-dependent ROC for survival", cex.lab = 1.4)
abline(a = 0, b = 1, lty = 2, col = "red")
text(x = 0.85, y = 0.12, cex = 1.4, label = paste0("AUC = ", sprintf("%.4f",roc_v$AUC)))

##External Validation##
x11 <- read.csv("cptac3_unique_file_id_tumor.csv")
rownames(x11) <- x11[,1]
x11 <- x11[,-1]
x11 <- x11[selected_gene,]
x11_scaled <- sweep(x11, 1, train_means, "-")
x11_scaled <- sweep(x11_scaled, 1, train_sds, "/")
x11_scaled <- t(x11_scaled)
rs1 <- (-0.11376153*(x11_scaled[,"ATP5F1D"])) + 
  (0.15883126*(x11_scaled[,"ATP5PB"])) + 
  (0.30887920*(x11_scaled[,"ATP5PF"])) +
  (-0.15810955*(x11_scaled[,"NDUFA2"])) +
  (-0.06622706*(x11_scaled[,"NDUFA4L2"])) +
  (0.20485381*(x11_scaled[,"NDUFB5"])) +
  (0.11144366*(x11_scaled[,"TCIRG1"])) +
  (-0.29045146*(x11_scaled[,"UQCRH"])) +
  (0.14899363*(x11_scaled[,"ALDOC"])) +
  (0.13617248*(x11_scaled[,"G6PD"])) +
  (0.12304353*(x11_scaled[,"GPI"])) +
  (0.11039748*(x11_scaled[,"PFKM"])) +
  (-0.17829467*(x11_scaled[,"PRPS2"]))
risk_score <- cbind(x11_scaled, rs1)
range(rs1)
hrisk_i <- which(rs1 > 0)
lrisk_i <- which(rs1 < 0)
rs_genes <- c("ATP5F1D", "ATP5PB", "ATP5PF","NDUFA2", "NDUFA4L2", "NDUFB5", "TCIRG1", "UQCRH", "ALDOC",
              "G6PD","GPI", "PFKM", "PRPS2")
rs_gn_i <- c()
for (i in 1:length(rs_genes)) {
  rs_gn_i[i] = which(colnames(risk_score) == rs_genes[i])
}
ex_vald <- risk_score[,c(rs_gn_i,ncol(risk_score))]
lr_test_i_v <- which(ex_vald[,14] < 0)
hr_test_i_v <- which(ex_vald[,14] > 0)
lr_test_v <- as.data.frame(sort((ex_vald[,14])[lr_test_i_v], decreasing = F))
hr_test_v <- as.data.frame(sort((ex_vald[,14])[hr_test_i_v], decreasing = F))
colnames(lr_test_v) <- c("risk_score")
colnames(hr_test_v) <- c("risk_score")
order_v <- rbind(lr_test_v, hr_test_v)
order_index_v <- rownames(order_v)
val_data <- ex_vald[order_index_v,]
val_data_col <- c(rep("orchid2", times = length(lr_test_i_v )), rep("dodgerblue", times = length(hr_test_i_v )) )
val_data <- cbind(val_data, val_data_col)
colnames(val_data)[c(14,15)] <- c("risk_score", "colour")
plot((val_data[,14]), col = val_data[,15], ylab = "Risk Score", pch = 16,
     main = "Risk score distribution in External Validation set", cex.lab = 1.4)
abline(v = length(lr_test_i_v), lty = 2, col = "red", lwd = 2)
legend("bottomright", legend = c("LR", "HR"),              
       col = c("orchid2", "dodgerblue"),
       pch = 16, cex = 1.4,                              
       title = "Risk Group")
library(pheatmap)
heat_val_data <- val_data[,-c(14:15)]
heat_val_data <- t(heat_val_data)
#heat_val_data <- heat_val_data[rownames(val_data),]
negbreak <- seq(-3,0, by =0.05)
posbreak <- seq(0,3, by =0.05)
cbreak <- c(negbreak, posbreak[-1])
risk_annot <- c(rep("LR", times = length(lr_test_i_v)), rep("HR", times = length(hr_test_i_v)))
risk_annot <- as.data.frame(risk_annot)
rownames(risk_annot) <- colnames(heat_val_data)
colnames(risk_annot) <- c("Risk Group")
heat_val_data1 <- apply(heat_val_data, 2, as.numeric)
rownames(heat_val_data1) <- rownames(heat_val_data)
pheatmap((heat_val_data1), color = colorRampPalette(c("#d7f0f9","#76bce8", "#1c3f9a", "#000000", "#7a0e0e", "#e04b3a", "#fca283"))(length(cbreak)),
         fontsize_col = 3, fontfamily = "Times", fontsize = 11, cluster_rows = F, breaks = cbreak, fontsize_row = 11, cluster_cols = F,
         annotation_col = risk_annot, gaps_col = length(lr_test_i_v))

##tdROC curve in external validation##
validation <- ex_vald[,c(14)][rownames(val_data)]
y1 <- read.csv("survival_cptac.csv")
rownames(y1) <- y1[,3]
val_os <- y1[rownames(ex_vald),]
val_os <- na.omit(val_os)
validation <- validation[rownames(val_os)]
validation_data <- cbind(validation, val_os)
validation_data <- as.data.frame(validation_data)
colnames(validation_data) <- c("rs", "status", "time", "patient_id")
validation_data$status <- gsub("Dead", "1", validation_data$status)
validation_data$status <- gsub("Alive", "0", validation_data$status)
roc_v <- survivalROC(Stime = c(validation_data$time), status = c(validation_data$status), marker = c(validation_data$rs),
                     predict.time = 365, method = "KM")
plot(roc_v$FP, roc_v$TP, type = "l", col = "blue", xlab = "False Positive Rate",
     ylab = "True Positive Rate", main = "Time-dependent ROC for survival", cex.lab = 1.4)
abline(a = 0, b = 1, lty = 2, col = "red")
text(x = 0.85, y = 0.12, cex = 1.4, label = paste0("AUC = ", sprintf("%.4f",roc_v$AUC)))

