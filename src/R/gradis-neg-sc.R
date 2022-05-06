library("caret")
library("data.table")
library("ggplot2")
library("here")
library("logger")
library("parallel")
library("readr")
library("ROCR")

here::i_am("src/R/gradis-neg-sc.R")

param.model_type <- "rf"

EXPfilename <- here::here("src/data/Expression_SC.txt")
GENEfilename <- here::here("src/data/Genes_SC.csv")
NETfilename <- here::here("src/data/Network_SC.csv")

log_info("Loading data.  This may take some time...")
T <- read_csv(EXPfilename, col_names = FALSE, show_col_types = FALSE)
X <- array(T)
EXPData <- read_csv(GENEfilename, col_names = FALSE, show_col_types = FALSE)
NETData <- read_csv(NETfilename, col_names = FALSE, show_col_types = FALSE)

ExpGenes <- EXPData$X1
NetTF <- NETData$X1
NetGenes <- NETData$X2
PosOrNeg <- NETData$X3
log_info("Done loading data.")

log_info("Initializing additional data structions...")
UniNetTF <- unique(NetTF)
UniGenes <- unique(NetGenes)
UniNetGene <- sort(union(UniGenes, UniNetTF))

C1 <- intersect(UniNetGene, ExpGenes)
iNetGenes <- match(C1, UniNetGene)
iExpGenes <- match(C1, ExpGenes)

XX <- X[iExpGenes,]

GeneNumber <- nrow(XX)
TimeNumber <- ncol(XX)

k <- floor((1 + sqrt(4 * TimeNumber + 1)) / 2)
XX_Clustered <- kmeans(t(XX), k)
XX_Clustered <- t(XX_Clustered$centers)

rowmax <- apply(XX_Clustered, 1, max)
XX_scaled <- XX_Clustered / rowmax

ExpGenes <- ExpGenes[iExpGenes]

C2 <- intersect(ExpGenes, UniNetTF)
iTF_Genes <- match(C2, ExpGenes)
iTF <- match(C2, UniNetTF)

TFsInExp <- ExpGenes[iTF_Genes]

inonTF_Genes <- setdiff(seq_along(iExpGenes), iTF_Genes)
inonTF_Genes <- t(inonTF_Genes)
nonTFsInExp <- ExpGenes[inonTF_Genes]
num_TF <- length(iTF_Genes)
num_nonTF <- length(inonTF_Genes)

log_info("Network initialization...")
## This could probably be vectorized
n <- length(ExpGenes)

NetPos <- matrix(0, n, n)
for(i in seq_along(NetTF))
  if(any(ExpGenes == NetTF[i]))
   if(any(ExpGenes == NetGenes[i]))
     if(PosOrNeg[i] == 1)
       NetPos[match(TRUE, ExpGenes == NetTF[i]), match(TRUE, ExpGenes == NetGenes[i])] <- 1

NetNeg <- matrix(0, n, n)
for(i in seq_along(NetTF))
  if(any(ExpGenes == NetTF[i]))
   if(any(ExpGenes == NetGenes[i]))
     if(PosOrNeg[i] == 0)
       NetNeg[match(TRUE, ExpGenes == NetTF[i]), match(TRUE, ExpGenes == NetGenes[i])] <- 1

NetworkPos <- NetPos[iTF_Genes,]
NetworkNeg <- NetNeg[iTF_Genes,]

r_Pos <- which(NetworkPos == 1, arr.ind = TRUE)[,1]
c_Pos <- which(NetworkPos == 1, arr.ind = TRUE)[,2]

r_Neg <- which(NetworkNeg == 1, arr.ind = TRUE)[,1]
c_Neg <- which(NetworkNeg == 1, arr.ind = TRUE)[,2]

Num_Features <- (ncol(XX_scaled) - 1) * ncol(XX_scaled) / 2

log_info("Feature construction...     (takes a while)")
Data_Pos <- matrix(0, length(r_Pos), Num_Features)
for(r in seq_along(r_Pos)) {
  X_pos <- XX_scaled[iTF_Genes[r_Pos[r]],]
  Y_pos <- XX_scaled[c_Pos[r],]

  x_pos_len <- length(X_pos)
  Dist_pos <- matrix(0, x_pos_len, x_pos_len)

  for(i in seq.int(1, length(X_pos) - 1))
      for(j in seq.int(i + 1, length(X_pos))) {
          Dist_pos[i, j] <- ((X_pos[i] - X_pos[j])^2 + (Y_pos[i] - Y_pos[j])^2)^(0.5)
          Dist_pos[j, i] <- ((X_pos[i] - X_pos[j])^2 + (Y_pos[i] - Y_pos[j])^2)^(0.5)
      }

  inx <- 1
  for(i in seq_along(ncol(Dist_pos))) {
    Data_Pos[r, seq.int(inx, inx + ncol(Dist_pos) - i - 1)] <-
      Dist_pos[i, seq.int(i + 1, ncol(Dist_pos))]
  }
}

Data_Neg <- matrix(0, length(r_Neg), Num_Features)
for(r in seq_along(r_Neg)) {
  X_neg <- XX_scaled[iTF_Genes[r_Neg[r]],]
  Y_neg <- XX_scaled[c_Neg[r],]

  x_neg_len <- length(X_neg)
  Dist_neg <- matrix(0, x_neg_len, x_neg_len)

  for(i in seq.int(1, length(X_neg) - 1))
      for(j in seq.int(i + 1, length(X_neg))) {
          Dist_neg[i, j] <- ((X_neg[i] - X_neg[j])^2 + (Y_neg[i] - Y_neg[j])^2)^(0.5)
          Dist_neg[j, i] <- ((X_neg[i] - X_neg[j])^2 + (Y_neg[i] - Y_neg[j])^2)^(0.5)
      }

  inx <- 1
  for(i in seq_along(ncol(Dist_neg))) {
    Data_Neg[r, seq.int(inx, inx + ncol(Dist_neg) - i - 1)] <-
      Dist_neg[i, seq.int(i + 1, ncol(Dist_neg))]
  }
}

log_info('Feature construction...     (Done)')

## several models for predicting negative labels
log_info('Prediction of negative samples...     (takes a while)')

rand_num_Pos <- sample.int(nrow(Data_Pos))
rand_num_Neg <- sample.int(nrow(Data_Neg))

model_count <- floor(length(rand_num_Neg) / length(rand_num_Pos))

log_info(model_count, ' Models are going to be trained for predicting negative samples')
Start_train_neg <- rep(0, model_count)
end_train_neg <- rep(0, model_count)
for(i in seq.int(model_count)) {
    Start_train_neg[i] <- 1 + length(rand_num_Pos) * (i - 1)
    end_train_neg[i] <- length(rand_num_Pos) * i
}

Labels <- matrix(0, length(rand_num_Neg),1)

fitControl <- trainControl(method = "cv", allowParallel = TRUE, verboseIter = TRUE)
for(i in seq_along(Start_train_neg)) {
  log_info('Model number ', i, ' for predicting negative samples.')

  Data_Train_A <- Data_Pos
  Data_Train_B <- Data_Neg[rand_num_Neg[Start_train_neg[i]:end_train_neg[i]],]
  Data_Train <- rbind(Data_Train_A, Data_Train_B)
  colnames(Data_Train) <- paste0("X", seq_len(ncol(Data_Train)))

  Group_Train <- as.factor(c(rep(1, length(rand_num_Pos)), rep(0, length(rand_num_Pos))))

  model <- train(Data_Train, Group_Train, method = param.model_type, preProcess = c("center", "scale"), trControl = fitControl)

  range_test_neg <- setdiff(seq_along(rand_num_Neg), rand_num_Neg[Start_train_neg[i]:end_train_neg[i]])
  Data_Test_Neg <- Data_Neg[range_test_neg,]
  colnames(Data_Test_Neg) <- paste0("X", seq_len(ncol(Data_Test_Neg)))

  label <- ifelse(predict(model, newdata = Data_Test_Neg, type = "raw") == 1, 1, 0)
  score <- predict(model, newdata = Data_Test_Neg, type = "prob")

  Labels[range_test_neg] <- label + Labels[range_test_neg]
}

log_info('Start training final model')
Neg_labels <- which(Labels == 0, arr.ind = FALSE)

Train_number <- floor(9 * length(r_Pos) / 10)
rand_num_Pos <- sample.int(nrow(Data_Pos))
rand_Neg <- sample.int(length(Neg_labels))

Data_Train_A <- Data_Pos[seq_len(Train_number),]
Data_Train_B <- Data_Neg[Neg_labels[rand_Neg[seq_len(Train_number)]],]
Data_Train <- rbind(Data_Train_A, Data_Train_B)
colnames(Data_Train) <- paste0("X", seq_len(ncol(Data_Train)))

Group_Train <- as.factor(c(rep(1, Train_number), rep(0, Train_number)))

model <- train(Data_Train, Group_Train, method = param.model_type, preProcess = c("center", "scale"), trControl = fitControl)

log_info('Test final model')
Data_Test_A <- Data_Pos[rand_num_Pos[seq.int(Train_number+1, length(rand_num_Pos))],]
Data_Test_B <- Data_Neg[Neg_labels[rand_Neg[seq.int(Train_number+1, length(r_Pos))]],]
Data_Test <- rbind(Data_Test_A, Data_Test_B)
colnames(Data_Test) <- paste0("X", seq_len(ncol(Data_Train)))

Test_number <- length(r_Pos) - Train_number
Group_Test <- as.factor(c(rep(1, Test_number), rep(0, Test_number)))

label <- ifelse(predict(model, newdata = Data_Test, type = "raw") == 1, 1, 0)
score <- predict(model, newdata = Data_Test, type = "prob")[,2]

accuracy <- sum(label == Group_Test) / length(Group_Test) * 100

pred.model <- ROCR::prediction(score, Group_Test)
perf.model <- ROCR::performance(pred.model, measure = "tpr", x.measure = "fpr")
auc <- ROCR::performance(pred.model, measure = "auc")
auc.model <- auc@y.values[[1]]

log_info("Plotting results")
model.label.model <- paste0(param.model_type, " (auc: ", round(auc.model, 2), ")")

df <- data.frame(x = 0:1 , y = 0:1)
roc.model <- data.frame(pfa = unlist(perf.model@x.values), pd = unlist(perf.model@y.values), model = model.label.model)

p <- ggplot() +
    geom_line(data = roc.model, aes(x=pfa, y=pd, color = model)) +
    xlab("False Positive Rate") + ylab("True Positive Rate") +
    geom_line(data = df, aes(x = x, y = y), linetype = "dotted") +
    labs(color = "Model") + theme(legend.position="bottom")
print(p)
ggsave(here::here("tmp/auc-perf-gradis-neg-sc-rf.pdf"), p, width = 30, height = 30, units = "cm")

