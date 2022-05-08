library("caret")
library("data.table")
library("ggplot2")
library("here")
library("logger")
library("readr")
library("ROCR")

here::i_am("src/R/glmmer-ec.R")

log_info("Loading data positive and negative data")
genepos.dt <- data.table(read_csv(here::here("src/data/datapos.csv"), col_names = FALSE))
geneneg.dt <- data.table(read_csv(here::here("src/data/dataneg.csv"), col_names = FALSE))

log_info("Munging datasets")
genepos.dt[, y := 1]
geneneg.dt[, y := 0]
gene.dt <- rbind(genepos.dt, geneneg.dt)
gene.dt[, y := as.factor(y)]

log_info("Creating test/train split")
trainIndex <- createDataPartition(gene.dt$y, p = .8, list = FALSE, times = 1)
gene.trn.dt <- gene.dt[trainIndex]
gene.tst.dt <- gene.dt[-trainIndex]

fitControl <- trainControl(method = "cv", allowParallel = TRUE, verboseIter = TRUE)

model1 <- glm(as.factor(y) ~ ., data = gene.trn.dt, family = binomial())

log_info("Creating logit model")
model.name.model1 <- "glm"
model1 <- train(as.factor(y) ~ ., data = gene.trn.dt, method = model.name.model1, trControl = fitControl)

log_info("Creating randomforest model")
model.name.model2 <- "rf"
model2 <- train(as.factor(y) ~ ., data = gene.trn.dt, method = model.name.model2, trControl = fitControl)

log_info("Creating naive Bayes model")
model.name.model3 <- "naive_bayes"
model3 <- train(as.factor(y)  ~ ., data = gene.trn.dt, method = model.name.model3, trControl = fitControl)

log_info("Testing models against holdout data")
fitted.model1 <- predict(model1, newdata = gene.tst.dt, type = "prob")
fitted.model2 <- predict(model2, newdata = gene.tst.dt, type = "prob")
fitted.model3 <- predict(model3, newdata = gene.tst.dt, type = "prob")

pred.model1 <- ROCR::prediction(fitted.model1[, 2], gene.tst.dt[, y])
pred.model2 <- ROCR::prediction(fitted.model2[, 2], gene.tst.dt[, y])
pred.model3 <- ROCR::prediction(fitted.model3[, 2], gene.tst.dt[, y])

perf.model1 <- ROCR::performance(pred.model1, measure = "tpr", x.measure = "fpr")
perf.model2 <- ROCR::performance(pred.model2, measure = "tpr", x.measure = "fpr")
perf.model3 <- ROCR::performance(pred.model3, measure = "tpr", x.measure = "fpr")

auc <- ROCR::performance(pred.model1, measure = "auc")
auc.model1 <- auc@y.values[[1]]
auc <- ROCR::performance(pred.model2, measure = "auc")
auc.model2 <- auc@y.values[[1]]
auc <- ROCR::performance(pred.model3, measure = "auc")
auc.model3 <- auc@y.values[[1]]

log_info("Plotting results")
model.label.model1 <- paste(model.name.model1, " (auc: ", round(auc.model1, 2), ")", sep = "")
model.label.model2 <- paste(model.name.model2, " (auc: ", round(auc.model2, 2), ")", sep = "")
model.label.model3 <- paste(model.name.model3, " (auc: ", round(auc.model3, 2), ")", sep = "")

df <- data.frame(x = 0:1 , y = 0:1)
roc.model1 <- data.frame(pfa = unlist(perf.model1@x.values), pd = unlist(perf.model1@y.values), model = model.label.model1)
roc.model2 <- data.frame(pfa = unlist(perf.model2@x.values), pd = unlist(perf.model2@y.values), model = model.label.model2)
roc.model3 <- data.frame(pfa = unlist(perf.model3@x.values), pd = unlist(perf.model3@y.values), model = model.label.model3)

p <- ggplot() +
    geom_line(data = roc.model1, aes(x=pfa, y=pd, color = model)) +
    geom_line(data = roc.model2, aes(x=pfa, y=pd, color = model)) +
    geom_line(data = roc.model3, aes(x=pfa, y=pd, color = model)) +
    xlab("False Positive Rate") + ylab("True Positive Rate") +
    geom_line(data = df, aes(x = x, y = y), linetype = "dotted") +
    labs(color = "Model") + theme(legend.position="bottom")
print(p)
ggsave(here::here("tmp/auc-perf-ec.pdf"), p, width = 30, height = 30, units = "cm")

log_info("Saving data image")
save.image("tmp/glmmer-ec.Rdata")