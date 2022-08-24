library(randomForest)
library(e1071)
library(caret)
library(ggplot2)
set.seed(1234)

tt.path <- "mutations/HA5"
path <- "mutations/HA2"
ff <- "hotreg_3sd_tox_diffcspDPS.txt"
g1 <- read.table(file.path(path, ff), header=T, stringsAsFactors=F)
nn <- paste0(gsub(".txt", "", ff), "_allpars")
print(nn)
dir.create(file.path(tt.path, nn), showWarnings = F)
setwd(file.path(tt.path, nn))
ntrees <- seq(from = 50, to = 3000, by = 50)
mtrys <- seq(from = 1, to = 15, by = 1)
model3 <- matrix(NA,length(ntrees), length(mtrys))
model6 <- matrix(NA,length(ntrees), length(mtrys))
rownames(model3) <- rownames(model6) <- ntrees
colnames(model3) <- colnames(model6) <- mtrys
model3.list <- list()
model6.list <- list()
tt <- read.table("mutations/HA/single_wt_comb_3sd.txt", header = T)
stop <- F
for (ntree.index in 1:length(ntrees)) {
	for (mtry.index in 1:length(mtrys)) {
		ntree <- ntrees[ntree.index]
		mtry <- mtrys[mtry.index]
		print(ntree)
		print(mtry)
		set.seed(123)
		rf.model3 <- randomForest(TOX ~ DCS + DpDP + DS, data = g1, ntree = ntree, mtry = mtry, importance = T)
		set.seed(123)
		rf.model6 <- randomForest(TOX ~ DCS + DpDP + DS + CS + pDP + S, data = g1, ntree = ntree, mtry = mtry, importance = T)
		print(cor(predict(rf.model3, g1), g1$TOX))
		print(cor(predict(rf.model3, tt), tt$TOX))
		print(cor(predict(rf.model6, g1), g1$TOX))
                print(cor(predict(rf.model6, tt), tt$TOX))
		if ((cor(predict(rf.model3, tt), tt$TOX) > 0.8) && (cor(predict(rf.model6, tt), tt$TOX) > 0.8)) {
			stop = T
			break
		}
	#	model3[ntree.index, mtry.index] <- cor(predict(rf.model3, tt), tt$TOX)
	#	model6[ntree.index, mtry.index] <- cor(predict(rf.model6, tt), tt$TOX)
		#model3.list[paste0(ntree, "_", mtry.index)] <- rf.model3
		#model6.list[paste0(ntree, "_", mtry.index)] <- rf.model6
		#rf.model1 <- randomForest(TOX ~ DCS + DpDP + DS + CS + pDP + S, data = g1, ntree = 500, mtry = 6, importance = T)
	}
	if (stop) {
		break
	}
}
save.image()
ntree <- 300
mtry <- 2
                set.seed(123)
                rf.model3 <- randomForest(TOX ~ DCS + DpDP + DS, data = g1, ntree = ntree, mtry = mtry, importance = T)
                #set.seed(123)
                #rf.model6 <- randomForest(TOX ~ DCS + DpDP + DS + CS + pDP + S, data = g1, ntree = ntree, mtry = mtry, importance = T)
                print(cor(predict(rf.model3, g1), g1$TOX))
                print(cor(predict(rf.model3, tt), tt$TOX))
                #print(cor(predict(rf.model6, g1), g1$TOX))
                #print(cor(predict(rf.model6, tt), tt$TOX))

f <- "mutations/HA2/hotreg_3sd_tox_als_diffcspDPS.txt"
tt2 <- read.table(f, header=T, stringsAsFactors=F)
print(cor(predict(rf.model3, tt2), tt2$TOX))
#print(cor(predict(rf.model6, tt2), tt2$TOX))

f3 <- "mutations/HA/single_wt_comb_3sd_als.txt"
tt3 <- read.table(f3, header=T, stringsAsFactors=F)
print(cor(predict(rf.model3, tt3), tt3$TOX)) 
#print(cor(predict(rf.model6, tt3), tt3$TOX)) 
nn <- gsub(".txt", "", ff)



save(rf.model3, file = paste0(nn, "_rf_ntree300.RData"))

gg.cor.rf.model3 <- ggplot(data.frame(predicted_TOX = predict(rf.model3, g1), real_TOX = g1$TOX), aes(real_TOX, predicted_TOX)) +  geom_point(size=1, alpha = 0.6) + theme_bw()
ggsave(gg.cor.rf.model3, filename = paste0(nn,"_gg.rf3.corplot.pdf"))


#gg.cor.rf.model6 <- ggplot(data.frame(predicted_TOX = predict(rf.model6, g1), real_TOX = g1$TOX), aes(real_TOX, predicted_TOX)) +  geom_point(size=1, alpha = 0.6) + theme_bw()
#ggsave(gg.cor.rf.model6, filename = paste0(nn,"_gg.rf6.corplot.pdf"))

gg.cor.rf.model3.tt <- ggplot(data.frame(predicted_TOX = predict(rf.model3, tt), real_TOX = tt$TOX), aes(real_TOX, predicted_TOX)) +  geom_point(size=1, alpha = 0.6) + theme_bw()
ggsave(gg.cor.rf.model3.tt, filename = paste0(nn,"_gg.rf3.corplot_on_single_wt_comb_3sd.pdf"))

#gg.cor.rf.model6.tt <- ggplot(data.frame(predicted_TOX = predict(rf.model6, tt), real_TOX = tt$TOX), aes(real_TOX, predicted_TOX)) +  geom_point(size=1, alpha = 0.6) + theme_bw()
ggsave(gg.cor.rf.model6.tt, filename = paste0(nn,"_gg.rf6.corplot_on_single_wt_comb_3sd.pdf"))


########################3


gg.cor.rf.model3.tt2 <- ggplot(data.frame(predicted_TOX = predict(rf.model3, tt2), real_TOX = tt2$TOX), aes(real_TOX, predicted_TOX)) +  geom_point(size=1, alpha = 0.6) + theme_bw()
ggsave(gg.cor.rf.model3.tt2, filename = paste0(nn,"_gg.rf3.corplot_on_hotreg_3sd_tox_als_diffcspDPS.pdf"))

gg.cor.rf.model6.tt2 <- ggplot(data.frame(predicted_TOX = predict(rf.model6, tt2), real_TOX = tt2$TOX), aes(real_TOX, predicted_TOX)) +  geom_point(size=1, alpha = 0.6) + theme_bw()
ggsave(gg.cor.rf.model6.tt2, filename = paste0(nn,"_gg.rf6.corplot_on_hotreg_3sd_tox_als_diffcspDPS.pdf"))



rf.model3.imp1 <- rf.model3.imp2 <- as.data.frame(importance(rf.model3))
rf.model3.imp1 <- rf.model3.imp1[order(rf.model3.imp1[,1], decreasing = T),]
rf.model3.imp1$Variable <- rownames(rf.model3.imp1)
rf.model3.imp1$Variable <- factor(rf.model3.imp1$Variable, levels = rf.model3.imp1$Variable)
colnames(rf.model3.imp1) <- c("perc_IncMSE", "IncNodePurity", "Variable")
gg.rf.model3.IncMSE <- ggplot(rf.model3.imp1, aes(Variable, perc_IncMSE)) + geom_bar(stat = "identity") + theme_bw()
ggsave(gg.rf.model3.IncMSE, filename = paste0(nn, "_gg.rf.model3.incMSE.pdf"))

rf.model3.imp2 <- rf.model3.imp2[order(rf.model3.imp2[,2], decreasing = T),]
rf.model3.imp2$Variable <- rownames(rf.model3.imp2)
rf.model3.imp2$Variable <- factor(rf.model3.imp2$Variable, levels = rf.model3.imp2$Variable)
colnames(rf.model3.imp2) <- c("perc_IncMSE", "IncNodePurity", "Variable")
gg.rf.model3.IncNodePurity <- ggplot(rf.model3.imp2, aes(Variable, IncNodePurity)) + geom_bar(stat = "identity")  + theme_bw()
ggsave(gg.rf.model3.IncNodePurity, filename = paste0(nn, "_gg.rf.model3.incNodePurity.pdf"))

rf.model6.imp1 <- rf.model6.imp2 <- as.data.frame(importance(rf.model6))
rf.model6.imp1 <- rf.model6.imp1[order(rf.model6.imp1[,1], decreasing = T),]
rf.model6.imp1$Variable <- rownames(rf.model6.imp1)
rf.model6.imp1$Variable <- factor(rf.model6.imp1$Variable, levels = rf.model6.imp1$Variable)
colnames(rf.model6.imp1) <- c("perc_IncMSE", "IncNodePurity", "Variable")
gg.rf.model6.IncMSE <- ggplot(rf.model6.imp1, aes(Variable, perc_IncMSE)) + geom_bar(stat = "identity") + theme_bw()
ggsave(gg.rf.model6.IncMSE, filename = paste0(nn, "_gg.rf.model6.incMSE.pdf"))

rf.model6.imp2 <- rf.model6.imp2[order(rf.model6.imp2[,2], decreasing = T),]
rf.model6.imp2$Variable <- rownames(rf.model6.imp2)
rf.model6.imp2$Variable <- factor(rf.model3.imp2$Variable, levels = rf.model6.imp2$Variable)
colnames(rf.model6.imp2) <- c("perc_IncMSE", "IncNodePurity", "Variable")
gg.rf.model6.IncNodePurity <- ggplot(rf.model6.imp2, aes(Variable, IncNodePurity)) + geom_bar(stat = "identity")  + theme_bw()
ggsave(gg.rf.model6.IncNodePurity, filename = paste0(nn, "_gg.rf.model6.incNodePurity.pdf"))
