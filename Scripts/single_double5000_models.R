library(randomForest)
library(caret)
library(ggplot2)

tt.path <- "/Users/ahorvath/Monika/mutations/HA4"
path <- "/Users/ahorvath/Monika/mutations/HA"
ff <- "single_wt_comb_3sd.txt"
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
tt <- read.table("/Users/ahorvath/Monika/mutations/HA2/hotreg_3sd_tox_diffcspDPS.txt", header = T)
stop <- F
for (ntree.index in 1:length(ntrees)) {
	for (mtry.index in 1:length(mtrys)) {
		ntree <- ntrees[ntree.index]
		mtry <- mtrys[mtry.index]
		print(ntree)
		print(mtry)
		set.seed(123)
		rf.model3 <- randomForest(TOX ~ DCS + DpDP + DS, data = merged.tt.train, ntree = ntree, mtry = mtry, importance = T)
		print(cor(predict(rf.model3, merged.tt.test), merged.tt.test$TOX))
		print(cor(predict(rf.model3, g1), g1$TOX))
		print(cor(predict(rf.model3, tt), tt$TOX))
		#if ((cor(predict(rf.model3, tt), tt$TOX) > 0.72) && (cor(predict(rf.model3, ), tt$TOX) > 0.82)) {
		#	stop = T
		#	break
		#}
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

ntree <- 300
mtry <- 2
                set.seed(123)
		merged.tt <- rbind(g1[,c("TOX","DCS","DpDP","DS")],tt[sample(1:nrow(tt), 5000), c("TOX","DCS","DpDP","DS")])
		group.size <- nrow(merged.tt)
		train.group.size<- group.size * 0.7
		train.indices <- sample(1:group.size, train.group.size)
		merged.tt.train <- merged.tt[train.indices,]
		merged.tt.test <- merged.tt[-train.indices,]
                rf.model3 <- randomForest(TOX ~ DCS + DpDP + DS, data = merged.tt.train, ntree = ntree, mtry = mtry, importance = T)
                #set.seed(123)
                #rf.model6 <- randomForest(TOX ~ DCS + DpDP + DS + CS + pDP + S, data = g1, ntree = ntree, mtry = mtry, importance = T)
                print(cor(predict(rf.model3, merged.tt.train), merged.tt.train$TOX))
                print(cor(predict(rf.model3, merged.tt.test), merged.tt.test$TOX))
                print(cor(predict(rf.model3, tt), tt$TOX))
                print(cor(predict(rf.model3, g1), g1$TOX))
               # print(cor(predict(rf.model6, g1), g1$TOX))
               # print(cor(predict(rf.model6, tt), tt$TOX))

f <- "/Users/ahorvath/Monika/mutations/HA2/hotreg_3sd_tox_als_diffcspDPS.txt"
tt2 <- read.table(f, header=T, stringsAsFactors=F)
print(cor(predict(rf.model3, tt2), tt2$TOX))
#print(cor(predict(rf.model6, tt2), tt2$TOX))

f3 <- "/Users/ahorvath/Monika/mutations/HA/single_wt_comb_3sd_als.txt"
tt3 <- read.table(f3, header=T, stringsAsFactors=F)
print(cor(predict(rf.model3, tt3), tt3$TOX)) 
#print(cor(predict(rf.model6, tt3), tt3$TOX)) 
nn <- paste0(gsub(".txt", "", ff), "_allpars")
nn <- "single_double5000"
save(rf.model3, file = paste0(nn, "_rf_ntree300.RData"))

gg.cor.rf.model3 <- ggplot(data.frame(predicted_TOX = predict(rf.model3, merged.tt.test), real_TOX = merged.tt.test$TOX), aes(real_TOX, predicted_TOX)) +  geom_point(size=1, alpha = 0.6) + theme_bw()
ggsave(gg.cor.rf.model3, filename = paste0(nn,"_gg.rf3.corplot.pdf"))

gg.cor.rf.model3.merged.test <- ggplot(data.frame(predicted_TOX = predict(rf.model3, merged.tt.test), real_TOX = merged.tt.test$TOX)) + geom_point(aes(real_TOX, predicted_TOX, colour = real_TOX), size = 1, alpha = 0.85) + theme_bw() +theme(axis.text.x = element_text(face="bold", size=16), axis.text.y = element_text(face="bold", size=16), aspect.ratio=0.65) + theme (panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +scale_colour_gradient(low = "green", high = "red") + xlim(-0.3,0.3) + ylim(-0.3,0.3)
ggsave(gg.cor.rf.model3.merged.test, filename = paste0(nn,"_gg.rf3.merged.test.corplot.pdf"))

gg.cor.rf.model3.single <- ggplot(data.frame(predicted_TOX = predict(rf.model3, g1), real_TOX = g1$TOX)) + geom_point(aes(real_TOX, predicted_TOX, colour = real_TOX), size = 1, alpha = 0.85) + theme_bw() +theme(axis.text.x = element_text(face="bold", size=16), axis.text.y = element_text(face="bold", size=16), aspect.ratio=0.65) + theme (panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +scale_colour_gradient(low = "green", high = "red") + xlim(-0.3,0.3) + ylim(-0.3,0.3)
ggsave(gg.cor.rf.model3.single, filename = paste0(nn,"_gg.rf3.single.corplot.pdf"))

gg.cor.rf.model3.double <- ggplot(data.frame(predicted_TOX = predict(rf.model3, tt), real_TOX = tt$TOX)) + geom_point(aes(real_TOX, predicted_TOX, colour = real_TOX), size = 1, alpha = 0.85) + theme_bw() +theme(axis.text.x = element_text(face="bold", size=16), axis.text.y = element_text(face="bold", size=16), aspect.ratio=0.65,) + theme (panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +scale_colour_gradient(low = "green", high = "red") + xlim(-0.3,0.3) + ylim(-0.3,0.3)
ggsave(gg.cor.rf.model3.double, filename = paste0(nn,"_gg.rf3.double.corplot.pdf"))


gg.cor.rf.model3.single.als <- ggplot(data.frame(predicted_TOX = predict(rf.model3, tt3), real_TOX = tt3$TOX)) + geom_point(aes(real_TOX, predicted_TOX, colour = real_TOX), size = 1, alpha = 0.85) + theme_bw() +theme(axis.text.x = element_text(face="bold", size=16), axis.text.y = element_text(face="bold", size=16), aspect.ratio=0.65) + theme (panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +scale_colour_gradient(low = "green", high = "red") + xlim(-0.3,0.3) + ylim(-0.3,0.3)
ggsave(gg.cor.rf.model3.single.als, filename = paste0(nn,"_gg.rf3.single.als.corplot.pdf"))

gg.cor.rf.model3.double.als <- ggplot(data.frame(predicted_TOX = predict(rf.model3, tt2), real_TOX = tt2$TOX)) + geom_point(aes(real_TOX, predicted_TOX, colour = real_TOX), size = 1, alpha = 0.85) + theme_bw() +theme(axis.text.x = element_text(face="bold", size=16), aspect.ratio=0.65, axis.text.y = element_text(face="bold", size=16)) + theme (panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +scale_colour_gradient(low = "green", high = "red") + xlim(-0.3,0.3) + ylim(-0.3,0.3)
ggsave(gg.cor.rf.model3.double.als, filename = paste0(nn,"_gg.rf3.double.als.corplot.pdf"))




f.single.drop <- "/Users/ahorvath/Monika/mutations/HA6/single_3sd_drop_label.txt"
t.single.drop <- read.table(f.single.drop, header=F, stringsAsFactors=F)
colnames(t.single.drop) <- c("ID", "MUT", "TOX", "SD", "DCS",  "DpDP",  "DS",  "CS",  "pDP", "S", "Label")
print(cor(predict(rf.model3, t.single.drop), t.single.drop$TOX))

f.double.drop <- "/Users/ahorvath/Monika/mutations/HA6/double_3sd_drop_label.txt"
t.double.drop <- read.table(f.double.drop, header=F, stringsAsFactors=F)
colnames(t.double.drop) <- c("ID1", "ID2", "MUT", "TOX", "SD", "M2", "DCS", "DpDP", "DS", "CS", "pDP","S", "Label")
print(cor(predict(rf.model3, t.double.drop), t.double.drop$TOX))


