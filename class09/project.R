library(dplyr)

wisc.df <- read.csv("https://bioboot.github.io/bimm143_S18/class-material/WisconsinCancer.csv")

sum(wisc.df$diagnosis=="M")
grepl("_mean", colnames(wisc.df)) %>% sum(.)


wisc.data <- as.matrix(wisc.df[,3:32])
row.names(wisc.data) <- wisc.df$id
head(wisc.data)

diagnosis <- wisc.df$diagnosis

apply(wisc.data,2,sd)#sd
colMeans(wisc.data)#mean

wisc.pr <- prcomp(wisc.data,scale. = TRUE)
summary(wisc.pr)

biplot(wisc.pr)
plot(wisc.pr$x[,c(1)],wisc.pr$x[,2],col = diagnosis,
     xlab = "PCA1",ylab = "PCA2")

plot(wisc.pr$x[,c(1)],wisc.pr$x[,3],col = diagnosis,
     xlab = "PCA1",ylab = "PCA3")

pr.var <- wisc.pr$sdev^2
head(pr.var)
pve <- pr.var/sum(pr.var)
plot(pve,type = "o",xlab = "Principle Component",
     ylab = "Proportion of variance explained")

barplot(pve,names.arg = paste0("PC",c(1:length(pve))),las = 2,
        ylab = "Percent of variance explained",axes = FALSE)+
  axis(2, at=pve, labels=round(pve,2)*100 )
#wisc.pr
pve
?paste0
?barplot

require(factoextra)
fviz_eig(wisc.pr, addlabels = TRUE, ylim = c(0, 50))

wisc.pr$rotation[1,] %>% max(.)
wisc.pr$rotation["radius_mean",1]
sort( abs(wisc.pr$rotation[,1]) )


#hierachical clustering
data.scaled <- scale(wisc.data)
data.dist <- dist(data.scaled)
hclust(data.dist) %>% plot(.)+abline(h = 18,col = "red")
cutree(hclust(data.dist),k = 4)%>% plot(.)
wisc.hclust.clusters <- cutree(hclust(data.dist),k = 6)

table(wisc.hclust.clusters,diagnosis)

####combining method (PCA+hisrachical clusteirng)

#start with the PCs that capture 90% of the variance 
wisc.pr.hclust <- hclust(dist(wisc.pr$x[,1:7]),method="ward.D2")
plot(wisc.pr.hclust)
grps <- cutree(wisc.pr.hclust,k = 2)
table(grps)
table(greps,diagnosis)

plot(wisc.pr$x[,1],wisc.pr$x[,2],col = grps,
     xlab = "PCA1",ylab = "PCA2")
plot(wisc.pr$x[,c(1)],wisc.pr$x[,2],col = diagnosis,
     xlab = "PCA1",ylab = "PCA2")

#sensitivity
table(wisc.pr,diagnosis)
table(wisc.hclust.clusters,diagnosis)
table(greps,diagnosis)
#164/212 = 0.7736 Sensitivity
#(337)/(569-212) = 0.9439776 Specificity

# doing prediction
new_data <- read.csv("https://tinyurl.com/new-samples-CSV")
pred <- predict(wisc.pr,new_data)

plot(wisc.pr$x[,c(1)],wisc.pr$x[,2],col = diagnosis,
     xlab = "PCA1",ylab = "PCA2")+points(pred,col = "green",pch=16)
