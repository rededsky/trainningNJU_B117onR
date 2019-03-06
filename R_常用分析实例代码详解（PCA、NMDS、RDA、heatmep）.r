## 养成好的习惯，像写英文文章一样，单词、标点规范使用，便于后续阅读。不及时敲空格是本人自学时偷懒的坏毛病，请不要follow！！！
## 画图推荐用ggplot2及在该系统下的包与函数，其利用图层的方法很好，注意掌握数据绘图和非数据绘图的分割。R的入门内容多，需要记忆的内容
##很杂，这里就不细讲ggplot2的方法。建议慢慢自学。

######### 杜京京Environmental Science Nano的部分数据为例，展示NMDS和heatmap的做法 #############

setwd("D:/南大/za/dujj")	## 设置工作目录。注意windows的“\”要更改。
library(vegan)			## 加载用于群落分析的vegan包。其关键函数的详细用法可见帮助文件。
phlm<-read.csv('phlm.csv')	## 读入数据细菌门的群落数据文件
mds_phlm<-metaMDS(phlm[,-1],distance="bray")	## NMDS排序，距离系数采用bray-cirtis距离


plot(mds_phlm,display='site')	## 画图，采用默认参数。


scores(mds_phlm)	## 计算所得的排序值


plot(mds_phlm,display='site',type='n')		## 一步步画图，先画个空图
xz<-c(21,19,22,23,24,8)				## 自定义点的性状参数
points(mds_phlm,display="sites",pch=xz,cex=2)	## 画出“sites”的各点
legend("topright",pch=xz,legend=levels(phlm$trt),bty='n')	## 添加图例。（见Du et al. 2017的图7a）


phlm_corMDS<-data.frame(corMDS1=rep(0,5),P_corMDS1=rep(0,5),corMDS2=rep(0,5),P_corMDS2=rep(0,5),row.names=colnames(phlm[,-1]))
for(i in 1:5){
phlm_corMDS$corMDS1[i]<-cor(phlm[,-1][,i],scores(mds_phlm)[,1])
phlm_corMDS$P_corMDS1[i]<-cor.test(phlm[,-1][,i],scores(mds_phlm)[,1])$p.value
phlm_corMDS$corMDS2[i]<-cor(phlm[,-1][,i],scores(mds_phlm)[,2])
phlm_corMDS$P_corMDS2[i]<-cor.test(phlm[,-1][,i],scores(mds_phlm)[,2])$p.value
}	## 计算各分类单元（这里是某一个门）与两个排序轴的相关系数及显著性。（见图7a）
write.csv(phlm_corMDS,'phlm_corMDS.csv')


library(gplots)
library(RColorBrewer)
scalebulered<-colorRampPalette(c("deepskyblue","red"),space ="rgb")(100)	## 自定义一个颜色变化区间
row.names(phlm)<-phlm$trt
phlm<-phlm[,-1]
phlm.drp.dist<-vegdist(phlm,method="bray")		## 计算bray距离
phlmcol.clus<-hclust(phlm.drp.dist, "aver")		## 根据距离对处理进行组间聚类
phlm.drp.t.dist<-vegdist(t(phlm),method="bray")		## 转置数据矩阵后计算距离
phlmrow.clus<-hclust(phlm.drp.t.dist,"aver")		## 对物种进行聚类
heatmap.2(t(as.matrix(phlm)),Colv=as.dendrogram(phlmcol.clus),Rowv=as.dendrogram(phlmrow.clus),col=scalebulered,margins=c

(8,16),trace="none",density.info="none",lhei=c(2,8),scale='row')	## 绘制热图。（见Du et al. 2017的图6）


########## 孔第二篇文章的数据，水体酶活变异受水质的影响，展示PCA、RDA的做法 #######
library(vegan)
dat<-read.csv('RDA分析.csv')
dat$处理<-as.factor(dat$处理)	## 将处理变量定义为分类变量
Ez<-dat[,2:9]		## 酶活数据矩阵，可用于做PCA
Ev<-dat[,10:15]		## 水质数据矩阵，作为环境变量对酶活进行解释的RDA分析

decorana(Ez)		# 一般而言，最大值<3, 用PCA而非CA（继而限制性排序时用RDA非CCA）
EzPCA<-rda(Ez, scale=T)	# PCA分析。画图与前面NMDS类似。
plot(EzPCA)		# 默认绘图
biplot(EzPCA,display="sp",xlim=c(-2,2),ylim=c(-2,1))
points(EzPCA,display="sites",pch=dat$处理+14,cex=2,col=dat$处理+1)

rda_raw_scaled<-rda(Ez~.,Ev,scale=T)	## RDA分析，结果存储在“rda_raw_scaled”中
coef(rda_raw_scaled)			## 各系数展示
RsquareAdj(rda_raw_scaled)$r.squared	## 环境因子对酶活变异的解释力
anova.cca(rda_raw_scaled, step=1000)	## RDA模型的显著性
anova.cca(rda_raw_scaled, by="axis", step=1000)	## 各轴的显著性
rda_raw_scaled$CA$eig[rda_raw_scaled$CA$eig > mean(rda_raw_scaled$CA$eig)] ## 相对重要的轴（可作为参考）
summary(rda_raw_scaled)			## RDA分析所得的所有数据

windows(title="RDA scaling 2 + wa")	## 画图
plot(rda_raw_scaled,main="xxxx",type='n')
points(rda_raw_scaled,display='wa',pch=dat$处理+20,cex=1.4,col='darkgrey')
text(rda_raw_scaled,display='sp',col=1,cex=.8)
text(rda_raw_scaled,display='cn',col=1,cex=.8,lwd=2)
arrows(0, 0, spe2.sc[, 1], spe2.sc[, 2], length=0, lty=1, col=1)


scores(rda_raw_scaled)	## 一些可能有用的数据的提取
rda_raw_scaled$CCA$u.eig
rda_raw_scaled$CCA$u
rda_raw_scaled$CCA$wa
rda_raw_scaled$CCA$wa.eig
rda_raw_scaled$CCA$biplot







