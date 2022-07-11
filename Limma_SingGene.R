#https://www.bilibili.com/video/BV1p34y1o7rv?vd_source=b2134ebf72ef3c8ede259339d95f9f98

library(limma)#加载包
options(stringsAsFactors = F)#设置为因子

#构建三个矩阵
##1 表达矩阵：符合要求的表达矩阵准备，详见GEO下载全集
library(GEOquery)
  studyID="GSE197671"
  Sys.setenv("VROOM_CONNECTION_SIZE"=131072*6)
  gse=getGEO(studyID, destdir="./", GSEMatrix=T, getGPL = F)
  exprSet <- exprs(gse[[1]])
  pdata <- pData(gse[[1]])
exprSet2 <- openxlsx::read.xlsx('./GSE197671_RNAseq_raw_counts_matrix.xlsx',rowNames = T)
colnames(exprSet2) <- gsub('.*: ','',pdata$title)

##2分组矩阵
grouplist <- factor(gsub(" ","_",pdata$`genotype/variation:ch1`))
design <- model.matrix(~0+grouplist)
rownames(design) = colnames(exprSet2)
colnames(design) <- levels(grouplist)
colnames(design)

##3差异比较矩阵
cont.matrix <- makeContrasts(contrasts = c("LINC00881_knockdown-control_overexpression"), levels = design)


#归一化
##方法1通过boxplot查看数据整齐与否再决定是否归一化
library(RColorBrewer)
dge <- DGEList(counts = exprSet2)
col <- brewer.pal(ncol(dge$counts), "Paired")
par(mfrow=c(2,2))
boxplot(dge$counts,outline=F, col=col)
title(main="A. Unnormalised ",ylab="raw count")
boxplot(calcNormFactors(dge, method = "TMM")$counts,outline=F,col=col)
title(main="B. TMM ",ylab="raw count")
boxplot(cpm(dge$counts),outline=F, col=col)
title(main="C. CPM ",ylab="cpm")
boxplot(cpm(dge$counts,log=TRUE),outline=F, col=col)
title(main="D. Log-CPM ",ylab="log-cpm")
#由上图可以发现log-cpm归一化效果最好，如果直接log-cpm化效果还不是很好，可以先通过TMM标准化再进行log-cpm化

##方法2：通过voom进行归一化
dge <- DGEList(counts = exprSet2)
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=TRUE) #会自动计算log(cpm)值


##方法3：normalizeBetweenArrays归一化
keep <- rowSums(cpm(exprSet2) > 1 ) >= 2 #过滤
exprSet2 <- exprSet2[keep, ]
norm<-normalizeBetweenArrays(exprSet2)
par(las=2)
boxplot(norm,outline=F, col=col)
par()

#差异分析
##拟合线性模型
fit <- lmFit(v, design)
##针对给定的对比计算估计系数和标准误差
fit2 <- contrasts.fit(fit, cont.matrix)
##计算出t统计量，F统计量和差异表达倍数的对数
fit2 <- eBayes(fit2)


options(digits = 4)

#卡阈值提取差异基因
allDEG <- topTable(fit2, coef =1, n = Inf) #与makecontrast对应
allDEG <- na.omit(allDEG)

##方法1
padj = 0.05
foldChange= 1
diff_signif = allDEG[(allDEG$adj.P.Val < padj & abs(allDEG$logFC)>foldChange),]                    
diff_signif = diff_signif[order(diff_signif$logFC),]
save(diff_signif, file = 'limma_diff.Rdata')

allDEG$gene <- row.names(allDEG)

#单基因表达量
library(ggplot2)
library(ggpubr) #pvalue
data <- data.frame(gene=norm['NAT10',],group=grouplist )

#https://zhuanlan.zhihu.com/p/496799780
#https://blog.51cto.com/u_15069485/3982986
#https://www.jianshu.com/p/d8e002de9fc3

ggboxplot(data=data,x="group",y="gene",color="group",add="jitter",legend="none",palette = "nejm")+
  rotate_x_text(angle = 45)+
  geom_hline(yintercept=mean(norm['NAT10',]),linetype=2)+
  stat_compare_means(method="anova",label.y = max(norm['NAT10',]*1.2))+ #add global annova p_value
  stat_compare_means(ref.group=".all.",
                   method = "t.test",
                   label = "p.signif") #pairwise comparsion against all
  
