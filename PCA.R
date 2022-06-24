#学习网站：https://www.jianshu.com/p/2f71bc493042；#https://qastack.cn/stats/35561/imputation-of-missing-values-for-pca
#正常seq矩阵的处理：确保matrix的numeric
rownames(hp_DIFF) <- hp_DIFF[,1]
hp_DIFF <- hp_DIFF[,-1]

#缺失值的处理
library(missMDA)
pc_res <- imputePCA(hp_DIFF,method="EM",ncp=1)

# 计算中值绝对偏差 (MAD, median absolute deviation)度量基因表达变化幅度
# 在基因表达中，尽管某些基因很小的变化会导致重要的生物学意义，
# 但是很小的观察值会引入很大的背景噪音，因此也意义不大。
# 可以选择根据mad值过滤，取top 50， 500等做分析
pc_res <- pc_res$completeObs
mads <- apply(pc_res, 1, mad)
pc_res <- pc_res[rev(order(mads)),]
pc_res <- pc_res[1:500,]

#matrix的转置，保证观测为样本，变量为打到的靶点
adj_pc_res=as.data.frame(t(pc_res))

#grouplist的构建
adj_pc_res$group=c("SXSM","SXSM","SXSM","CTRL","CTRL","CTRL")#group_list

#方法一
com1 <- prcomp(adj_pc_res[,1:(ncol(adj_pc_res)-1)], center = TRUE, scale. = TRUE)
summary(com1)
suppressPackageStartupMessages(library(ggfortify))
autoplot(com1, data=adj_pc_res,colour = 'group')

#方法二
suppressPackageStartupMessages(library(FactoMineR))#画主成分分析图需要加载这两个包
suppressPackageStartupMessages(library(factoextra)) 
  
dat.pca <- PCA(adj_pc_res[,1:(ncol(adj_pc_res)-1)], graph = FALSE)#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的
    pc <- fviz_pca_ind(dat.pca,
                       # geom.ind = "point", # show points only (nbut not "text")
                       col.ind = adj_pc_res$group, # color by groups
                       # palette = c("#00AFBB", "#E7B800"),
                       addEllipses = TRUE, # Concentration ellipses
                       legend.title = "Groups",
                       ellipse.level=0.95,
                       ellipse.type="confidence",
                       mean.point=F,
    )
    # 根据分组上色并绘制95%置信区间
    print(pc)
