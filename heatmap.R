#热图通常有两大作用：数据质量控制和直观展示重点研究对象的差异变化情况
#https://blog.csdn.net/qq_43210428/article/details/120020284
library(dplyr)
library(pheatmap)

#差异物聚类就可以了，无差异的聚类的背景噪音太大
map_DIFF_1 <- DIFF_1[which(DIFF_1$Type!='insig'),-c(2:3,4:12,16,20:23)] 
map_DIFF_1[,1] <- paste0(map_DIFF_1[,1],"_pro1")
map_DIFF_2 <- DIFF_2[-c(2:3,4:12,19)]
map_DIFF_2[,1] <- paste0(map_DIFF_2[,1],"_pro2")
map_DIFF <- rbind(map_DIFF_1,map_DIFF_2)

#规整为只有numeric的，sample为变量，靶点为观测
rownames(map_DIFF) <- map_DIFF[,1]
map_DIFF <- map_DIFF[,-1]
colnames(map_DIFF) <- c("SXSM1","SXSM2","SXSM3","CTRL1","CTRL2","CTRL3")

#处理缺失值
map_DIFF <- imputePCA(map_DIFF,method="EM",ncp=1)

#列名注释
annotation_col = data.frame(Sample=factor(c(rep("SXSM",3),rep("CTRL",3))))
row.names(annotation_col) = colnames(map_DIFF$completeObs)
ann_color <- list(a<-c(group1="blue", group2="red", group3="green"))

#画图
pheatmap(map_DIFF$completeObs, 
         # annotation_row=dfGene, # （可选）指定行分组文件
         annotation_col=annotation_col, # （可选）指定列分组文件
         annotation_colors = ann_color,
         show_colnames = TRUE, # 是否显示列名
         show_rownames=F,  # 是否显示行名
         #fontsize=2, # 字体大小
         color = colorRampPalette(c('#0000ff','#ffffff','#ff0000'))(50), # color 设置渐变的颜色，通常借助于colorRampPalette函数，比如说设置红黄蓝渐变，并在这之间分成50个等级，我们可以设置color=colorRampPalette(c("red","yellow","blue"))(50)
         annotation_legend=TRUE, # 是否显示图例
         border_color=NA,  # 边框颜色 NA表示没有
         scale="row",  # 指定归一化的方式。"row"按行归一化，"column"按列归一化，"none"不处理
         cluster_rows = TRUE, # 是否对行聚类
         cluster_cols = TRUE # 是否对列聚类
)
