set.seed(2020213)
library(dplyr)
library(Hmisc)
library(reshape2)
#代谢组准备
#差异代谢物矩阵
daixiezu_kegg <- read.csv("./dxz_keggout.path.csv")
daixiezu_kegg_ID <- daixiezu_kegg$Genes[1:20]
daixiezu_kegg_ID <- strsplit(daixiezu_kegg_ID, split = ";")
daixiezu_kegg_ID <- unlist(daixiezu_kegg_ID)
daixiezu_kegg_ID <- unique(daixiezu_kegg_ID)
daixiezu_kegg_ID <- as.data.frame(daixiezu_kegg_ID)
colnames(daixiezu_kegg_ID) <- 'ID'
select <- ALL[c( "ID" ,"Name","CTRL-1","CTRL-2","CTRL-3","SXSM-1","SXSM-2","SXSM-3")]
daixie_net <- merge(daixiezu_kegg_ID,select, by='ID')
daixie_net <- daixie_net[-1]
row.names(daixie_net)<-daixie_net$Name
daixie_net <- daixie_net[-1]
colnames(daixie_net) <- c("CTRL1","CTRL2","CTRL3","SXSM1","SXSM2","SXSM3")
#代谢物全用，蛋白组要挑

#代谢node文件
node <- row.names(daixie_net)
hz_dxz_node <- data.frame(node,class=NA)
i=1
 for(i in 1:nrow(hz_dxz_node)){
  hz_dxz_node$class[i] <- ALL$Class[which(ALL$Name==hz_dxz_node$node[i])]
  }
hz_dxz_node$class[is.na(hz_dxz_node$class)] <- 'unknown'

library(org.Hs.eg.db)
library(clusterProfiler)

#蛋白node文件
danbai_kegg <- read.csv("./dbz_keggout.path.csv")
danbai_kegg_ID <- danbai_kegg[c(2,4:6,9,10,12,14:20),c("Genes","KEGG_B_class")]
danbai_kegg_ID$Genes <- strsplit(danbai_kegg_ID$Genes,';')
danbai_kegg_ID[2,2] <- "Energy metabolism"
danbai_kegg_ID[5,2] <- "Energy metabolism"
danbai_kegg_ID[9,2] <- "Cardiac muscle contraction"

hz_dbz_node <- data.frame(node=unique(unlist(danbai_kegg_ID[c("4","5","9"),1])),class="Energy metabolism")
t <- data.frame(node=unique(unlist(danbai_kegg_ID[c("10","16","17","19"),1])),class="Amino acid metabolism")
hz_dbz_node <- rbind(hz_dbz_node, t)

t <- data.frame(node=unique(unlist(danbai_kegg_ID[c("15"),1])),class="Cardiac muscle contraction")
hz_dbz_node <- rbind(hz_dbz_node, t)

t <- data.frame(node=unique(unlist(danbai_kegg_ID[c("20"),1])),class="Lipid metabolism")
hz_dbz_node <- rbind(hz_dbz_node, t)

t <- data.frame(node=unique(unlist(danbai_kegg_ID[c("6","12","14"),1])),class="Translation")
hz_dbz_node <- rbind(hz_dbz_node, t)

t <- data.frame(node=unique(unlist(danbai_kegg_ID[c("2"),1])),class="Transcription")
hz_dbz_node <- rbind(hz_dbz_node, t)

t <- data.frame(node=unique(unlist(danbai_kegg_ID[c("18"),1])),class="Signal transduction")
hz_dbz_node <- rbind(hz_dbz_node, t)

#去重，保留关注度大的
hz_dbz_node_v <- hz_dbz_node[-c(3,4,5,34,51,53),]
duplicated(hz_dbz_node_v)

node <- bitr(hz_dbz_node_v$node, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
hz_dbz_node_v[1] <- node$SYMBOL

danbai_diff_col6 <- read.csv('./dbz_diff_col6.csv')
colnames(danbai_diff_col6)[1] <- "node" 
danbai_net <- merge(hz_dbz_node_v,danbai_diff_col6, by='node')
row.names(danbai_net) <- danbai_net$node
danbai_net <- danbai_net[-c(1:2)]


#1 蛋白自身互作
library(reshape2)
danbai_net_t <- t(danbai_net)
db_db <- Hmisc::rcorr(danbai_net_t, type="spearman") 
db_db_c <- db_db$r
db_db_p <- db_db$P
db_db_c <-  melt(data = db_db_c, id.vars = c("X")) ## 转化r矩阵
db_db_p <-  melt(data = db_db_p, id.vars = c("X")) ## 转化r矩阵
db_db_net <- db_db_c
db_db_net$p <-db_db_p$value
colnames(db_db_net) <- c("node_1","node_2","r","p")
IS_node_differ <-  db_db_net$node_1 != db_db_net$node_2 ### 选出左右不相等的行
db_db_net <-  db_db_net[IS_node_differ,] ### 提取出左右不相等的行

#2 蛋白代谢互作
library(psych)
daixie_net_t <- t(daixie_net)
db_dx <-  corr.test(danbai_net_t, daixie_net_t, use = "pairwise", method = "spearman", adjust = "BH")

#获取关联cor,pvalue,FDR值
db_dx_c <- melt(db_dx$r)
db_dx_p <- melt(db_dx$p)
db_dx_fdr <- melt(db_dx$p.adj)

db_dx_net <- db_dx_c
db_dx_net$p <- db_dx_p$value
colnames(db_dx_net) <- c("node_1","node_2","r","p")

#总图
net_all <- rbind(db_db_net,db_dx_net)

#边的颜色
net_all$color <- "blue" ### 将默认颜色设定为蓝色

net_all$color[net_all$r>= 0] <- "red"  ### 将r值>= 0的边的color属性赋值为"red"，而r值< 0的边的color属性赋值依然为"blue"，这样可以通过颜色区分正相关和负相关。

#线宽
net_all$linewidth <-  -1*log(net_all$p) ### 使用-logP来代表线的宽度
### 计算一个 fdr校正后的p值
library(dplyr)
net_all$fdr <-  net_all$p %>% p.adjust(.,method = "fdr")
### 根据选定标准对数据集进行筛选，使用“&”符号对 p、r的判定结果取交集。

IS_edge_qualified <- net_all$p < 0.05 & abs(net_all$r) >= 0.7## 若想要用fdr来进行筛选，可以把这里的 edge_df$p < p_cutoff替换为：edge_df$fdr < p_cutoff

net_all <-  net_all[IS_edge_qualified,]
which(abs(net_all$r)==1)
net_all_fn <- net_all[(!(abs(net_all$r)==1)),]
str(net_all_fn)
net_all_fn$node_2 <-  net_all_fn$node_2 %>% as.character() ## 将node_2列转换为字符形式
net_all_fn$node_1 <-  net_all_fn$node_1 %>% as.character() ## 将node_1列转换为字符形式
write.csv(net_all_fn,"02-1 edge file.csv") ## row.names = F

#节点的处理
#蛋白组
hz_dbz_node_v$Shape <- 'Ellipse'

library(RColorBrewer)
display.brewer.all() #显示所有调色板
display.brewer.all(type = "qual")
hz_dbz_node_v$class %>% unique() -> N
hz_dbz_node_v$class %>% unique() %>% length() -> n
color <- colorRampPalette(brewer.pal(8,"Accent"))(n)
hz_dbz_node_v$color <- hz_dbz_node_v$class %>% factor(., levels=N, labels = color) ### 将3种分类的代谢物颜色属性分别设置为橙色、青绿色和砖红色。

#绘制供图注的颜色图形
height<- c(rep(1,n))
barplot(height,
          col=color,
          names.arg=N,
        legend.text=N
        )

#代谢组
hz_dxz_node$Shape <- 'Rectangle'
hz_dxz_node$class %>% unique() -> M
hz_dxz_node$class %>% unique() %>% length() -> m
COLOR <- brewer.pal(m, "Greys")
hz_dxz_node$color <- hz_dxz_node$class %>% factor(., levels=M, labels = COLOR) ### 将3种分类的代谢物颜色属性分别设置为橙色、青绿色和砖红色。

node_all <- rbind(hz_dbz_node_v,hz_dxz_node)

#计算SIZE
## 以各个节点连接的边的数量作为权重控制节点的大小。一般而言，网络图中节点和其他节点的连接数越多，可以认为这个节点对于该网络越重要，各个节点的连接边的数量可以通过统计节点在边文件（edge_df）数据两个端点（node_1列和node_2列）出现的次数来进行统计。

all_counts <-  c(net_all_fn$node_1,net_all_fn$node_2) ##将边文件中所有节点信息合并，用于统计每个节点的出现次数。

## 对每个数据框中的节点，依次查看其在all_nodes中出现的次数（在edge file 中出现的次数），并存放到“size”列中。

for(i in 1:nrow(node_all)){

    # if(i == 1){node_frequency_v <- NULL }

    each_node_frequency <-  (all_counts  == node_all[i,"node"] )%>% sum(.,na.rm = T)

    node_all[i,"size"] <-  each_node_frequency}
write.csv(node_all,"02-2 node file.csv")


height<- c(rep(1,m))
barplot(height,
          col=COLOR,
          names.arg=M,
        legend.text=M
        )
