#https://blog.csdn.net/weixin_54000907/article/details/118503051
library(RColorBrewer)
display.brewer.all() #显示所有调色板
display.brewer.all(type = "qual")
colorRampPalette(brewer.pal(8,"Accent"))(n)
