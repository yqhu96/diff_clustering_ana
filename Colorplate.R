library(RColorBrewer)
display.brewer.all() #显示所有调色板
display.brewer.all(type = "qual")
colorRampPalette(brewer.pal(8,"Accent"))(n)
