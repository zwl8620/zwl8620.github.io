#--------------------------------------------->
# AIM：ChIPseeker注释单样本和多样本的bed文件并以饼图可视化注释结果的代码：
# Author: lianzhiwei
# Time: 2021.01.06
#--------------------------------------------->


#--------------------------------------------->
### 加载相关R包
### 加载物种基因注释文件，注意要和输入文件的物种相对应
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene

library(ChIPseeker)

# 设置工作路径
setwd("path/to/work/")
#--------------------------------------------->



#--------------------------------------------->
### 1.单个样本的注释
# 1.1 指定样本的bed文件路径，注意ChIPseeker不需要将输入文件读入为变量
peak1_file <-  "input_peak1.bed"



# 1.2 annotatePeak函数
peak <-annotatePeak(peak1_file, tssRegion=c(-200,2000), TxDb=txdb, annoDb="org.Hs.eg.db")
# 参数说明：
# tssRegion=c(-200, 2000)：定义启动子的TSS上下游区域为TSS上游2kb，下游0.2kb
# TxDb=txdb：指定注释物种的数据库



# 1.3 plotAnnoPie函数以饼图可视化注释结果，有时候默认的画布不够大会导致图例被部分遮住，调整画布大小即可
# pdf("plotAnnoPie.pdf", width = 8, height = 8)
plotAnnoPie(peak)
# dev.off()
#--------------------------------------------->








#--------------------------------------------->
### 2. 多个样本同时注释
# 2.1 指定样本的bed文件路径，注意ChIPseeker不需要将输入文件读入为变量
peak1_file <-  "input_peak1.bed"
peak2_file <-  "input_peak2.bed"

# 将多个样本的文件路径整合为列表，方便函数使用
peak_file_List <-  list(peak1 = peak1_file, peak2 = peak2_file)



# 2.2 lapply循环调用annotatePeak函数进行注释
peakAnnoList <- lapply(peak_file_List, annotatePeak,  TxDb=txdb, tssRegion=c(-200, 2000),  verbose=FALSE, annoDb="org.Hs.eg.db")



# 2.3 plotAnnoPie函数以饼图可视化注释结果，有时候默认的画布不够大会导致图例被部分遮住，调整画布大小即可
# pdf("plotAnnoPie_multi_Sample.pdf")   
plotAnnoPie(peakAnnoList[[1]])
plotAnnoPie(peakAnnoList[[2]])
# dev.off()
#--------------------------------------------->








#--------------------------------------------->
### 3.ChIPseeker其他功能的拓展分析
# 3.1 plotAnnoBar函数以比例图可视化注释结果
# pdf("plotAnnoBar_multi_Sample.pdf")
plotAnnoBar(peakAnnoList)
# dev.off()



# 3.2 plotDistToTSS函数
genes= lapply(peakAnnoList, function(i)  as.data.frame(i)$geneId)
vennplot(genes)
plotDistToTSS(peakAnnoList)
#--------------------------------------------->


# ChIPseeker包还有更多的绘图和分析函数，尝试去了解