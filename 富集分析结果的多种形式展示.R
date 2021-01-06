#---------------------------------------->
# Type: 纯搬运代码
# AIM: Y叔的富集分析结果的热图-GSEA图等多种形式展示.R，参考：https://www.jianshu.com/p/c45cc2e3890f
#Y叔的一本神奇的书：https://yulab-smu.github.io/clusterProfiler-book/chapter12.html#bar-plot
# Time: 2021.01.06
# Author: lianzhiwei
#---------------------------------------->
# rm(list = ls())



# 设置工作路径
setwd("D:/language_learn/R_universal/2021/2021.01--纯代码搬运Y叔的富集分析教程")


# 1.安装和加载包

if(!require(clusterProfiler)) BiocManager::install("clusterProfiler")
if(!require(dplyr)) install.packages("dplyr")
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following object is masked from 'package:AnnotationDbi':
#> 
#>     select
#> The following objects are masked from 'package:IRanges':
#> 
#>     collapse, desc, intersect, setdiff, slice, union
#> The following objects are masked from 'package:S4Vectors':
#> 
#>     first, intersect, rename, setdiff, setequal, union
#> The following object is masked from 'package:Biobase':
#> 
#>     combine
#> The following objects are masked from 'package:BiocGenerics':
#> 
#>     combine, intersect, setdiff, union
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
if(!require(ggplot2)) install.packages("ggplot2")
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(cowplot)







# 2.示例数据
# 使用limma差异分析的结果，在生信星球公众号回复“富集输入”即可获得

# 提取差异表达基因的logFC，设置logFC阈值为1.5，即差异表达倍数为
load("step4output.Rdata")
geneList=deg$logFC
names(geneList)=deg$ENTREZID
geneList=sort(geneList,decreasing = T)
geneList2 = geneList[abs(geneList) > 1.5]
gene <- names(geneList)[abs(geneList) > 1.5]







# 3.GO富集分析
# (1)富集分析
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = "org.Hs.eg.db",
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)



# (2)富集分析可视化
#1)条带图
p1=barplot(ego,showCategory=20)

#2)点图
p2=dotplot(ego)



#(3)展示top5通路的共同基因
#Gene-Concept Network
p3=cnetplot(ego, categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)
p4=cnetplot(ego, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
gg <- ggdraw() +     
  draw_plot(p1, x=0, y=0.5, width=0.5, height=0.5) +  
  draw_plot(p2, 0, 0, 0.5, 0.5) +   
  draw_plot(p3, 0.5, 0, 0.5, 0.5) + 
  draw_plot(p4,0.5, 0.5, 0.5, 0.5)
gg


#Enrichment Map
emapplot(ego)



#(4)展示通路关系
goplot(ego)



#(5)Heatmap-like functional classification
heatplot(ego,foldChange = geneList)
# 这张热图调整画布大小才好看
# pdf("heatplot.pdf",width = 14,height = 5)
# heatplot(ego,foldChange = geneList)
# dev.off()








# 4.gsea-GO分析
ego2 <- gseGO(geneList     = geneList2,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

# 4.1 ridgeplot
ridgeplot(ego2)



# gesa图
# 4.2 一图只展示一条通路
gseaplot2(ego2, geneSetID = 1, title = ego2$Description[1])


# 4.3 多条通路同时展示
gseaplot2(ego2, geneSetID = 1:3)



#上下三联
pp <- lapply(1:3, function(i) {
  anno <- ego2[i, c("NES", "pvalue", "p.adjust")]
  lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
  
  gsearank(ego2, i, ego2[i, 2]) + xlab(NULL) +ylab(NULL) +
    annotate("text", 0, ego2[i, "enrichmentScore"] * .9, label = lab, hjust=0, vjust=0)
})
plot_grid(plotlist=pp, ncol=1)

# ggsave defaults to saving the last plot that you displayed
# 保存前记得调比例
ggsave("try.pdf",width = 15,height = 18)
