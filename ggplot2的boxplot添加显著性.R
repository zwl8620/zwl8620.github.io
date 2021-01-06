#---------------------------------------->
# 初次理解有显著性的箱线图
# 关键词：ggplot的boxplot添加显著性 | Add P-values and Significance Levels to ggplots | 方差分析
#---------------------------------------->











#---------------------------------------->
# copy教程一：https://www.cnblogs.com/leezx/p/10440050.html
#---------------------------------------->

# 多组比较，挑选感兴趣的显示显著性。
library(ggpubr)

my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )

options(repr.plot.width=4, repr.plot.height=4)

ggplot(ToothGrowth, aes(x=as.character(dose), y=len, fill=dose)) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA) +
  geom_jitter(width = 0.3, size=0.01) +# , aes(color=supp) +
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50, label.x = 1.5)     # Add global p-value







#---------------------------------------->
# copy教程二：https://www.r-bloggers.com/2017/06/add-p-values-and-significance-levels-to-ggplots/
#---------------------------------------->

# Install and load required R packages

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")

library(ggpubr)




data("ToothGrowth")
head(ToothGrowth)




### Compare two independent groups
p <- ggboxplot(ToothGrowth, x = "supp", y = "len",
               color = "supp", palette = "jco",
               add = "jitter")

#  Add p-value
p + stat_compare_means()
# Change method
p + stat_compare_means(method = "t.test")


p + stat_compare_means( aes(label = ..p.signif..), 
                        label.x = 1.5, label.y = 40)
# 调整Y轴刻度范围
# p + stat_compare_means( label = "p.signif", label.x = 1.5, label.y = 40)








### Compare two paired samples
ggpaired(ToothGrowth, x = "supp", y = "len",
         color = "supp", line.color = "gray", line.size = 0.4,
         palette = "jco")+
  stat_compare_means(paired = TRUE)









### Compare more than two groups
#  Global test
compare_means(len ~ dose,  data = ToothGrowth, method = "anova")


## Plot with global p-value:
# Default method = "kruskal.test" for multiple groups
ggboxplot(ToothGrowth, x = "dose", y = "len",
          color = "dose", palette = "jco")+
  stat_compare_means()


# Change method to anova
ggboxplot(ToothGrowth, x = "dose", y = "len",
          color = "dose", palette = "jco")+
  stat_compare_means(method = "anova")






# Pairwise comparisons. 
# 如果分组变量包含两个以上的级别，则将自动执行成对测试。默认的方法是wilcox。测试。您可以将其更改为t.test。

# Add p-values and significance levels to ggplots

## Perorm pairwise comparisons
compare_means(len ~ dose,  data = ToothGrowth)

# Visualize: Specify the comparisons you want
my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )

ggboxplot(ToothGrowth, x = "dose", y = "len",
          color = "dose", palette = "jco")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value


# 如果要指定bar的精确y位置，请使用参数label.y，简单来说就是调整Y轴刻度范围
ggboxplot(ToothGrowth, x = "dose", y = "len",
          color = "dose", palette = "jco")+ 
  stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+
  stat_compare_means(label.y = 45)







## Multiple pairwise tests against a reference group:

# Pairwise comparison against reference
compare_means(len ~ dose,  data = ToothGrowth, ref.group = "0.5",
              method = "t.test")


# Visualize
ggboxplot(ToothGrowth, x = "dose", y = "len",
          color = "dose", palette = "jco")+
  stat_compare_means(method = "anova", label.y = 40)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "0.5")                    # Pairwise comparison against reference






# Multiple pairwise tests against all (base-mean):

# Comparison of each group against base-mean
compare_means(len ~ dose,  data = ToothGrowth, ref.group = ".all.",
              method = "t.test")

# Visualize
ggboxplot(ToothGrowth, x = "dose", y = "len",
          color = "dose", palette = "jco")+
  stat_compare_means(method = "anova", label.y = 40)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")                  # Pairwise comparison against all



#------------------------------------->这种分析组别的举例
# 此处使用survminer软件包中的骨髓瘤数据集说明了一种典型的情况，即与“所有”的成对比较可能有用。
# 
# 我们将根据患者的分子类型绘制DEPDC1基因的表达谱。我们想知道组之间是否有任何区别。如果是，区别在哪里？
# 
# 要回答此问题，您可以在所有7个组之间进行成对比较。这将导致所有可能组合之间的大量比较。如果您有很多小组（如此处），可能很难解释。
# 
# 另一个简单的解决方案是将七个组中的每个与“全部”（即基本均值）进行比较。当测试显着时，则可以得出结论：与所有组相比，xxx组中DEPDC1显着过表达或过表达。
# Load myeloma data from survminer package
if(!require(survminer)) install.packages("survminer")

data("myeloma", package = "survminer")

# Perform the test
compare_means(DEPDC1 ~ molecular_group,  data = myeloma,
              ref.group = ".all.", method = "t.test")

# Visualize the expression profile
ggboxplot(myeloma, x = "molecular_group", y = "DEPDC1", color = "molecular_group", 
          add = "jitter", legend = "none") +
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(myeloma$DEPDC1), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 1600)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")                      # Pairwise comparison against all







## Multiple grouping variables

## Two independent sample comparisons after grouping the data by another variable:

# Box plot facetted by "dose"
p <- ggboxplot(ToothGrowth, x = "supp", y = "len",
               color = "supp", palette = "jco",
               add = "jitter",
               facet.by = "dose", short.panel.labs = FALSE)
# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format")
# Add p-values and significance levels to ggplots

# Or use significance symbol as label
p + stat_compare_means(label =  "p.signif", label.x = 1.5)
# 想要隐藏"ns"标签，使用参数hide.ns = TRUE



# Visualize (2/2). Create one single panel with all box plots. Plot y = “len” by x = “dose” and color by “supp”:
p <- ggboxplot(ToothGrowth, x = "dose", y = "len",
               color = "supp", palette = "jco",
               add = "jitter")
p + stat_compare_means(aes(group = supp))
# Show only p-value
p + stat_compare_means(aes(group = supp), label = "p.format")
# Use significance symbol as label
p + stat_compare_means(aes(group = supp), label = "p.signif")


## Paired sample comparisons after grouping the data by another variable:

# Perform the test:
compare_means(len ~ supp, data = ToothGrowth, 
              group.by = "dose", paired = TRUE)


# Visualize. Create a multi-panel box plots facetted by group (here, “dose”):
# Box plot facetted by "dose"
p <- ggpaired(ToothGrowth, x = "supp", y = "len",
              color = "supp", palette = "jco", 
              line.color = "gray", line.size = 0.4,
              facet.by = "dose", short.panel.labs = FALSE)
# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", paired = TRUE)





















#-------------------------------->Other plot types
# Bar and line plots (one grouping variable):
# Bar plot of mean +/-se
ggbarplot(ToothGrowth, x = "dose", y = "len", add = "mean_se")+
  stat_compare_means() +                                         # Global p-value
  stat_compare_means(ref.group = "0.5", label = "p.signif",
                     label.y = c(22, 29))                   # compare to ref.group

# Line plot of mean +/-se
ggline(ToothGrowth, x = "dose", y = "len", add = "mean_se")+
  stat_compare_means() +                                         # Global p-value
  stat_compare_means(ref.group = "0.5", label = "p.signif",
                     label.y = c(22, 29))     


# Bar and line plots (two grouping variables):
ggbarplot(ToothGrowth, x = "dose", y = "len", add = "mean_se",
          color = "supp", palette = "jco", 
          position = position_dodge(0.8))+
  stat_compare_means(aes(group = supp), label = "p.signif", label.y = 29)

ggline(ToothGrowth, x = "dose", y = "len", add = "mean_se",
       color = "supp", palette = "jco")+
  stat_compare_means(aes(group = supp), label = "p.signif", 
                     label.y = c(16, 25, 29))


