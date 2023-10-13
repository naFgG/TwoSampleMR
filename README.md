# 两样本MR
## 孟德尔随机化概念
孟德尔随机化在统计学上的本质实际是利用工具变量（Instrumental Variables, IV）来研究因果性，这一方法常用在经济学研究中。

工具变量简单来说就是，一个与`X`相关，但与被忽略的混淆因素以及`Y`不相关的变量。在经济学研究中工具变量可以是政策改革，自然灾害等等，而在遗传学中，这个变量就是基因。

如果一个基因变异`Z` 是某个暴露因素`X`的因果变量，并且对结果`Y`没有直接因果关系，那么这个基因变异`Z`与结果`Y`的关联，只能通过`X`对`Y`的因果关系而被观察到（X -> Y）。

## 方法概览
使用`TwoSampleMR`R包运行的两样本MR简单流程   
![alt text](https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "flows")

+ 主要步骤
1. 读取暴露因素GWAS数据

2. 选取合适的工具变量【如有必要，需要进行clumping。对于一个标准的两样本MR分析来说，我们需要确保工具变量之间是互相独立的，即不存在显著的连锁不平衡（LD），在读取完数据后应当对其进行LD Clumping】

3. 读取结果GWAS数据，提取上述的工具变量的SNP

4. 对暴露因素与结果的GWAS数据进行预处理，使其格式统一

5. 采用选中的MR方法进行分析

6. 将结果可视化：  
   1. 散点图  
   2. 森林图  
   3. 漏斗图  

## 数据格式
对于暴露因素的GWAS数据，TwoSampleMR需要一个工具变量的data frame，每行对应一个SNP，至少需要4列，分别为：  
+ SNP – rs ID  
+ beta – The effect size. If the trait is binary then log(OR) should be used  
+ se – The standard error of the effect size  
+ effect_allele – The allele of the SNP which has the effect marked in beta（effect_allele指基因组上与特定表型相关的等位基因，对表型和疾病有贡献的等位基因）
  
<br>

其他可能有助于MR预处理或分析的列包括：  
+ other_allele – The non-effect allele
+ eaf – The effect allele frequency
+ Phenotype – The name of the phenotype for which the SNP has an effect
  
<br>

可以提供额外的信息（可选，对于分析非必须）：  
+ chr – Physical position of variant (chromosome)
+ position – Physical position of variant (position)
+ samplesize – Sample size for estimating the effect size
+ ncase – Number of cases
+ ncontrol – Number of controls
+ pval – The P-value for the SNP’s association with the exposure
+ units – The units in which the effects are presented
+ gene – The gene or other annotation for the the SNP  

<br>

*参考资料*：  
[1] https://gwaslab.org/2021/06/24/mr/  
[2] https://mrcieu.github.io/TwoSampleMR/articles/introduction.html
