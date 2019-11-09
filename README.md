# 2019DCIC_fengchang
2019DCI-数字中国创新大赛-海上风场SCADA数据缺失智能修复-三蹦子-第4名

金风科技——团队三蹦子代码简述

# 运行环境说明
本代码运行环境为windows 7和10,
请确保至少有16G运行内存和至少20G硬盘空间。
本代码运行基于R 3.5, 建议用Rstudio作为IDE。
因为需要编译C++代码，所以还需要安装Rtools 3.5。
我们用到以下R packages，请确保更新至当前最新版。
tidyverse  
magrittr  
lubridate  
gridExtra  
imputeTS  
pracma  
zeallot  
Rcpp  
mgcv  
myCfun  

其中最后这个myCfun的R package是我们自己写的package, 用来存一些C++写的函数。
R可以调用它们来做weighted K neareast neigbhour和合并不同机组数据。
你可以在根目录下在找到源文件myCfun_0.1.0.tar.gz。

为了方便安装这些pacakge, 我们在根目录下提供了代码step0_installPackage.R
用以安装pacakges.


# 代码执行顺序
需要先把所给数据001到033这33个文件夹放在文件夹dataset下，
然后把template_submit_result.csv放在SCADA_sanbengzi根目录下。

实际上应用的时候是依次执行：
step1_fillNArow2rds.R -- 补充整行缺失的观察值，并输出新数据到文件夹addNA。  
step2_VIS_SCADA.R -- (辅助分析，可跳过)对每个机组每个变量做时间序列图以及可视化ACF和做差分后的ACF，输出到文件夹vis。我们在这个目录下提供了一个例子。  
step3_singleVarImpute.R -- 单变量（单条序列）做缺失值填补，输出预测结果到output_step3 并输出log到log_step3。 
step4a_clusterVar.R -- 聚类分析，利用变量间的相关性做Hierarchical clustering，将68个变量分成28个cluster。  
step4b_genClusterData.R -- 根据上一步得到的cluster, 把不同的机组在同一个cluster里的变量合并， 输出新的cluster-based data到文件夹clusterData。  
step5_clusterFill -- 这个是个文件夹， 里面包含代码用于执行多变量（多条序列,单cluster）做缺失值填补。输出预测结果到output_step5和log到log_step5。  
step6_merge.R -- 以单变量预测结果为主，与多变量预测的结果合并。哪些多变量值得合并主要看本地分数（log_step3和log_step5）和线上成绩是否提升，并输出结果sub_fullcluster.csv到文件夹submit。

其中做完step3_singleVarImpute.R的话，a榜有0.679+的分数(这里本地分数也有0.679+)。  
step5_clusterFill做完以下4个主要的cluster, 结合step3就有0.703+的a榜分数。  
(当用上其他变量后，我们虽然没去具体算本地增加的总的分数，不过感觉很容易overfit，特别是如果用上了其他机组的数据)  
1_5_38  
19_53_59_65_68  
4_27_34_42_43_46  
3_24_35_40_45_52  
之后结合其他其他cluster的话能达到我们目前最好的a榜 0.7051, b榜0.7056。

为了避免频繁设置代码工作路径，我们提供了allstep.R这个代码文件可以依次执行以上步骤。  
还有其他两个辅助代码文件  
FE_utilityMar11.R -- 用于单/多序列缺失值填补  
localApproximateV2fix.R -- 用于多序列缺失值填补  


# 运行所需时间
由于我们的方法是对的单机组单变量，或者多机组多变量（单cluster)执行，
所以实际运行的时候可以很简单地进行并行计算。
这里我们提供的代码为了方便调试，都不是并行的写法。
不过要实施并行也很简单，特别是针对最耗时的step5_clusterFill。
只需要开多几个Rstudio，分别运行不同cluster填补缺失值的代码，
或者不同机组的缺失值填补分开来运行。

以下给出CPU i5-2500k单核的大概运行时间，  
step3_singleVarImpute.R -- 15 hours  
step4b_genClusterData.R -- 需要写数据进去硬盘， 大概5小时。  
step5_clusterFill -- 所有cluster跑完估计至少7天， 只跑那几个提升多的cluster估计3天。  

