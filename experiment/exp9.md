# 实验九 制作人类基因剪接位点的GT-AG序列标识

一、实验目的

对于真核基因来说，几乎所有内含子DNA序列5'端起始的两个核苷酸总是5'-GT-3'，而3'端的最后两个核苷酸始终是5'-AG-3'。由于这两个碱基序列的高度保守性和广泛存在性，将其称为GT-AG法则，即5'-GT---AG-3'。利用序列标识（sequence logo），可以把GT-AG法则通过图形表示出来。而WebLogo便是一个灵活方便的制作序列标识的工具。针对人类基因组中22号染色体（chr22）上的所有基因，利用集成在Galaxy中的WebLogo和网络版的WebLogo，制作剪接位点的GT-AG序列标识。

1. 学习和掌握剪接位点的GT-AG法则。
2. 学习和掌握Galaxy的基本使用方法。
3. 学习和掌握WebLogo的使用方法。
4. 学习和掌握序列标识的含义。

二、实验内容

1. 打开Galaxy网站。通过搜索引擎搜索“Galaxy UCSC”，或者直接在浏览器的地址栏中输入网址，打开[https://main.g2.bx.psu.edu/](https://main.g2.bx.psu.edu/)即可。
2. 获取chr22上基因的内含子数据。打开Get Data工具集中的UCSC Main工具，调整参数提取人类hg19基因组中22号染色体（chr22）上的所有内含子信息，以BED格式进行存储。
3. 提取内含子上的剪接位点信息。
	1. 提取供体位点的信息。打开Operate on Genomic Intervals工具集中的Get flanks工具，利用刚刚获取的内含子数据，调整Region为Around Start，调整Location of the flanking region/s为Upstream，设置Offset为17，设置Length of flanking region(s)为32，提取出包含供体位点2bp在内、同时上下游各延伸15bp的坐标。
	2. 提取受体位点的信息。打开Operate on Genomic Intervals工具集中的Get flanks工具，利用刚刚获取的内含子数据，调整Region为Around End，调整Location of the flanking region/s为Downstream，设置Offset为-17，设置Length of flanking region(s)为32，提取出包含受体位点2bp在内、同时上下游各延伸15bp的坐标。
4. 获取剪接位点附近的序列。
	1. 提取供体位点的序列。打开Fetch Sequences工具集中的Extract Genomic DNA工具，利用供体位点的坐标信息，提取出以FASTA格式保存的供体位点的序列。
	2. 提取受体位点的序列。打开Fetch Sequences工具集中的Extract Genomic DNA工具，利用受体位点的坐标信息，提取出以FASTA格式保存的受体位点的序列。
5. 对序列进行多序列比对。此处提取的序列已经是根据坐标比对好的序列，没有必要再单独进行多序列比对了。
6. 制作序列标识。
	1. 制作供体位点的序列标识。打开Motif Tools中的Sequence Logo工具，利用供体位点的序列，调整参数，制作供体位点的序列标识。
	2. 制作受体位点的序列标识。打开Motif Tools中的Sequence Logo工具，利用受体位点的序列，调整参数，制作受体位点的序列标识。
7. WebLogo的使用。
	1. 下载剪接位点的序列。从Galaxy中下载提取出的供体位点和受体位点的序列，保存到本地计算机上。
	2. 打开WebLogo3网站。通过搜索引擎搜索“WebLogo”，或者直接在浏览器的地址栏中输入地址，打开[http://weblogo.threeplusone.com](http://weblogo.threeplusone.com)即可。
	3. 制作剪接位点的序列标识。在creat页面，上传供体位点或受体位点的序列，适当修改参数，就可制作出精美的序列标识。
8. 尝试对其他染色体或全基因组上的基因的剪接位点制作序列标识，进一步熟悉在Galaxy和WebLogo的使用方法。

