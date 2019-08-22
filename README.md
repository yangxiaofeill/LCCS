# LCCS
Detect LCCS (Local Clustering of Cytosine Sites) from WGBS/RRBS cytosine position data

#################################################################################### 

This program is used to detect Local Clustering of Cytosine Sites (LCCS) of 
single nucletide resolution methylation sequencing data 
Usage: `Rscript --vanilla LCCS.R minCNum maxSize posFile` 
Parameters: 
minCNum, the minimal Cytosine Number in a LCCS, default value is 3 
maxSize, the maxmize size of a LCCS, default value is 500bp 
posFile, the file includes positions of all cytosine, and at lest two column, like 
chr1 1 
chr1 5 
... 

#####################################################################################
Reference: doi: 10.1093/hmg/ddv172 copyright @ Xiaofei Yang, xfyang@xjtu.edu.cn
