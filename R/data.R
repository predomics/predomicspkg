#-----------------------------------------------------
# TEST DATASETS
#-----------------------------------------------------

#' @name cir_train
#' @title Cirhosis stage 1 (frequencies)
#' @docType data
#' @author Qin, Nan, Fengling Yang, Ang Li, Edi Prifti, Yanfei Chen, Li Shao, Jing Guo, et al “Alterations of the human gut microbiome in liver cirrhosis.” Nature 513, no. 7516 (July 23, 2014): 59–64 _.
#' @keywords liver cirrhosis, microbiome, species
#' @description This dataset consists of frequency abundance files as downloaded from http://waldronlab.io/curatedMetagenomicData/. 
#' This is a list containing two elements: (i) the X data matrix with 1045 species and 181 observations and (ii) patient class = -1 (n=98) and healthy controls (n=83)
NULL

#' @name cir_test
#' @title Cirhosis stage 2 (frequencies)
#' @docType data
#' @author Qin, Nan, Fengling Yang, Ang Li, Edi Prifti, Yanfei Chen, Li Shao, Jing Guo, et al “Alterations of the human gut microbiome in liver cirrhosis.” Nature 513, no. 7516 (July 23, 2014): 59–64.
#' @keywords liver cirrhosis, microbiome, species
#' @description This dataset consists of frequency abundance files as downloaded from http://waldronlab.io/curatedMetagenomicData/
#' This is a list containing two elements: (i) the X data matrix with 1045 species and 56 observations and (ii) patient class = -1 (n=25) and healthy controls (n=31)
NULL

#' @name ibd
#' @title Inflammatory Bowel Disease (frequencies) from the MetaHIT study
#' @docType data
#' @author Nielsen, H Bjørn, Mathieu Almeida, Agnieszka Sierakowska Juncker, Simon Rasmussen, Junhua Li, Shinichi Sunagawa, Damian R Plichta, et al “Identification and assembly of genomes and genetic elements in complex metagenomic samples without using reference genomes.” Nature biotechnology (July 6, 2014): 1–11.
#' @keywords inflamatory bowel disease, microbiome, species
#' @description This dataset consists of frequency abundance files as downloaded from http://waldronlab.io/curatedMetagenomicData/
#' This is a list containing two elements: (i) the X data matrix with 1045 species and 396 observations and (ii) patient class = -1 (n=148) and healthy controls (n=248)
NULL

#' @name obesity
#' @title Obesity (frequencies) from the MetaHIT study
#' @docType data
#' @author Le Chatelier, Emmanuelle, Trine Nielsen, Junjie Qin, Edi Prifti, Falk Hildebrand, Gwen Falony, Mathieu Almeida, et al “Richness of human gut microbiome correlates with metabolic markers.” Nature 500, no. 7464 (April 9, 2014): 541–546.
#' @keywords obesity, microbiome, species
#' @description This dataset consists of frequency abundance files as downloaded from http://waldronlab.io/curatedMetagenomicData/
#' This is a list containing two elements: (i) the X data matrix with 1045 species and 292 observations and (ii) patient class = -1 (n=167) and healthy controls (n=96).
#' Caution, this dataset has also a class 0 with overweight patients, which needs to be omited from both X and y
NULL

#' @name t2d
#' @title Type 2 diabetes (frequencies) BGI
#' @docType data
#' @author Qin, Junjie, Yingrui Li, Zhiming Cai, Shenghui Li, Jianfeng Zhu, Fan Zhang, Suisha Liang, et al “A metagenome-wide association study of gut microbiota in type 2 diabetes.” Nature (September 26, 2012).
#' @keywords type 2 diabetes, microbiome, species
#' @description This dataset consists of frequency abundance files as downloaded from http://waldronlab.io/curatedMetagenomicData/
#' This is a list containing two elements: (i) the X data matrix with 1045 species and 344 observations and (ii) patient class = -1 (n=170) and healthy controls (n=174)
NULL

#' @name t2dw
#' @title Type 2 diabetes (frequencies) Women Sweden
#' @docType data
#' @author Karlsson, Fredrik H, Valentina Tremaroli, Intawat Nookaew, Göran Bergström, Carl Johan Behre, Björn Fagerberg, Jens Nielsen, and Fredrik Bäckhed. “Gut metagenome in European women with normal, impaired and diabetic glucose control.” Nature (May 29, 2013): 1–7.
#' @keywords type 2 diabetes, microbiome, species
#' @description This dataset consists of frequency abundance files as downloaded from http://waldronlab.io/curatedMetagenomicData/
#' This is a list containing two elements: (i) the X data matrix with 1045 species and 145 observations and (ii) patient class = -1 (n=53) and healthy controls (n=43)
#' Caution, this dataset has also a class 0 with IG patients, which needs to be omited from both X and y
NULL

