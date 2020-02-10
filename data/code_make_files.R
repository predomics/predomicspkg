# # =====================================================================================
# # create the databases
# # =====================================================================================
# 
# # (1) cirrhosis stage 1
# load("/data/projects/predomics_testing/data/segata_2017/data.bugs/qinn_bug_stage1.rda")
# cir_train <- list()
# cir_train$X <- as.data.frame(data.qinn.bug.stage1.freq.species); dim(cir_train$X) # 1045  181
# cir_train$y <- data.qinn.bug.y
# save(cir_train, file="cir_train.rda", compress = TRUE, compression_level = 9)
# rm(list=ls()); gc()
# 
# # (2) cirrhosis stage 2
# load("/data/projects/predomics_testing/data/segata_2017/data.bugs/qinn_bug_stage2.rda")
# cir_test <- list()
# cir_test$X <- as.data.frame(data.qinn.bug.stage2.freq.species); dim(cir_test$X) # 1045  56
# cir_test$y <- data.qinn.bug.y
# save(cir_test, file="cir_test.rda", compress = TRUE, compression_level = 9)
# rm(list=ls()); gc()
# 
# # (3) ibd
# load("/data/projects/predomics_testing/data/segata_2017/data.bugs/nielsen_bug.rda")
# ibd <- list()
# ibd$X <- as.data.frame(data.nielsen.bug.freq.species); dim(ibd$X) # 1045  396
# ibd$y <- data.nielsen.bug.y
# save(ibd, file="ibd.rda", compress = TRUE, compression_level = 9)
# rm(list=ls()); gc()
# 
# # (4) obesity
# load("/data/projects/predomics_testing/data/segata_2017/data.bugs/lechat_bug.rda")
# obesity <- list()
# obesity$X <- as.data.frame(data.lechat.bug.freq.species); dim(obesity$X) # 1045  292
# obesity$y <- data.lechat.bug.y
# save(obesity, file="obesity.rda", compress = TRUE, compression_level = 9)
# rm(list=ls()); gc()
# 
# # (5) t2d
# load("/data/projects/predomics_testing/data/segata_2017/data.bugs/qinj_bug.rda")
# t2d <- list()
# t2d$X <- as.data.frame(data.qinj.bug.freq.species); dim(t2d$X) # 1045  344
# t2d$y <- data.qinj.bug.y
# save(t2d, file="t2d.rda", compress = TRUE, compression_level = 9)
# rm(list=ls()); gc()
# 
# # (6) t2dw
# load("/data/projects/predomics_testing/data/segata_2017/data.bugs/karlsson_bug.rda")
# t2dw <- list()
# t2dw$X <- as.data.frame(data.karlsson.bug.freq.species); dim(t2dw$X) # 1045  344
# t2dw$y <- data.karlsson.bug.y
# save(t2dw, file="t2dw.rda", compress = TRUE, compression_level = 9)
# rm(list=ls()); gc()
# 
# # (7) cirrhosis stage 1 counts
# load("/data/projects/predomics_testing/data/segata_2017/data.bugs/qinn_bug_stage1.rda")
# cir_train_count <- list()
# cir_train_count$X <- as.data.frame(data.qinn.bug.stage1.counts.species); dim(cir_train_count$X) # 1045  181
# cir_train_count$y <- data.qinn.bug.y
# save(cir_train_count, file="cir_train_count.rda", compress = TRUE, compression_level = 9)
# rm(list=ls()); gc()

