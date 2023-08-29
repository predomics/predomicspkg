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

# # # =====================================================================================
# # # 29/08/2023: create the pre-computed results for testing new visualization functions
# # # =====================================================================================
# # # (1) cirrhosis stage 1
# # Prepare the data
# X <- cir_train$X; y <- cir_train$y # set global variables
# X <- X[rowSums(X)!=0,]; dim(X) # filter out variables with only zero values
# X <- filterNoSignal(X = X, side = 1, threshold = "auto", verbose = FALSE); dim(X) 
# # Run prediction model; Terga2 + terinter
# clf.terga2.terinter <- terga2(nCores = 1, 
#                               seed = 1, 
#                               plot = TRUE, language = "terinter"
# )
# res_clf_terinter <- fit(X = X, y = y, clf = clf.terga2.terinter, cross.validate = TRUE, nfolds = 10)
# # Run prediction model; Terga2 + bininter
# clf.terga2.bininter <- terga2(nCores = 1, 
#                               seed = 1, 
#                               plot = TRUE, language = "bininter"
# )
# res_clf_bininter <- fit(X = X, y = y, clf = clf.terga2.bininter, cross.validate = TRUE, nfolds = 10)
# # Save results in list
# list.results_cir_train <- list()
# list.results_cir_train[["terinter"]] <- res_clf_terinter
# list.results_cir_train[["bininter"]] <- res_clf_bininter
# save(list.results_cir_train, file="data/list.results_cir_train.rda")
# 
# # # (2) t2d
# # Prepare data
# X <- t2d$X; y <- t2d$y # set global variables
# X <- X[rowSums(X)!=0,]; dim(X) # filter out variables with only zero values
# X <- filterNoSignal(X = X, side = 1, threshold = "auto", verbose = FALSE); dim(X) 
# # Run prediction model; Terga2 + terinter
# clf.terga2.terinter <- terga2(nCores = 1, 
#                               seed = 1, 
#                               plot = TRUE, language = "terinter"
# )
# res_clf_terinter <- fit(X = X, y = y, clf = clf.terga2.terinter, cross.validate = TRUE, nfolds = 10)
# # Run prediction model; Terga2 + bininter
# clf.terga2.bininter <- terga2(nCores = 1, 
#                               seed = 1, 
#                               plot = TRUE, language = "bininter"
# )
# res_clf_bininter <- fit(X = X, y = y, clf = clf.terga2.bininter, cross.validate = TRUE, nfolds = 10)
# # Save results in list
# list.results.t2d <- list()
# list.results.t2d[["terinter"]] <- res_clf_terinter
# list.results.t2d[["bininter"]] <- res_clf_bininter
# save(list.results.t2d, file="data/list.results.t2d.rda")
# 
# # # (3) ibd
# # Prepare data
# X <- ibd$X; y <- ibd$y # set global variables
# X <- X[rowSums(X)!=0,]; dim(X) # filter out variables with only zero values
# X <- filterNoSignal(X = X, side = 1, threshold = "auto", verbose = FALSE); dim(X) 
# # Run prediction model; Terga2 + terinter
# clf.terga2.terinter <- terga2(nCores = 1, 
#                               seed = 1, 
#                               plot = TRUE, language = "terinter"
# )
# res_clf_terinter <- fit(X = X, y = y, clf = clf.terga2.terinter, cross.validate = TRUE, nfolds = 10, parallelize.folds = TRUE)
# # run prediction model; Terga2 + bininter
# clf.terga2.bininter <- terga2(nCores = 1, 
#                               seed = 1, 
#                               plot = TRUE, language = "bininter"
# )
# res_clf_bininter <- fit(X = X, y = y, clf = clf.terga2.bininter, cross.validate = TRUE, nfolds = 10)
# # Save results in list
# list.results.ibd <- list()
# list.results.ibd[["terinter"]] <- res_clf_terinter
# list.results.ibd[["bininter"]] <- res_clf_bininter
# save(list.results.ibd, file="data/list.results.ibd.rda")
# 
# # Clean workspace
# rm(list=ls()); gc()
