#' @title predict
#'
#' @description Calculates predicted values from a fitted \code{Hmsc} model.
#'
#' @param object a fitted \code{Hmsc} model object
#' @param post a list of posterior samples of the HMSC model. By default uses all samples from the pooled
#' posterior of the hM object.
#' @param XData a dataframe specifying the unpreprocessed covariates for the predictions to be made.
#' Works only if the \code{XFormula} argument was specified in the \code{Hmsc()} model constructor call.
#' Requirements are similar to those in the \code{Hmsc} model constructor.
#' @param X a matrix specifying the covariates for the predictions to be made.
#' Only one of XData and X arguments may be provided.
#' @param XRRRData a dataframe of covariates for reduced-rank regression
#' @param XRRR a matrix of covariates for reduced-rank regression
#' @param studyDesign a matrix, specifying the structure of the study design for the prediction.
#' Requirements are similar to those of the \code{Hmsc} constructor. By default this argument is
#' assigned the study design of the training data in the fitted Hmsc model.
#' @param ranLevels a list of \code{HmscRandomLevel} objects, further specifying the structure of
#' random levels. Requirements are similar to those of the \code{Hmsc} constructor.
#' Each level must cover all units, specified in the correspondingly named column of \code{studyDesign}
#' argument. By default this argument is assigned the list of \code{HmscRandomLevel} objects
#' specified for fitting Hmsc model.
#'
#' @param Gradient an object returned by
#'     \code{\link{constructGradient}}. Providing \code{Gradient} is
#'     an alternative for providing \code{XData}, \code{studyDesign}
#'     and \code{ranLevels}. Cannot be used together with \code{Yc}.
#'
#' @param Yc a matrix of the outcomes that are assumed to be known for
#'     conditional predictions. Cannot be used together with
#'     \code{Gradient}.
#'
#' @param mcmcStep the number of extra mcmc steps used for updating the random effects
#' @param expected boolean flag indicating whether to return the location parameter of the observation
#' models or sample the values from those.
#' @param predictEtaMean boolean flag indicating whether to use the estimated mean values of posterior
#' predictive distribution for random effects corresponding for the new units.
#' @param predictEtaMeanField boolean flag indicating whether to use draws from the mean-field of the
#' posterior predictive distribution for random effects corresponding for the new units.
#'
#' @param nParallel Number of parallel processes. Parallel processing
#'     is only useful with new \code{Yc} data and extra
#'     \code{mcmcStep}.
#' @param useSocket (logical) Use socket clusters in parallel
#'     processing; these are the only alternative in Windows, but in
#'     other systems this should be usually set \code{FALSE} for
#'     forking.
#'
#' @param \dots other arguments passed to functions.
#'
#' @details In \code{mcmcStep,the number of extra mcmc steps used for updating the random effects
#' for the Eta parameters, starting from the samples of the fitted Hmsc model in order to
#' account for the conditional information provided in the Yc argument. The higher this number is,
#' the more the obtained updated samples are unaffected by the posterior estimates of latent factors
#' in the model fitted to the training data and more resembles the true conditional posterior. However,
#' the elapsed time for conditional prediction grows approximately linearly as this parameter increases.
#' The exact number for sufficient is problem-dependent and should be assessed by e.g. gradually
#' increasing this parameter till the stationarity of the produced predictions.}
#'
#' @return A list of length \code{length(post)}, each element of which contains a sample from the posterior
#' predictive distribution (given the sample of the Hmsc model parameters in the corresponding element of
#' the \code{post} argument)
#'
#' @seealso \code{\link{predictLatentFactor}}
#'
#'
#' @importFrom stats model.matrix rnorm pnorm rpois
#' @importFrom parallel detectCores mclapply makeCluster stopCluster
#'     clusterExport clusterEvalQ parLapply
#'
#' @export

predict.Hmsc = function(
      object,
      # post=poolMcmcChains(object$postList),
      Loff=NULL, XData=NULL, X=NULL, XRRRData=NULL, XRRR=NULL, # this has to be updated to cov-dependent associations
      # studyDesign=object$studyDesign, ranLevels=object$ranLevels,
      Gradient=NULL, Yc=NULL, mcmcStep=1, expected=FALSE,
      predictEtaMean=FALSE, predictEtaMeanField=FALSE,
      nParallel = 1,
      useSocket = TRUE, TempDir = "TEMP2Pred",
      PredDir = "Model_Predictions", Evaluate = FALSE, Evaluate_Name = NULL,
      EvalDir = "Eval",
      ModelName = "Train", RC = FALSE, Pred_PA = NULL, Pred_XY = NULL, ...) {


   fs::dir_create(c(TempDir, PredDir, EvalDir))

   if (inherits(object, "character")) {
      object <- IASDT.R::LoadAs(object)
   }
   post=poolMcmcChains(object$postList)
   studyDesign=object$studyDesign
   ranLevels=object$ranLevels


   ## check valid nParallel
   nParallel <- min(nParallel, detectCores())
   if (nParallel > 1) {
      if (.Platform$OS.type == "windows" && !useSocket) {
         useSocket <- TRUE
         message("setting useSocket=TRUE; the only choice in Windows")
      }
   }

   if(!is.null(Gradient)) {
      ## don't know what to do if there is also Yc, and spatial models
      ## will trigger an error in updateEta (github issue #135)
      if(!is.null(Yc)) stop("predict with arguments 'Yc' and 'Gradient' jointly is not implemented (yet)")
      XData=Gradient$XDataNew
      studyDesign=Gradient$studyDesignNew
      ranLevels=Gradient$rLNew
   }

   if(!is.null(XData) && !is.null(X)){
      stop("only one of XData and X arguments can be specified")
   }
   if(!is.null(XRRRData) && !is.null(XRRR)){
      stop("only one of XRRRData and XRRR arguments can be specified")
   }
   if(predictEtaMean==TRUE && predictEtaMeanField==TRUE)
      stop("predictEtaMean and predictEtaMeanField arguments cannot be TRUE simultaneously")

   if(!is.null(XData)){
      switch(
         class(XData)[1L],
         list={
            if(any(unlist(lapply(XData, is.na)))) {
               stop("NA values are not allowed in 'XData'")
            }
            xlev = lapply(Reduce(rbind,object$XData), levels)[unlist(lapply(Reduce(rbind,object$XData), is.factor))]
            X = lapply(XData, function(a) model.matrix(object$XFormula, a, xlev=xlev))
         },
         data.frame={
            if(any(is.na(XData))) stop("NA values are not allowed in 'XData'")
            xlev = lapply(object$XData, levels)[unlist(lapply(object$XData, is.factor))]
            X = model.matrix(object$XFormula, XData, xlev=xlev)
         }
      )
   } else{
      if(is.null(X))
         X = object$X
   }

   if(!is.null(XRRRData)){
      xlev = lapply(object$XRRRData, levels)[unlist(lapply(object$XRRRData, is.factor))]
      XRRR = model.matrix(object$XRRRFormula, XRRRData, xlev=xlev)
   } else{
      if(is.null(object$ncRRR)) object$ncRRR=0
      if(is.null(XRRR) && object$ncRRR>0)
         XRRR=object$XRRR
   }

   switch(
      class(X)[1L],
      list={ nyNew = nrow(X[[1]]) },
      matrix={ nyNew = nrow(X) }
   )

   if(!is.null(Yc)){
      if(ncol(Yc) != object$ns) stop("number of columns in Yc must be equal to ns")
      if(nrow(Yc) != nyNew) stop("number of rows in Yc and X must be equal")
   }
   if(!is.null(Loff)){
      if(ncol(Loff) != object$ns) stop("number of columns in Loff must be equal to ns")
      if(nrow(Loff) != nyNew) stop("number of rows in Loff and X must be equal")
   }
   if(!all(object$rLNames %in% colnames(studyDesign))){
      stop("dfPiNew does not contain all the necessary named columns")
   }
   if(!all(object$rLNames %in% names(ranLevels))){
      stop("rL does not contain all the necessary named levels")
   }

   if(!is.null(studyDesign)){
      dfPiNew = studyDesign[,object$rLNames,drop=FALSE]
   } else
      dfPiNew = matrix(NA,nyNew,0)
   rL = ranLevels[object$rLNames]

   if(!is.null(Yc)){
      ## object can have pre-computed data parameters, but not
      ## necessarily. These are needed only in updateEta(), but get it
      ## here anyway...
      if (is.null(object$rLPar)) {
         rLPar = computeDataParameters(object)$rLPar
      } else {
         rLPar = object$rLPar
      }
   } else{
      rLPar = NULL
   }

   predN = length(post)

   # free some memory
   object$postList <- object$YScaled <- object$X <- object$XScaled <- NULL

   predPostEta = vector("list", object$nr)
   PiNew = matrix(NA,nrow(dfPiNew),object$nr)

   for(r in seq_len(object$nr)) {
      # Do not use `predictLatentFactor` when predicting values for response
      # curves
      NewUnits <- all(levels(dfPiNew[,r]) == "new_unit")
      InfCoords <- tibble::as_tibble(ranLevels[[r]]$s, rownames = "ID") %>%
         tail() %>%
         dplyr::pull(ID) %>%
         magrittr::equals("new_unit") %>%
         all()

      if (NewUnits && InfCoords) {
         nLF <- length(post[[1]]$Alpha[[1]])

         predPostEta[[r]] <- replicate(
            n = predN,
            expr = structure(
               rep(0, nLF), dim = c(1L, nLF),
               dimnames = list("new_unit", NULL)),
            simplify = FALSE)
      } else {
         postEta = lapply(post, function(c) c$Eta[[r]])
         postAlpha = lapply(post, function(c) c$Alpha[[r]])

         if (r == object$nr){
            # free some memory
            post <- lapply(post, function(x) {
               x$Eta <- x$Psi <- x$V <- x$Delta <- x$Gamma <- x$rho <- NULL
               x
            })

            # Save post to file and load it later
            post_file <- file.path(TempDir, paste0(ModelName, "_post.qs"))
            qs::qsave(post, file = post_file, preset = "fast")
            rm(post)
            invisible(gc())
         }

         # Save postEta to file and load it from predictLatentFactor
         # This helps to avoid the unnecessary copying of the postEta object to all cores
         postEta_file <- file.path(TempDir, paste0(ModelName, "_postEta.qs"))
         qs::qsave(postEta, file = postEta_file, preset = "fast")
         rm(postEta)
         invisible(gc())


         predPostEta[[r]] = predictLatentFactor(
            unitsPred = levels(dfPiNew[,r]),
            units = levels(object$dfPi[,r]),
            postEta = postEta_file, postAlpha = postAlpha,
            rL = rL[[r]],
            predictMean=predictEtaMean,
            predictMeanField=predictEtaMeanField,
            nParallel = nParallel,
            TempDir = TempDir, ModelName = ModelName)

         fs::file_delete(postEta_file)
      }

      rowNames = rownames(predPostEta[[r]][[1]])
      # PiNew[,r] = sapply(dfPiNew[,r], function(s) which(rowNames==s))
      PiNew[,r] = fastmatch::fmatch(dfPiNew[, r], rowNames)
   }

   if(object$nr > 0){
      ppEta <- simplify2array(predPostEta)
   } else{
      ppEta <- matrix(list(),predN,0)
   }

   # free some memory
   try(rm(predPostEta), silent = TRUE)

   if (!exists("post")) {
      post <- qs::qread(post_file, nthreads = 5)
   }

   invisible(gc())

   # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

   if (RC) {
      lapply(seq_len(predN), function(pN, ...){
         get1prediction(object, X, XRRR, Yc, Loff, rL, rLPar, post[[pN]],
                        ppEta[pN,], PiNew, dfPiNew, nyNew, expected, mcmcStep)
      })
   } else {
      if (nParallel == 1) {  # non-Parallel
         pred <- lapply(seq_len(predN), function(pN, ...){
            get1prediction(
               object, X, XRRR, Yc, Loff, rL, rLPar, post[[pN]], ppEta[pN,],
               PiNew, dfPiNew, nyNew, expected, mcmcStep)
         })

      } else if (useSocket) { # socket cluster (Windows, mac, Linux)
         # Save ppEta / post as small chunks
         chunk_size <- 25
         ChunkIDs <- ceiling(seq_along(post) / chunk_size)
         Chunks <- purrr::map_chr(
            .x = seq_len(max(ChunkIDs)),
            .f = ~ {
               IDs <- which(ChunkIDs == .x)
               Ch <- list(ppEta = ppEta[IDs], post = post[IDs])
               ChunkFile <- file.path(
                  TempDir, paste0(ModelName, "_ppEta_ch", .x, ".qs"))
               qs::qsave(Ch, file = ChunkFile, preset = "fast")
               return(ChunkFile)
            })

         rm(ChunkIDs, post, ppEta)
         invisible(gc())




         seeds <- sample.int(.Machine$integer.max, predN)

         c1 <- snow::makeSOCKcluster(nParallel)
         on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
         snow::clusterExport(
            cl = c1,
            list = c(
               "object", "X", "XRRR", "Yc", "Loff", "rL", "rLPar", "PiNew",
               "dfPiNew", "nyNew", "expected", "mcmcStep", "seeds",
               "get1prediction", "chunk_size", "Chunks", "TempDir", "ModelName"),
            envir = environment())

         invisible(snow::clusterEvalQ(
            cl = c1,
            expr = {
               library(dplyr)
               library(Rcpp); library(RcppArmadillo); library(Matrix);
               library(purrr); library(Hmsc); library(float); library(qs)
            }))

         pred <- snow::parLapply(
            cl = c1,
            x = seq_len(length(Chunks)),
            fun = function(Chunk) {

               ChunkFile <- Chunks[Chunk]
               Ch <- qs::qread(ChunkFile)
               ppEta <- Ch$ppEta
               post <- Ch$post
               rm(Ch)
               Seed <- (Chunk - 1) * chunk_size
               Seed <- seq(Seed + 1, Seed + chunk_size)
               Seed <- seeds[Seed]

               PredChunk <- purrr::map(
                  .x = seq_len(chunk_size),
                  .f = function(pN) {
                     get1prediction(
                        object = object, X = X, XRRR = XRRR, Yc = Yc,
                        Loff = Loff, rL = rL, rLPar = rLPar, sam = post[[pN]],
                        predPostEta = ppEta[pN], PiNew = PiNew,
                        dfPiNew = dfPiNew, nyNew = nyNew, expected = expected,
                        mcmcStep = mcmcStep, seed = Seed[pN])
                  })

               ChunkSp <- purrr::map_dfr(
                  .x = seq_len(length(object$spNames)),
                  .f = function(Sp) {
                     SpD <- purrr::map(PredChunk, ~ .x[, Sp], ncol = 1) %>%
                        simplify2array() %>%
                        float::fl()
                     dimnames(SpD) <- NULL

                     ChunkSpFile <- file.path(
                        TempDir,
                        paste0("Pred_", ModelName, "_ch", Chunk, "_Sp", Sp, ".qs"))
                     qs::qsave(SpD, file = ChunkSpFile, preset = "fast")

                     return(
                        cbind.data.frame(
                           Chunk = Chunk, Sp = Sp, IAS_ID = object$spNames[Sp],
                           ChunkSpFile = ChunkSpFile))
                  })

               fs::file_delete(ChunkFile)
               rm(PredChunk)

               return(ChunkSp)
            })

         snow::stopCluster(c1)
         invisible(gc())

         pred <- tibble::tibble(dplyr::bind_rows(pred))

      } else { # fork (mac, Linux)
         seed <- sample.int(.Machine$integer.max, predN)
         pred <- mclapply(seq_len(predN), function(pN, ...)
            get1prediction(object, X, XRRR, Yc, Loff, rL, rLPar, post[[pN]],
                           ppEta[pN,], PiNew, dfPiNew, nyNew, expected,
                           mcmcStep, seed = seed[pN]),
            mc.cores=nParallel)
      }

      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      c1 <- snow::makeSOCKcluster(NCores)
      on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
      future::plan("future::cluster", workers = c1, gc = TRUE)
      on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)

      Eval_DT <- tibble::tibble(pred) %>%
         dplyr::select(-Chunk) %>%
         dplyr::group_nest(Sp, IAS_ID) %>%
         dplyr::mutate(data = purrr::map(data, unlist))

      Eval_DT <- future.apply::future_lapply(
         X = seq_len(nrow(Eval_DT)),
         FUN = function(ID) {
            Sp <- Eval_DT$Sp[[ID]]
            IAS_ID <- Eval_DT$IAS_ID[[ID]]
            data <- Eval_DT$data[[ID]]

            SpDT <- purrr::map(data, qs::qread) %>%
               do.call(cbind, .) %>%
               as.double()

            SpDT_Mean <- Rfast::rowmeans(SpDT)
            SpDT_SD <- Rfast::rowVars(SpDT, std = TRUE)
            SpDT_Cov <- SpDT_SD / SpDT_Mean
            rm(SpDT)

            if (is.null(Pred_XY)) {
               Pred_XY <- object$rL$sample$s
            }

            PredSummary <- tibble::tibble(
               as.data.frame(Pred_XY),
               Mean = SpDT_Mean, SD = SpDT_SD, Cov = SpDT_Cov) %>%
               stats::setNames(
                  c("x", "y", paste0(IAS_ID, "_mean"),
                  paste0(IAS_ID, "_sd"), paste0(IAS_ID, "_cov"))) %>%
               sf::st_as_sf(coords = c("x", "y"), crs = 3035, remove = FALSE)

            PredSummaryFile <- file.path(
               PredDir, paste0("Pred_", ModelName, "_Sp", Sp, ".qs"))
            qs::qsave(PredSummary, file = PredSummaryFile, preset = "fast")
            fs::file_delete(data)

            if (Evaluate) {

               if (is.null(Pred_PA)) {
                  PresAbs <- object$Y[, Sp]
               } else {
                  PresAbs <- Pred_PA[, Sp]
               }

               if(length(unique(PresAbs)) == 2) {
                  # Calculate evaluation metrics if there are two both presence
                  # and absence info in (testing) data
                  RMSE <- caret::RMSE(PresAbs, SpDT_Mean)
                  MeanPres <- mean(SpDT_Mean[which(PresAbs == 1)])
                  MeanAbs <- mean(SpDT_Mean[which(PresAbs == 0)])
                  TjurR2 <- MeanPres - MeanAbs
                  AUC <- pROC::auc(
                     response = PresAbs, predictor = SpDT_Mean,
                     levels = c(0, 1), direction = "<") %>%
                     as.numeric()
                  Boyce <- ecospat::ecospat.boyce(
                     fit = SpDT_Mean, obs = SpDT_Mean[PresAbs==1],
                     PEplot = FALSE)$cor
               } else {
                  RMSE <- TjurR2 <- AUC <- Boyce <- NA_real_
               }
            } else {
               RMSE <- TjurR2 <- AUC <- Boyce <- NA_real_
            }

            return(tibble::tibble(
               Sp = Sp, IAS_ID = IAS_ID, Path_pred = PredSummaryFile,
               RMSE = RMSE, AUC = AUC, Boyce = Boyce, TjurR2 = TjurR2))

         },
         future.scheduling = Inf, future.seed = TRUE,
         future.globals = c(
            "Eval_DT", "Evaluate", "object", "PredDir",
            "ModelName", "Pred_PA", "Pred_XY"),
         future.packages = c(
            "dplyr", "Matrix", "purrr", "tibble", "Hmsc", "float", "qs",
            "Rfast", "caret", "pROC", "ecospat", "sf"))

      snow::stopCluster(c1)
      future::plan("future::sequential", gc = TRUE)

      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      Eval_DT <- dplyr::bind_rows(Eval_DT)

      # Save predictions for all species in a single file
      Predictions <- dplyr::select(Eval_DT, Path_pred, Sp, IAS_ID) %>%
         dplyr::mutate(
            Sp_data = purrr::map(
               .x = Path_pred,
               .f = ~ {
                  qs::qread(.x) %>%
                     tidyr::pivot_longer(
                        cols = starts_with("Sp_"),
                        names_to = "Species", values_to = "Prediction")
               })) %>%
         dplyr::pull(Sp_data) %>%
         dplyr::bind_rows() %>%
         tidyr::pivot_wider(
            names_from = "Species", values_from = "Prediction") %>%
         dplyr::relocate(gtools::mixedsort(names(.))) %>%
         dplyr::select(x, y, geometry, tidyselect::everything())

      Pred_File <- file.path(PredDir, paste0("Preds_", ModelName, ".qs"))
      qs::qsave(Predictions, file = Pred_File, preset = "fast")
      fs::file_delete(Eval_DT$Path_pred)

      if (Evaluate) {
         Eval_DT <- dplyr::select(Eval_DT, -Path_pred)
         if (is.null(Evaluate_Name)) {
            save(Eval_DT, file = file.path(EvalDir, "Eval_DT.RData"))
         } else {
            save(
               Eval_DT,
               file = file.path(EvalDir, paste0(Evaluate_Name, ".RData")))
         }
      }
      fs::file_delete(post_file)

      invisible(return(NULL))
   }
}

## ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
## ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

## internal function to get one prediction
##
##  Needs following variables or arguments that must be passed:
##  PiNew X XRRR Yc dfPiNew expected mcmcStep nyNew object pN post
##  predPostEta rL rLPar
#' @export

get1prediction <- function(
      object, X, XRRR, Yc, Loff, rL, rLPar, sam, predPostEta, PiNew, dfPiNew,
      nyNew, expected, mcmcStep, seed = NULL) {

   if (!is.null(seed))
      set.seed(seed)
   if(object$ncRRR>0){
      XB=XRRR%*%t(sam$wRRR)
   }

   switch(
      class(X)[1L],
      matrix = {
         X1=X
         if(object$ncRRR>0){
            X1=cbind(X1,XB)
         }
         LFix = X1 %*% sam$Beta
      },
      list = {
         LFix = matrix(NA, nyNew, object$ns)
         for(j in 1:object$ns){
            X1=X[[j]]
            if(object$ncRRR>0){
               X1=cbind(X1,XB)
            }
            LFix[,j] = X1%*%sam$Beta[,j]
         }
      }
   )

   LRan = vector("list",object$nr)
   Eta = vector("list",object$nr)

   for(r in seq_len(object$nr)){
      Eta[[r]] = predPostEta[[r]]

      if(rL[[r]]$xDim == 0){
         LRan[[r]] = Eta[[r]][as.character(dfPiNew[,r]),] %*% sam$Lambda[[r]]
      } else{
         LRan[[r]] = matrix(0,object$ny,object$ns)
         for(k in 1:rL[[r]]$xDim)
            LRan[[r]] = LRan[[r]] + (Eta[[r]][as.character(dfPiNew[,r]),]*rL[[r]]$x[as.character(dfPiNew[,r]),k]) %*% sam$Lambda[[r]][,,k]
      }
   }

   L = Reduce("+", c(list(LFix), LRan))
   if(!is.null(Loff)) L = L + Loff

   ## predict can be slow with Yc and especially with high mcmcStep
   if(!is.null(Yc) && any(!is.na(Yc))){
      Z = L
      Z = updateZ(Y=Yc, Z=Z, Beta=sam$Beta, iSigma=1/sam$sigma, Eta=Eta,
                  Lambda=sam$Lambda, Loff=Loff, X=X, Pi=PiNew, dfPi=dfPiNew,
                  distr=object$distr, rL=rL)
      ## species CV from computePredictedValues runs this innermost
      ## loop nfolds * nfolds.sp * predN * mcmcStep times
      for(sN in seq_len(mcmcStep)){
         Eta = updateEta(
            Y=Yc, Z=Z, Beta=sam$Beta, iSigma=1/sam$sigma,
            Eta=Eta, Lambda=sam$Lambda, Alpha=sam$Alpha,
            rLPar=rLPar, Loff=Loff, X=X, Pi=PiNew, dfPi=dfPiNew, rL=rL)
         Z = updateZ(Y=Yc, Z=Z, Beta=sam$Beta, iSigma=1/sam$sigma, Eta=Eta,
                     Lambda=sam$Lambda, Loff=Loff, X=X, Pi=PiNew, dfPi=dfPiNew,
                     distr=object$distr, rL=rL)
      }
      for(r in seq_len(object$nr)){
         if(rL[[r]]$xDim == 0){
            LRan[[r]] = Eta[[r]][as.character(dfPiNew[,r]),] %*%
               sam$Lambda[[r]]
         } else{
            LRan[[r]] = matrix(0,object$ny,object$ns)
            for(k in 1:rL[[r]]$xDim)
               LRan[[r]] = LRan[[r]] +
                  (Eta[[r]][as.character(dfPiNew[,r]),] *
                     rL[[r]]$x[as.character(dfPiNew[,r]),k]) %*%
                  sam$Lambda[[r]][,,k]
         }
      }
      L = Reduce("+", c(list(LFix), LRan))
   }

   if(!expected){
      Z = L + matrix(sqrt(sam$sigma),nrow(L),object$ns,byrow=TRUE) * matrix(rnorm(nrow(L)*object$ns),nrow(L),object$ns)
   } else{
      Z = L
   }

   for(j in 1:object$ns){
      if(object$distr[j,"family"] == 2){ # probit
         if(expected){
            Z[,j] = Hmsc::fast_pnorm(Z[,j])
         } else{
            Z[,j] = as.numeric(Z[,j]>0)
         }
      }
      if(object$distr[j,"family"] == 3){ # poisson
         if(expected){
            Z[,j] = exp(Z[,j] + sam$sigma[j]/2)
         } else{
            Z[,j] = rpois(nrow(Z),exp(Z[,j]))
         }
      }
   }
   colnames(Z) = object$spNames

   for(i in 1:object$ns){
      m = object$YScalePar[1,i]
      s = object$YScalePar[2,i]
      if(m!=0 || s!=1){
         Z[,i] = Z[,i]*s + m
      }
   }

   Z
}
