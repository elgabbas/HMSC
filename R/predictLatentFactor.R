#' @title predictLatentFactor
#'
#' @description Draws samples from the conditional predictive
#'     distribution of latent factors
#'
#' @param unitsPred a factor vector with random level units for which
#'     predictions are to be made
#' @param units a factor vector with random level units that are
#'     conditioned on
#' @param postEta a list containing samples of random factors at
#'     conditioned units
#' @param postAlpha a list containing samples of range (lengthscale)
#'     parameters for latent factors
#' @param rL a \code{HmscRandomLevel}-class object that describes the
#'     random level structure
#' @param predictMean a boolean flag indicating whether to return the
#'     mean of the predictive Gaussian process distribution
#' @param predictMeanField a boolean flag indicating whether to return
#'     the samples from the mean-field distribution of the predictive
#'     Gaussian process distribution
#'
#' @return a list of length \code{length(postEta)} containing samples
#'     of random factors at \code{unitsPred} from their predictive
#'     distribution conditional on the values at \code{units}
#'
#' @details Length of \code{units} vector and number of rows in
#'     \code{postEta} matrix shall be equal. The method assumes that
#'     the i-th row of \code{postEta} correspond to i-th element of
#'     \code{units}.
#'
#'   This method uses only the coordinates \code{rL$s} field of the
#'   \code{rL$s} argument. This field shall be a matrix with rownames
#'   covering the union of \code{unitsPred} and \code{units}
#'   factors. Alternatively, it can use distance matrix
#'   \code{rL$distMat} which is a symmetric square matrix with similar
#'   row names as the coordinate data (except for the GPP models that
#'   only can use coordinates).
#'
#'   In case of spatial random level, the computational complexity of
#'   the generic method scales cubically as the number of unobserved
#'   units to be predicted. Both \code{predictMean=TRUE} and
#'   \code{predictMeanField=TRUE} options decrease the asymptotic
#'   complexity to linear. The \code{predictMeanField=TRUE} option
#'   also preserves the uncertainty in marginal distribution of
#'   predicted latent factors, but neglects the inter-dependence
#'   between them.
#'
#' @importFrom stats rnorm dist
#' @importFrom methods is
#' @importFrom FNN knnx.index
#' @importFrom sp spDists
#'
#' @export

predictLatentFactor = function(
      unitsPred, units, postEta, postAlpha, rL, predictMean=FALSE,
      predictMeanField=FALSE, nParallel = 1, TempDir = "TEMP2Pred")    {

   fs::dir_create(TempDir)

   if (inherits(postEta, "character")) {
      postEta <- qs::qread(postEta, nthreads = 5)
   }

   if(predictMean && predictMeanField)
      stop("predictMean and predictMeanField arguments cannot be simultaneously TRUE")

   predN = length(postEta)
   indOld = (unitsPred %in% units)
   indNew = !(indOld)
   n = length(unitsPred)
   np = length(units)
   nn = sum(indNew)

   # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

   # Save postEta as small chunks
   chunk_size <- 25
   ChunkIDs <- ceiling(seq_along(postEta) / chunk_size)
   Chunks <- purrr::map_chr(
      .x = seq_len(max(ChunkIDs)),
      .f = ~ {
         IDs <- which(ChunkIDs == .x)
         Ch <- postEta[IDs]
         ChunkFile <- file.path(TempDir, paste0("postEta_ch", .x, ".qs"))
         qs::qsave(Ch, file = ChunkFile, preset = "fast")
         return(ChunkFile)
      }
   )
   rm(postEta, ChunkIDs)
   invisible(gc())

   # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

   calc_eta_pred <- function(pN, postEta) {

      mm <- getMethod("%*%", "Matrix")

      eta = postEta[[pN]]
      nf = ncol(eta)
      etaPred = matrix(NA, n, nf)
      rownames(etaPred) = unitsPred
      # etaPred[indOld,] = eta[match(unitsPred[indOld],units),]
      etaPred[indOld,] = eta[fastmatch::fmatch(unitsPred[indOld], units), ]

      if(nn > 0){
         if(rL$sDim == 0){
            if(predictMean){
               etaPred[indNew,] = 0
            } else{
               etaPred[indNew,] = matrix(rnorm(sum(indNew)*nf), sum(indNew), nf)
            }
         } else{
            alpha = postAlpha[[pN]]
            alphapw = rL$alphapw

            if(predictMean || predictMeanField){
               if(!is.null(rL$s)){
                  s1 = rL$s[units,, drop=FALSE]
                  s2 = rL$s[unitsPred[indNew],,drop=FALSE]
                  if (is(s1, "Spatial")) {
                     D11 <- spDists(s1)
                     D12 <- spDists(s1, s2)
                  } else {
                     dim = NCOL(s1)
                     # D11 = as.matrix(dist(s1))
                     D11 = Rfast::Dist(s1)
                     invisible(gc())
                     # D12 = sqrt(
                     #    Reduce(
                     #       "+",
                     #       Map(function(i) outer(s1[,i], s2[,i], "-")^2, seq_len(dim))))
                     D12 <- Rfast::dista(s1, s2)
                     invisible(gc())
                  }
               } else {
                  D11 = rL$distMat[units, units, drop=FALSE]
                  D12 = rL$distMat[units, unitsPred[indNew], drop=FALSE]
               }

               for(h in 1:nf){
                  if(alphapw[alpha[h],1] > 0){
                     K11 = exp(-D11/alphapw[alpha[h],1])
                     K12 = exp(-D12/alphapw[alpha[h],1])

                     # m = crossprod(K12, solve(K11, eta[,h]))
                     m = Matrix::crossprod(K12, as.vector(Solve2vect(K11, eta[,h])))

                     if(predictMean)
                        etaPred[indNew,h] = m

                     if(predictMeanField){
                        # LK11 = t(chol(K11))
                        LK11 = Matrix::t(Matrix::chol(K11))

                        # iLK11K12 = solve(LK11, K12)
                        iLK11K12 = Solve2(LK11_2, K12)

                        v = 1 - colSums(iLK11K12^2)
                        etaPred[indNew,h] = m + rnorm(nn, sd=sqrt(v))
                     }
                  } else{
                     if(predictMean)
                        etaPred[indNew,h] = 0
                     if(predictMeanField)
                        etaPred[indNew,h] = rnorm(nn)
                  }
               }

            } else{
               switch(
                  rL$spatialMethod,
                  'Full' = {
                     unitsAll = c(units,unitsPred[indNew])
                     if(!is.null(rL$s)){
                        s = rL$s[unitsAll,,drop=FALSE]
                        if (is(s, "Spatial"))
                           D <- spDists(s)
                        else
                           D = as.matrix(dist(s))
                     } else {
                        D = rL$distMat[unitsAll,unitsAll]
                     }

                     for(h in 1:nf){
                        if(alphapw[alpha[h],1] > 0){
                           K = exp(-D/alphapw[alpha[h],1])
                           K11 = K[1:np,1:np]
                           K12 = K[1:np,np+(1:nn)]
                           K22 = K[np+(1:nn),np+(1:nn)]
                           m = crossprod(K12, solve(K11, eta[,h]))
                           W = K22 - crossprod(K12, solve(K11, K12))
                           L = try(t(chol(W)))
                           if (inherits(L, "try-error")) # assume sd is zero
                              etaPred[indNew,h] <- m
                           else
                              etaPred[indNew,h] = m + L%*%rnorm(nn)
                        } else{
                           etaPred[indNew,h] = rnorm(nn)
                        }
                     }
                  },
                  'NNGP' = {
                     unitsAll = c(units,unitsPred[indNew])
                     indices = list()
                     dist12 = matrix(NA,nrow=rL$nNeighbours,ncol=nn)
                     dist11 = array(NA, c(rL$nNeighbours,rL$nNeighbours,nn))
                     if (is.null(rL$s)) { # distMat instead of coordinates s
                        d <- rL$distMat[unitsAll, unitsAll]
                        indNN <- t(apply(d[np+(1:nn), 1:np], 1, order)[seq_len(rL$nNeighbours),])
                        for(i in 1:nn) {
                           ind <- indNN[i,]
                           indices[[i]] <- rbind(i*rep(1, length(ind)), ind)
                           dist12[,i] <- d[np+i, ind]
                           dist11[,,i] <- d[ind, ind]
                        }
                     } else { # spatial coordinates
                        s = rL$s[unitsAll,,drop=FALSE]
                        sOld = s[1:np,, drop=FALSE]
                        sNew = s[np+(1:nn),, drop=FALSE]
                        ## In Euclidean coordinates we use fast
                        ## FNN::knnx.index, but for Spatial coordinates we
                        ## need to first calculate spatial distances
                        if (is(sOld, "Spatial")) {
                           ## if we use NNGP, full distance matrix can be
                           ## too big, and we loop over sOld rows: this is
                           ## slow but needs less memory
                           nnabo <- rL$nNeighbours
                           indNN <- matrix(0, nn, nnabo)
                           for (i in seq_len(nn)) {
                              indNN[i,] <-
                                 order(spDists(sOld, sNew[i,, drop=FALSE]))[seq_len(nnabo)]
                           }
                        } else {
                           sNew <- as.matrix(sNew)
                           indNN = knnx.index(sOld,sNew,k=rL$nNeighbours)
                        }
                        for(i in 1:nn){
                           ind = indNN[i,]
                           indices[[i]] = rbind(i*rep(1,length(ind)),ind)
                           if (is(sOld, "Spatial")) {
                              dist12[,i] <- spDists(sOld[ind,,drop=FALSE], sNew[i,])
                              dist11[,,i] = spDists(sOld[ind,, drop=FALSE])
                           } else {
                              das <- 0
                              for (dim in seq_len(rL$sDim))
                                 das <- das + (sOld[ind, dim] - sNew[i, dim])^2
                              dist12[,i] <- sqrt(das)
                              dist11[,,i] <- as.matrix(dist(sOld[ind,]))
                           }
                        }
                     }
                     # BgA = list()
                     # FgA = list()

                     for(h in 1:nf){
                        if(alphapw[alpha[h],1] > 0){
                           K12 = exp(-dist12/alphapw[alpha[h],1])
                           ind1 = t(matrix(unlist(indices),nrow=2))[,2]
                           ind2 = t(matrix(unlist(indices),nrow=2))[,1]
                           K11 = exp(-dist11/alphapw[alpha[h],1])
                           K21iK11 = matrix(NA,ncol=nn,nrow=rL$nNeighbours)
                           for(i in 1:nn){
                              iK11 = solve(K11[,,i])
                              K21iK11[,i] = K12[,i]%*%iK11
                           }
                           B = Matrix(0,nrow=nn, ncol=np,sparse=TRUE)
                           B[cbind(ind2,ind1)] = as.vector(K21iK11)
                           Fmat = 1 - colSums(K21iK11*K12)
                           m = B%*%eta[,h]
                           etaPred[indNew,h] = as.numeric(m + sqrt(Fmat)*rnorm(nn))
                        } else{
                           etaPred[indNew,h] = rnorm(nn)
                        }
                     }
                  },
                  "GPP" = {
                     sKnot = rL$sKnot
                     unitsAll = c(units,unitsPred[indNew])
                     s = rL$s[unitsAll,,drop=FALSE]
                     if (is(s, "Spatial")) {
                        das <- spDists(s, sKnot)
                        dss <- spDists(sKnot)
                     } else {
                        dim = NCOL(s)
                        dss = as.matrix(dist(sKnot))
                        das = sqrt(
                           Reduce(
                              "+",
                              Map(function(i) outer(s[,i], sKnot[,i], "-")^2, seq_len(dim))))
                     }
                     dns = das[np+(1:nn),]
                     dnsOld = das[1:np,]

                     for(h in 1:nf){
                        ag = alpha[h]
                        if(alphapw[ag,1]>0){
                           Wns = exp(-dns/alphapw[ag,1])
                           W12 = exp(-dnsOld/alphapw[ag,1])
                           Wss = exp(-dss/alphapw[ag,1])

                           # iWss = solve(Wss)
                           iWss <- Solve1(Wss) #####

                           # WnsiWss = Wns %*% iWss
                           WnsiWss = mm(Wns,  iWss) #####

                           dDn = 1 - rowSums(WnsiWss*Wns)
                           ## dDn can be numerically 0, but negative, say -2.2e-16
                           if (any(dDn < 0))
                              dDn[dDn < 0] <- 0

                           # D = W12 %*% iWss %*% t(W12)
                           D = Matrix::tcrossprod(mm(W12, iWss), W12)

                           dD = 1-diag(D)
                           idD = 1/dD
                           tmp0 = matrix(rep(idD,NROW(sKnot)),ncol=NROW(sKnot))
                           idDW12 = tmp0*W12

                           # FMat = Wss + t(W12) %*% idDW12
                           FMat = Wss + mm(Matrix::t(W12), idDW12)

                           # iF = solve(FMat)
                           iF <- Solve1(FMat)

                           # LiF = chol(iF)
                           LiF = Matrix::chol(iF)

                           # muS1 = iF %*% t(idDW12) %*% eta[, h]
                           muS1 = mm(mm(iF, Matrix::t(idDW12)), eta[, h])

                           # epsS1 = LiF %*% rnorm(nrow(Wss))
                           epsS1 = mm(LiF, rnorm(nrow(Wss)))

                           # m = Wns %*% (muS1 + epsS1)
                           m = mm(Wns, (muS1 + epsS1))

                           etaPred[indNew,h] = as.numeric(m + sqrt(dDn)*rnorm(nn))
                        } else{
                           etaPred[indNew,h] = rnorm(nn)
                        }
                     }
                  })
            }
         }
      }
      return(etaPred)
   }

   # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

   postEtaPred <- if (nParallel == 1) {
      lapply(1:predN, calc_eta_pred)
   } else {
      # withr::local_options(
      # future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)

      c1 <- snow::makeSOCKcluster(nParallel)
      on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
      snow::clusterExport(
         cl = c1,
         list = c(
            "calc_eta_pred", "indOld", "indNew", "n", "np", "nn",
            "unitsPred", "units", "postAlpha", "rL",
            "predictMean", "predictMeanField", "chunk_size", "Chunks"),
         envir = environment())

      invisible(snow::clusterEvalQ(
         cl = c1,
         expr = {
            library(Rcpp); library(RcppArmadillo); library(Matrix)
            library(purrr); library(qs); library(fs)
         }))

      result <- snow::parLapply(
         cl = c1,
         x = seq_len(length(Chunks)),
         fun = function(Chunk){
            ChunkFile <- Chunks[Chunk]
            Out <- purrr::map(
               .x = seq_len(chunk_size), .f = calc_eta_pred,
               postEta = qs::qread(ChunkFile))
            invisible(fs::file_delete(ChunkFile))
            return(Out)
         })
      result <- do.call(c, result)

      snow::stopCluster(c1)
      return(result)
   }

   return(postEtaPred)
}
