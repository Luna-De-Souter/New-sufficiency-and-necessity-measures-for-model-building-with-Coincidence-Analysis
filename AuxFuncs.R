# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Auxiliary functions for 
# Evaluating sufficiency and necessity in                                     # 
# model building with Coincidence Analysis                                    #
# Experiment series on fuzzy-set data                                         #
# [anonymized]                                                                #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

library(cna)
stopifnot(packageVersion("cna") >= "3.5.0")

# Generate environment opsEnv containing custom (fs-)logical operators:
# ->  &:  x&y essentially returns max(x, y)
# ->  |:  x|y essentially returns min(x, y)
# ->  !:  -x return 1-x
# these are used in the function evalConCov defined below
opsEnv <- as.list(cna:::logicalOperators)[c("&", "|", "-")]
names(opsEnv)[[3]] <- "!"
opsEnv <- as.environment(opsEnv)
parent.env(opsEnv) <- baseenv()

# Function evalConCov(): Calculation of con & cov in the customized framework
# 
# tbl is a data frame numeric column X and Y with n*nCond rows
# ccDef is a list of two expressions (definitions of con & cov)
# f are weights/frequencies of length n
# g is an optional grouping factor with nCond distinct values.
# 
# evalConCov() returns a 2-row matrix where rows correspond to con and cov
evalConCov <- function(tbl, ccDef, f, g = NULL){
  if ("x" %in% ccDef$requiredTransf) tbl$x <- evalq(!X, tbl)
  if ("y" %in% ccDef$requiredTransf) tbl$y <- evalq(!Y, tbl)
  if ("XY" %in% ccDef$requiredTransf) tbl$XY <- evalq(X&Y, tbl)
  if ("xy" %in% ccDef$requiredTransf) tbl$xy <- evalq(x&y, tbl)
  if ("Xy" %in% ccDef$requiredTransf) tbl$Xy <- evalq(X&y, tbl)
  if ("xY" %in% ccDef$requiredTransf) tbl$xY <- evalq(x&Y, tbl)
  if ("XYxy" %in% ccDef$requiredTransf) tbl$XYxy <- evalq(X&Y&x&y, tbl)
  if (is.null(g)) g <- rep(1, nrow(tbl))
  nCond <- nrow(tbl) / length(f)
  xySums <- rowsum(tbl*rep.int(f, nCond), g, reorder = FALSE)
  .env <- ccDef$params
  if (is.null(.env)) .env <- baseenv() 
  if (is.list(.env)) .env <- list2env(.env, parent = baseenv())
  rbind(eval(ccDef$con, xySums, .env), 
        eval(ccDef$cov, xySums, .env))
}
environment(evalConCov) <- opsEnv
# Stand-alone appication of evalConCov:
if (F){
  n <- 10
  tbl.in <- data.frame(X = runif(n), Y = runif(10))
  ccDefault <- list(con = quote(XY/X), cov = quote(XY/Y))
  evalConCov(tbl.in, 
             ccDef = ccDefault,
             f = rep(1, 10))
}    

augmentCCdef <- function(ccDef, call_env = parent.frame(2)){
  if (identical(names(ccDef), c("con", "cov", "requiredTransf", "params")))
    return(ccDef)
  if (is.null(ccDef)) 
    ccDef <- list(con = quote(XY/X),  # the default con/cov definitions
                  cov = quote(XY/Y)) 
  `%nin%` <- Negate(`%in%`)

  varsInCCdef <- union(all.vars(ccDef$con), all.vars(ccDef$cov))
  requiredTransf <- intersect(varsInCCdef, 
                              c("X", "Y", "x", "y", "XY", "xy", "Xy", "xY", "XYxy"))
  if ("x" %nin% requiredTransf && any(grepl("x", requiredTransf, fixed = TRUE))){
    requiredTransf <- c("x", requiredTransf)
  }
  if ("y" %nin% requiredTransf && any(grepl("y", requiredTransf, fixed = TRUE))){
    requiredTransf <- c("y", requiredTransf)
  }
  ccDef$requiredTransf <- requiredTransf
  
  parnms <- setdiff(varsInCCdef, requiredTransf)
  ccDef$params <- mget(parnms, call_env)

  ccDef
}

# ------------------------------------------------------------------------------

# clone cna package namespace 
.cna1Nmsp. <- rlang::env_clone(getNamespace("cna"))

# Replace (cna-internal) functions conj_conCov() and C_find_asf() in .cna1Nmsp.
# by 'customizable' versions accepting additional arg ccDef:

# -> conj_conCov:
.cna1Nmsp.$conj_conCov <- function(cols, x, y, f, ccDef = getOption("ccDef")){
  nCond <- nrow(cols)
  n <- nrow(x)
  xscores <- vapply(split.default(cols, row(cols)), function(r) rowMins(x[, r, drop = FALSE]), 
                    FUN.VALUE = numeric(n))
  dxy <- data.frame(Y = rep.int(y, nCond),
                    X = as.vector(xscores))
  conCov <- evalConCov(dxy, ccDef, f = f, 
                       g = rep(seq_len(nCond), each = length(y)))
  rownames(conCov) <- c("con", "cov")
  conCov
}

# -> C_find_asf:
.cna1Nmsp.$C_find_asf <- function (conjlen, x, y, f, con, cov, ..., ccDef = getOption("ccDef")){
  DIM <- function(x) c(NROW(x), NCOL(x))
  ma <- match(conjlen, unique(conjlen))
  tbls <- x[!duplicated(ma)]
  stopifnot(sapply(tbls, ncol)>0)
  n <- length(y)
  multiplicity <- tabulate(ma)
  nlev <- length(tbls)
  a <- mapply(combn, sapply(tbls, ncol), multiplicity, 
              SIMPLIFY = FALSE)
  
  pmlist <- vector("list", nlev)
  for (l in seq_len(nlev)){
    mlist <- apply(a[[l]], 1, function(r) as.vector(tbls[[l]][, r]), simplify = FALSE)
    pmlist[[l]] <- do.call(pmax, mlist)
    dim(pmlist[[l]]) <- c(n, ncol(a[[l]]))
  }

  expand_indices <- do.call(expand.grid, 
                            c(lapply(pmlist, function(x) seq_len(ncol(x))), 
                              list(KEEP.OUT.ATTRS = FALSE)))
  
  chunksize <- 1e5
  nchunks <- (nrow(expand_indices)-1) %/% chunksize + 1
  chunks <- rep(1:nchunks, each = chunksize, length.out = nrow(expand_indices))
  indsList <- vector("list", nchunks)
  
  for (ch in seq_len(nchunks)){
    # process in chunks to reduce memory usage
    .exp_ind <- expand_indices[chunks == ch, , drop = FALSE]
    nCond <- nrow(.exp_ind)
  
    sc <- array(0, dim = c(n, nCond))
    for (l in seq_len(nlev)){
      newsc <- pmlist[[l]][, .exp_ind[[l]]]
      stopifnot(identical(DIM(newsc), DIM(sc)))
      sc[] <- pmax(sc, newsc)
    }
    dxy <- data.frame(Y = rep.int(y, nCond),
                      X = as.vector(sc))
    conCov <- evalConCov(dxy, ccDef, f = f, 
                         g = rep(seq_len(nCond), each = length(y)))
    isSolution <- conCov[1, ] >= con & conCov[2, ] >= cov
    isSolution[is.na(isSolution)] <- FALSE   # case con is NA
    
    indsList[[ch]] <- .exp_ind[isSolution, , drop = FALSE]
  }
  inds <- do.call(rbind, indsList)
  
  out <- mapply(function(tbl, i) tbl[, i, drop = FALSE], a, inds, 
                SIMPLIFY = FALSE)
  out <- t(do.call(rbind, out)) - 1
  mode(out) <- "integer"
  dimnames(out) <- NULL
  out
}

# make.asf
.cna1Nmsp.$make.asf <- function (cti, zname, .sol, inus.only, details, asf.selection, 
                                 verbose = FALSE, ccDef = getOption("ccDef")){
  vnms <- colnames(cti$scores)
  .lhs <- hconcat(.sol, c("+", "*"), f = function(i) vnms[i])
  qcond0 <- qcond_asf(
    paste0(.lhs, "<->", zname), sc = cti$scores) 
  dxy <- data.frame(X = as.vector(qcond0[, 1, ]), 
                    Y = as.vector(qcond0[, 2, ]))
  gr <- rep(seq_len(dim(qcond0)[3]), each = length(cti$freq))
  
  if (asf.selection != "none"){
    varying <- mapply(C_varies, split.default(dxy$X, gr), split.default(dxy$Y, gr), 
                      MoreArgs = list(asfSelection = asf.selection), 
                      USE.NAMES = FALSE)
    if (!all(varying)){
      sel <- gr %in% which(varying)
      dxy <- dxy[sel, , drop = FALSE]
      gr <- gr[sel]
      qcond0 <- structure(
        qcond0[, , varying, drop = FALSE], 
        condition = attr(qcond0, "condition")[varying], 
        response = attr(qcond0, "response")[varying])
    }
  }
  cc <- evalConCov(dxy, ccDef = ccDef, f = cti$freq, g = gr)
  asf <- structure(data.frame(outcome = attr(qcond0, "response"),   # code taken from cna::cond2condTbl
                              condition = structure(attr(qcond0, "condition"), 
                                                    class = c("stdAtomic", "character")), 
                              consistency = cc[1, ], 
                              coverage = cc[2, ], 
                              row.names = NULL, 
                              stringsAsFactors = FALSE), 
                   class = c("condTbl", "data.frame"))

  if (verbose && nrow(asf) < length(.lhs)) {
    cat("    ", length(.lhs) - nrow(asf), " asf (of ", length(.lhs), 
      ") are removed based on criterion asf.selection=\"", 
      asf.selection, "\"\n", sep = "")
  }
  if (nrow(asf) == 0) 
    return(NULL)
  asf$condition[] <- as.character(lhs(asf$condition))
  n_plus_stars <- gregexpr("[\\*\\+]", asf$condition)
  asf$complexity <- lengths(n_plus_stars) + 1L - (vapply(n_plus_stars, 
    "[[", integer(1), 1L) == -1L)
  asf <- asf[order(asf$complexity, -asf$consistency * asf$coverage), 
    , drop = FALSE]
  n0 <- nrow(asf)
  if (length(details)) {
    .cond <- if (n0 == 0) 
      character(0)
    else paste0(asf$condition, "<->", asf$outcome)
    if (useCtiList(cti)) 
      cti <- ctiList(cti, .cond)
    asf <- cbind(asf, .det(cti, .cond, what = details, cycle.type = NULL))
  }
  if (inus.only) {
    asf <- asf[asf$inus, , drop = FALSE]
    if (verbose) 
      cat("    ", n0 - nrow(asf), " non-INUS asf (of ", 
        n0, ") are removed, as inus.only=TRUE\n", sep = "")
  }
  else if (verbose) {
    cat("    Keeping all asf, as inus.only=FALSE\n")
  }
  rownames(asf) <- NULL
  asf
}

# cond2condTbl
.cna1Nmsp.$cond2condTbl <- function(cond, freqs, ccDef = getOption("ccDef")){
  dxy <- data.frame(X = as.vector(cond[, 1, , drop = TRUE]), 
                    Y = as.vector(cond[, 2, , drop = TRUE]))
  gr <- rep(seq_len(dim(cond)[3]), each = dim(cond)[1])
  if (is.null(ccDef)){
    ccDef <- list(con = quote(XY/X), cov = quote(XY/Y)) 
  }
  cc <- evalConCov(dxy, ccDef = ccDef, f = freqs, g = gr)
  structure(
    data.frame(outcome = attr(cond, "response"),
               condition = structure(attr(cond, "condition"),
                                     class = c("stdAtomic", "character")),
               consistency = cc[1, , drop = TRUE],
               coverage = cc[2, , drop = TRUE],
               row.names = NULL,
               stringsAsFactors = FALSE),
    class = c("condTbl", "data.frame"))
}


# change the environment of all functions in .cna1Nmsp. to .cna1Nmsp.
# (in order to make it independent of the package namespace)
cna1_fns <- ls(.cna1Nmsp., all.names = TRUE)
isfn <- sapply(mget(cna1_fns, envir = .cna1Nmsp.), is.function)
cna1_fns <- cna1_fns[isfn]
for (nm in cna1_fns){
  environment(.cna1Nmsp.[[nm]]) <- .cna1Nmsp.
}
rm(cna1_fns, isfn, nm)

# Move evalConCov and augmentCCdef to .cna1Nmsp.
.cna1Nmsp.$evalConCov <- evalConCov
.cna1Nmsp.$augmentCCdef <- augmentCCdef
rm(evalConCov, opsEnv) # augmentCCdef is kept in globalenv


# write a copy of cna() to .cna1Nmsp.
.cna1 <- cna::cna
environment(.cna1) <- .cna1Nmsp.
.cna1Nmsp.$.cna1 <- .cna1
rm(.cna1)

# ------------------------------------------------------------------------------

# condTbl1: cna::condTbl() with additional ccDef argument
condTbl1 <- function(cond, data, ccDef = list(con = quote(XY/X), cov = quote(XY/Y)), 
                     ...){ 
  ccDef <- augmentCCdef(ccDef)
  oldopt <- options(ccDef = ccDef, call_env = parent.frame())
  on.exit(options(oldopt))
  cond0 <- condition(cond, data, ...)
  
  if (!all(attr(cond0, "info")$condTypes %in% c("stdAtomic", "atomic", "stdComplex", "complex"))){
    stop("Inadmissable input 'cond': only asfs and csfs are allowed.")
  }
  
  f <- attributes(cond0[[1]])$n
  out <- as.condTbl(cond0)
  if (nrow(out)) out[, c("consistency", "coverage")] <- NA_real_
  if (any(sel <- attr(cond0, "info")$condTypes %in% c("stdAtomic", "atomic"))){
    # Asfs
    l <- length(cond0)
    stopifnot(lengths(cond0) == 2)
    for (i in which(sel)){
      dxy <- as.data.frame(cond0[[i]])
      attributes(dxy)[c("type", "n", "cases", "info")] <- NULL
      names(dxy) <- c("X", "Y")
      out[i, c("consistency", "coverage")] <- evalConCov(dxy, ccDef, f = f)
    }
  }
  if (any(sel <- attr(cond0, "info")$condTypes %in% c("stdComplex", "complex"))){
    # Csfs
    suppressPackageStartupMessages(requireNamespace("collapse", quietly = TRUE))
    ll <- lengths(cond0[sel], use.names = FALSE)
    n <- nrow(data)
    asfs <- unlist(extract_asf(names(cond0)[sel]))
    condAsf <- array(unlist(cond0[sel], use.names = FALSE), 
                     dim = c(n, 2, sum(ll)))
    attr(condAsf, "condition") <- asfs
    attr(condAsf, "response") <- rhs(asfs)
    ccAsf <- cond2condTbl(condAsf, freqs = f, ccDef = ccDef)
    inds <- rep(seq_along(ll), ll)
    out[sel, "consistency"] <- collapse::fmin.default(ccAsf$consistency, inds, use.g.names = FALSE)
    out[sel, "coverage"] <- collapse::fmin.default(ccAsf$coverage, inds, use.g.names = FALSE)
  }
  out
}
environment(condTbl1) <- .cna1Nmsp.

# Function negateCond():
# Takes a vector of conds (asfs!) and replaces lhs and rhs by their respective negations.
# Uses minimalize(): reshapes lhs only where minimalize() finds a dnf representation
negateCond <- function(cond, try2minimalize = TRUE, ...){
  cond <- noblanks(cond)
  lhs0 <- paste0("!(", lhs(cond), ")")
  rhs0 <- paste0("!(", rhs(cond), ")")
  out <- paste0(lhs0, "<->", rhs0)
  if (try2minimalize){
    lhs1 <- try(minimalize(lhs0, ...), silent = TRUE)
    rhs1 <- try(minimalize(rhs0, ...), silent = TRUE)
    lhs1 <- ifelse(lengths(lhs1)>0, 
                   sapply(lhs1, head, 1), 
                   lhs0)
    rhs1 <- ifelse(lengths(rhs1)>0, 
                   sapply(rhs1, head, 1), 
                   rhs0)
    out <- paste0(lhs1, "<->", rhs1)
  }
  out
}

# ------------------------------------------------------------------------------

# cna1: cna::cna() with additional ccDef argument
cna1 <- function() NULL
formals(cna1) <- c(formals(cna::cna), 
                   list(ccDef = list(con = quote(XY/X), cov = quote(XY/Y))))
body(cna1) <- expression({
  ccDef <- augmentCCdef(ccDef)
  oldopt <- options(ccDef = ccDef, call_env = parent.frame())
  on.exit(options(oldopt))
  cl <- match.call()
  cl$ccDef <- NULL
  cl[[1]] <- as.name(".cna1")
  out <- eval(cl, environment(), parent.frame())
  out$call <- match.call()
  out$ccDef <- ccDef
  out
})
environment(cna1) <- .cna1Nmsp.


# cna2: copy of cna1 with different defaults
# cna2 <- cna1
# formals(cna2)$only.minimal.msc <- FALSE
# formals(cna2)$con.msc <- 0
# formals(cna2)$maxstep <- quote(c(2, 3, 6))
# body(cna2) <- expression({
#   oldopt <- options(ccDef = ccDef, call_env = parent.frame())
#     on.exit(options(oldopt))
#     cl <- match.call()
#   cl$ccDef <- NULL
#   cl[[1]] <- as.name(".cna1")
#   if (missing(only.minimal.msc)) cl$only.minimal.msc <- formals(sys.function())$only.minimal.msc
#   if (missing(con.msc)) cl$con.msc <- formals(sys.function())$con.msc
#   if (missing(maxstep)) cl$maxstep <- formals(sys.function())$maxstep
#   out <- eval(cl, environment(), parent.frame())
#   out$call <- match.call()
#   parnms <- setdiff(union(all.vars(ccDef$con), all.vars(ccDef$cov)), 
#                     c("X", "Y", "x", "y", "XY", "xy", "Xy", "xY", "XYxy"))
#   attr(ccDef, "params") <- mget(parnms, getOption("call_env"))
#   out$ccDef <- ccDef
#   out
# })

# ------------------------------------------------------------------------------

# asf1(): Extension of asf()
# Takes a cna and returns asf(x) with two additional columns:
# - min.con.msc:  minimum of consistency scores of the msc included in asf
# - asf_ok:       logical indicating whether all msc have the required consistency score
asf1 <- function(x, ccDef = x$ccDef, ...){
  ccDef <- augmentCCdef(ccDef)
  stopifnot(inherits(x, "cna"))
  asfs <- asf(x)
  asfCond <- asfs$condition
  if (!identical(ccDef, x$ccDef)){
    asfs[c("consistency", "coverage")] <- 
      condTbl1(asfCond, x$configTable, ccDef = ccDef, ...)[c("consistency", "coverage")]
  }
  extract_msc <- mapply(paste, 
                        strsplit(lhs(asfCond), "+", fixed = TRUE), 
                        rhs(asfCond), 
                        MoreArgs = list(sep = "->"), 
                        SIMPLIFY = FALSE)
  asfs_split <- rapply(extract_msc, 
    function(co) transform(condTbl1(co, x$configTable, ccDef = ccDef, ...), 
                           ok = consistency >= x$con),
    how = "replace")
  asfs$min.con.msc <- sapply(asfs_split, function(x) min(x$consistency))
  asfs$asf_ok <- sapply(asfs_split, function(x) all(x$ok))
  asfs
}

# ------------------------------------------------------------------------------

# csf1(): Extension of csf()
# Takes a cna and returns csf(x) with two additional columns:
csf1 <- function(x, ccDef = x$ccDef, ...){
  stopifnot(inherits(x, "cna"), !is.null(ccDef))
  ccDef <- augmentCCdef(ccDef)
  oldopt <- options(ccDef = ccDef, call_env = parent.frame())
  on.exit(options(oldopt))

  out <- csf(x, ...)
  out
}
environment(csf1) <- .cna1Nmsp.





generate_GTs <- function(num_factors=7, compl=list(sample(1:3,1),sample(3:5,1))){ 
  num_factors <- num_factors #could also become argument
    GT <- randomAsf(num_factors, compl = compl, outcome = "A") # standard list(sample(1:3,1),sample(3:5,1))
   return(GT)
}


generate_GTs_csfs <- function(num_factors=7){ 
  num_factors <- num_factors #could also become argument
  #change in generate_data as well!!!
  # set.seed(seed = seed)
  check <- FALSE
  while(check == FALSE){
    #draw GT
    GT <- randomCsf(num_factors, compl = list(sample(1:3,1),sample(3:5,1)), outcome = c("A", "B"))
    #check if GT is not a single cause structure
    if(getComplexity(GT) != 1){
      check <- TRUE
    }else{ #if GT is a causal structure, change seed to draw new GT
      # seed <- seed + 1
      # set.seed(seed = seed)
      check <- FALSE
    }
  }
  return(GT)
}



makeNoisy <- function(x, GT, ratio, type = c("fs", "cs")) # x: noise-free cs data
{
  x <- ct2df(x)
  rownames(x) <-1:nrow(x) 
  outcome <- cna::rhs(noblanks(GT))
  if(type == "cs"){
    r <- condition(GT, x)
    r <- ct2df(r[[1]])
    #select indices of rows to distort (to make into noisy observations)
    distort <- sample(1:nrow(r), round(nrow(r)*ratio))
    #change outcome values to noisy values
    if(length(distort)==0){noisyDat <- x
    }else{
    for(i in 1:length(distort)){r[distort[i],2] <- 1-r[distort[i],2]}
    #put the rows of r in the same order as x
    r <- r[order(as.numeric(rownames(r)), decreasing = F),]
    x[,which(names(x)==outcome)] <- r[,2]
    noisyDat <- x
    }
  }
  if(type == "fs"){
    exo <-  x[,which(names(x)!=outcome)] 
    exoFS <- ct2df(makeFuzzy(exo, fuzzvalues = seq(0,0.49, 0.02)))
    r <- condition(as.vector(lhs(GT)),exoFS)
    r <- ct2df(r[[1]])
    r <- round(r,2)
    distort <-  round(rnorm(nrow(r), mean = 0, sd = (exp(ratio)+0.2)*ratio ),2)
    outFS <- r+distort
    outFS <- apply(cbind(outFS,rep(1,nrow(outFS))), 1, FUN = min)
    outFS <- apply(cbind(outFS,rep(0,length(outFS))), 1, FUN = max)
    noisyDat <- cbind(exoFS,outFS)
    colnames(noisyDat) <- c(names(exoFS),outcome)  
  }
  noisyDat  
}


generate_data <- function(groundTruth, type= "cs", num_factors=7, sample_size = 200){
  num_factors <- num_factors #could also become argument of function
  sample_size <- sample_size #could also become argument of function
  # set.seed(seed = seed)
  frag_ratio <- sample(seq(0.2,0.5,0.01),1)
  noise_level <- sample(seq(0.01,0.3,0.01),1)
  factors <- rep(2, num_factors)
  fullData <- allCombs(factors) - 1
  idealData <- ct2df(selectCases(groundTruth, fullData)) #DB: groundTruth <- randomAsf(num_factors, compl = list(sample(1:3,1),sample(3:5,1)), outcome = "A")
  # Introduce frag_ratio fragmentation.
  #maybe fragmented data could have a constant outcome. in such cases, keep generating new fragmented 
  #datasets until we get a dataset in which the outcome does vary
  check1 <- FALSE
  while(check1 == FALSE){
      if (round(nrow(idealData)*frag_ratio) == 0){
        fragData <- idealData
      }else{
        fragData <- idealData[-sample(1:nrow(idealData), round(nrow(idealData)*frag_ratio)), ]
      }
      #check if data does not have a constant outcome (OLD)
      if( (sum(fragData$A) != nrow(fragData) & (sum(fragData$A != 0)) )){ #if and only if outcome is 0 all the time, sum(fragData$A != 0) is 0, such that the expression evaluates to False, which is what we want
        check1 <- TRUE
      }else{ #if outcome is constant, change seed to make noisy again
        # seed <- seed + 1
        # set.seed(seed = seed)
        check1 <- FALSE
      }
  }

  # Sample with replacement from the fragmented data. (Making sure every case in the fragmented data is included
  # at least once.)
   if(sample_size - nrow(fragData)>0){
    sampledData <- rbind(fragData, fragData[sample(1:nrow(fragData), sample_size - nrow(fragData), replace = TRUE),])
   }else{sampledData <- fragData}
    #sometimes, noisy data as made here has a constant outcome, in that case, keep generating until we get a dataset
    #in which the outcome does vary
    check <- FALSE
    while(check == FALSE){
      #draw data
      noisyData <- makeNoisy(x=sampledData,GT=groundTruth,ratio=noise_level,type=type)
      #check if data does not have a constant outcome
      #pre-Apr 5:
      # if( (sum(noisyData$A) != nrow(noisyData) & (sum(noisyData$A != 0)) )){
      #new version Apr 5: (I am requiring values on both sides of the 0.5 anchor)
      if( sum(noisyData$A > 0.5) != 0 & sum(noisyData$A < 0.5) != 0 ){
        check <- TRUE
      }else{ #if outcome is constant, change seed to make noisy again
        # seed <- seed + 1
        # set.seed(seed = seed)
        check <- FALSE
      }
    }

    #noisyData <- makeNoisy(sampledData,groundTruth,ratio=noise_level,type="cs")

                                        
  return(list(data.frame(noisyData)))
}



introInflation.with.constNoise <- function(datNoisy, groundTruth, inflate.ratio, t, constant = c("sample", "fragmentation"))
{
datInflate <- vector("list", length(datNoisy))
if(constant=="sample"){
for(i in 1:length(t)){
  if(t[[i]] == toupper(t[[i]])) {
    noise.natural <- round(checkNoise(GT=groundTruth[i],x=datNoisy[[i]]),2)
    # r <- condition(groundTruth[i], datNoisy[[i]])
    # r <- ct2df(r[[1]])
    # # apply(r, 1, function(y) abs(y[1]-y[2])) %>% table
    # difference <- apply(r, 1, function(y) abs(y[1]-y[2]))
    compatible <- ct2df(selectCases(groundTruth[[i]],datNoisy[[i]]))
    incompatible <- dplyr::setdiff(datNoisy[[i]],compatible)
    select1.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% t[[i]])]==1)
    select1.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% t[[i]])]==1)
    select1 <- rbind(select1.compatible,select1.incompatible)
    select0.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% t[[i]])]==0)
    select0.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% t[[i]])]==0)
    select0 <- rbind(select0.compatible,select0.incompatible)
    datA1.compatible <- some(select1.compatible, round(c*inflate.ratio*1-noise.natural), replace = T)
    datA1.incompatible <- some(select1.incompatible, round(c*inflate.ratio*(noise.natural)), replace = T)
    datA0.compatible <- some(select0.compatible, round(c*(1-inflate.ratio)*1-noise.natural), replace = T)
    datA0.incompatible <- some(select0.incompatible, round(c*(1-inflate.ratio)*(noise.natural)), replace = T)
    datInflate[[i]] <- rbind(datA1.compatible, datA1.incompatible, datA0.compatible, datA0.incompatible) 
    
  } else {
    noise.natural <- round(checkNoise(GT=groundTruth[[i]],x=datNoisy[[i]]),2)
    # r <- condition(groundTruth[i], datNoisy[[i]])
    # r <- ct2df(r[[1]])
    # difference <- apply(r, 1, function(y) abs(y[1]-y[2]))
    compatible <- ct2df(selectCases(groundTruth[[i]],datNoisy[[i]]))
    incompatible <- dplyr::setdiff(datNoisy[[i]],compatible)
    select1.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% toupper(t[[i]]))]==0)
    select1.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% toupper(t[[i]]))]==0)
    select1 <- rbind(select1.compatible,select1.incompatible)
    select0.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% toupper(t[[i]]))]==1)
    select0.incompatible  <-  subset(incompatible, incompatible[,which(names(incompatible)%in% toupper(t[[i]]))]==1)
    select0 <- rbind(select0.compatible,select0.incompatible)
    datA1.compatible <- some(select1.compatible, round(c*inflate.ratio*1-noise.natural), replace = T)
    datA1.incompatible <- some(select1.incompatible, round(c*inflate.ratio*(noise.natural)), replace = T)
    datA0.compatible <- some(select0.compatible, round(c*(1-inflate.ratio)*1-noise.natural), replace = T)
    datA0.incompatible <- some(select0.incompatible, round(c*(1-inflate.ratio)*(noise.natural)), replace = T)
    datInflate[[i]] <- rbind(datA1.compatible, datA1.incompatible, datA0.compatible, datA0.incompatible) 
    
  }
}
}
  if(constant=="fragmentation"){
    for(i in 1:length(t)){
      
      if(t[[i]] == toupper(t[[i]])) {
        noise.natural <- round(checkNoise(GT=groundTruth[i],x=datNoisy[[i]]),2)
        # r <- condition(groundTruth[i], datNoisy[[i]])
        # r <- ct2df(r[[1]])
        # # apply(r, 1, function(y) abs(y[1]-y[2])) %>% table
        # difference <- apply(r, 1, function(y) abs(y[1]-y[2]))
        # compatible <- datNoisy[[i]][which(difference==0),]
        # incompatible <- dplyr::setdiff(datNoisy[[i]],compatible)
        compatible <- ct2df(selectCases(groundTruth[[i]],datNoisy[[i]]))
        incompatible <- dplyr::setdiff(datNoisy[[i]],compatible)
        select1.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% t[[i]])]==1)
        select1.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% t[[i]])]==1)
        select1 <- rbind(select1.compatible,select1.incompatible)
        select0.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% t[[i]])]==0)
        select0.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% t[[i]])]==0)
        select0 <- rbind(select0.compatible,select0.incompatible)
        inflate.natural <- nrow(select1)/nrow(datNoisy[[i]])

        if(inflate.ratio>inflate.natural){
          draw <- ((nrow(datNoisy[[i]])*inflate.ratio) - nrow(select1))/(1-inflate.ratio)
          datA1.compatible <- some(select1.compatible, draw*(1-noise.natural), replace = T)
          datA1.incompatible <- some(select1.incompatible, draw*(noise.natural), replace = T)
          datInflate[[i]]  <- rbind(datA1.compatible, datA1.incompatible, datNoisy[[i]]) 
        }else{
          draw <- ((nrow(datNoisy[[i]])*(1-inflate.ratio)) - nrow(select0))/(1-(1-inflate.ratio))
          datA0.compatible <- some(select0.compatible, draw*(1-noise.natural), replace = T)
          datA0.incompatible <- some(select0.incompatible, draw*(noise.natural), replace = T)
          datInflate[[i]]  <- rbind(datA0.compatible,datA0.incompatible, datNoisy[[i]])        
        }
        
      } else {
        noise.natural <- round(checkNoise(GT=groundTruth[[i]],x=datNoisy[[i]]),2)
        # r <- condition(groundTruth[i], datNoisy[[i]])
        # r <- ct2df(r[[1]])
        # difference <- apply(r, 1, function(y) abs(y[1]-y[2]))
        compatible <- ct2df(selectCases(groundTruth[[i]],datNoisy[[i]]))
        incompatible <- dplyr::setdiff(datNoisy[[i]],compatible)
        select1.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% toupper(t[[i]]))]==0)
        select1.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% toupper(t[[i]]))]==0)
        select1 <- rbind(select1.compatible,select1.incompatible)
        select0.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% toupper(t[[i]]))]==1)
        select0.incompatible  <-  subset(incompatible, incompatible[,which(names(incompatible)%in% toupper(t[[i]]))]==1)
        select0 <- rbind(select0.compatible,select0.incompatible)

        inflate.natural <- nrow(select1)/nrow(datNoisy[[i]])

        if(inflate.ratio>inflate.natural){
          draw <- ((nrow(datNoisy[[i]])*inflate.ratio) - nrow(select1))/(1-inflate.ratio)
          datA1.compatible <- some(select1.compatible, draw*(1-noise.natural), replace = T)
          datA1.incompatible <- some(select1.incompatible, draw*(noise.natural), replace = T)
          datInflate[[i]]  <- rbind(datA1.compatible, datA1.incompatible, datNoisy[[i]]) 
        }else{
          draw <- ((nrow(datNoisy[[i]])*(1-inflate.ratio)) - nrow(select0))/(1-(1-inflate.ratio))
          datA0.compatible <- some(select0.compatible, draw*(1-noise.natural), replace = T)
          datA0.incompatible <- some(select0.incompatible, draw*(noise.natural), replace = T)
          datInflate[[i]]  <- rbind(datA0.compatible,datA0.incompatible, datNoisy[[i]]) 
          
        }
        
      }
    }
    
  }
datInflate
}

introInflation.with.constNoise.fs <- function(datNoisy, groundTruth, inflate.ratio, t, constant = c("sample", "fragmentation"))
{
datInflate <- vector("list", length(datNoisy))
if(constant=="sample"){
for(i in 1:length(t)){
  if(t[[i]] == toupper(t[[i]])) {
    noise.natural <- round(checkNoise(GT=groundTruth[i],x=datNoisy[[i]]),2)
    r <- condition(groundTruth[i], datNoisy[[i]])
    r <- ct2df(r[[1]])
    # apply(r, 1, function(y) abs(y[1]-y[2])) %>% table
    difference <- apply(r, 1, function(y) abs(y[1]-y[2]))
    compatible <- datNoisy[[i]][which(difference<=noise.natural+0.15),]
    incompatible <- dplyr::setdiff(datNoisy[[i]],compatible)
    select1.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% t[[i]])]>=0.5)
    select1.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% t[[i]])]>=0.5)
    select1 <- rbind(select1.compatible,select1.incompatible)
    select0.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% t[[i]])]<0.5)
    select0.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% t[[i]])]<0.5)
    select0 <- rbind(select0.compatible,select0.incompatible)
    datA1.compatible <- some(select1.compatible, round(c*inflate.ratio*1-noise.natural), replace = T)
    datA1.incompatible <- some(select1.incompatible, round(c*inflate.ratio*(noise.natural)), replace = T)
    datA0.compatible <- some(select0.compatible, round(c*(1-inflate.ratio)*1-noise.natural), replace = T)
    datA0.incompatible <- some(select0.incompatible, round(c*(1-inflate.ratio)*(noise.natural)), replace = T)
    datInflate[[i]] <- rbind(datA1.compatible, datA1.incompatible, datA0.compatible, datA0.incompatible) 
    
  } else {
    noise.natural <- round(checkNoise(GT=groundTruth[[i]],x=datNoisy[[i]]),2)
    r <- condition(groundTruth[i], datNoisy[[i]])
    r <- ct2df(r[[1]])
    difference <- apply(r, 1, function(y) abs(y[1]-y[2]))
    compatible <- datNoisy[[i]][which(difference<=noise.natural+0.15),]
    incompatible <- dplyr::setdiff(datNoisy[[i]],compatible)
    select1.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% toupper(t[[i]]))]<0.5)
    select1.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% toupper(t[[i]]))]<0.5)
    select1 <- rbind(select1.compatible,select1.incompatible)
    select0.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% toupper(t[[i]]))]>=0.5)
    select0.incompatible  <-  subset(incompatible, incompatible[,which(names(incompatible)%in% toupper(t[[i]]))]>=0.5)
    select0 <- rbind(select0.compatible,select0.incompatible)
    datA1.compatible <- some(select1.compatible, round(c*inflate.ratio*1-noise.natural), replace = T)
    datA1.incompatible <- some(select1.incompatible, round(c*inflate.ratio*(noise.natural)), replace = T)
    datA0.compatible <- some(select0.compatible, round(c*(1-inflate.ratio)*1-noise.natural), replace = T)
    datA0.incompatible <- some(select0.incompatible, round(c*(1-inflate.ratio)*(noise.natural)), replace = T)
    datInflate[[i]] <- rbind(datA1.compatible, datA1.incompatible, datA0.compatible, datA0.incompatible) 
    
  }
}
}
  if(constant=="fragmentation"){
    for(i in 1:length(t)){
      
      if(t[[i]] == toupper(t[[i]])) {
        noise.natural <- round(checkNoise(GT=groundTruth[i],x=datNoisy[[i]]),2)
        r <- condition(groundTruth[i], datNoisy[[i]])
        r <- ct2df(r[[1]])
        # apply(r, 1, function(y) abs(y[1]-y[2])) %>% table
        difference <- apply(r, 1, function(y) abs(y[1]-y[2]))
        compatible <- datNoisy[[i]][which(difference<=noise.natural+0.15),]
        incompatible <- dplyr::setdiff(datNoisy[[i]],compatible)
        select1.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% t[[i]])]>=0.5)
        select1.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% t[[i]])]>=0.5)
        select1 <- rbind(select1.compatible,select1.incompatible)
        select0.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% t[[i]])]<0.5)
        select0.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% t[[i]])]<0.5)
        select0 <- rbind(select0.compatible,select0.incompatible)
        inflate.natural <- nrow(select1)/nrow(datNoisy[[i]])

        if(inflate.ratio>inflate.natural){
          draw <- ((nrow(datNoisy[[i]])*inflate.ratio) - nrow(select1))/(1-inflate.ratio)
          datA1.compatible <- some(select1.compatible, draw*(1-noise.natural), replace = T)
          datA1.incompatible <- some(select1.incompatible, draw*(noise.natural), replace = T)
          datInflate[[i]]  <- rbind(datA1.compatible, datA1.incompatible, datNoisy[[i]]) 
        }else{
          draw <- ((nrow(datNoisy[[i]])*(1-inflate.ratio)) - nrow(select0))/(1-(1-inflate.ratio))
          datA0.compatible <- some(select0.compatible, draw*(1-noise.natural), replace = T)
          datA0.incompatible <- some(select0.incompatible, draw*(noise.natural), replace = T)
          datInflate[[i]]  <- rbind(datA0.compatible,datA0.incompatible, datNoisy[[i]])        
        }
        
      } else {
        noise.natural <- round(checkNoise(GT=groundTruth[[i]],x=datNoisy[[i]]),2)
        r <- condition(groundTruth[i], datNoisy[[i]])
        r <- ct2df(r[[1]])
        difference <- apply(r, 1, function(y) abs(y[1]-y[2]))
        compatible <- datNoisy[[i]][which(difference<=noise.natural+0.15),]
        incompatible <- dplyr::setdiff(datNoisy[[i]],compatible)
        select1.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% toupper(t[[i]]))]<0.5)
        select1.incompatible <-  subset(incompatible, incompatible[,which(names(incompatible)%in% toupper(t[[i]]))]<0.5)
        select1 <- rbind(select1.compatible,select1.incompatible)
        select0.compatible <-  subset(compatible, compatible[,which(names(compatible)%in% toupper(t[[i]]))]>=0.5)
        select0.incompatible  <-  subset(incompatible, incompatible[,which(names(incompatible)%in% toupper(t[[i]]))]>=0.5)
        select0 <- rbind(select0.compatible,select0.incompatible)

        inflate.natural <- nrow(select1)/nrow(datNoisy[[i]])

        if(inflate.ratio>inflate.natural){
          draw <- ((nrow(datNoisy[[i]])*inflate.ratio) - nrow(select1))/(1-inflate.ratio)
          datA1.compatible <- some(select1.compatible, draw*(1-noise.natural), replace = T)
          datA1.incompatible <- some(select1.incompatible, draw*(noise.natural), replace = T)
          datInflate[[i]]  <- rbind(datA1.compatible, datA1.incompatible, datNoisy[[i]]) 
        }else{
          draw <- ((nrow(datNoisy[[i]])*(1-inflate.ratio)) - nrow(select0))/(1-(1-inflate.ratio))
          datA0.compatible <- some(select0.compatible, draw*(1-noise.natural), replace = T)
          datA0.incompatible <- some(select0.incompatible, draw*(noise.natural), replace = T)
          datInflate[[i]]  <- rbind(datA0.compatible,datA0.incompatible, datNoisy[[i]]) 
          
        }
        
      }
    }
    
  }
datInflate
}

# Function to check noise (auxiliary for introInflation.with.constNoise).
checkNoise <- function(x, GT){
  r <- condition(GT, x)
  r <- ct2df(r[[1]])
  r <- round(r,2)
  out <- mean(apply(r, 1, function(y) abs(y[1]-y[2])))
  out
}


quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# Function to remove brackets.
remove_par <- function(string){
  gsub("^\\(|\\)$", "", string)
}

#function to run cna. Currently, a separate function like this will be
#needed for each pair of con cov scores. I would like to pass on the
#ccDef (here called def) as an argument to run_harm, but do not know
#how to make that work.  (see end of script for appendix with demonstration)
#(weirdly, this does only work when x is used to indicate the dataset)
#this function has harmonic means as con cov scores
run_harm <- function(x, con, cov, maxstep){
  def <- list(con = quote( ( (1 + 1^2) * XY/X*xy/y ) / ( 1^2 * XY/X + xy/y) ),
              cov = quote( ( (1 + 1^2) * XY/Y*xy/x ) / ( 1^2 * XY/Y + xy/x) ))
  
  out <- csf(cna1(x, outcome = "A", con = con, cov = cov, maxstep = maxstep,
                  rm.dup.factors = F, rm.const.factors = F,
                  ccDef = def))
  return(out)
}

# #careful! save under different name
# run_harm <- function(x, con, cov, maxstep){
#   #def <- list(con = quote( ( (1 + 1^2) * XY/X*xy/y ) / ( 1^2 * XY/X + xy/y) ),
#   #            cov = quote( ( (1 + 1^2) * XY/Y*xy/x ) / ( 1^2 * XY/Y + xy/x) ))
#   
#   out <- csf(cna(x, outcome = "A", con = con, cov = cov, maxstep = maxstep,
#                   rm.dup.factors = F, rm.const.factors = F))
#   return(out)
# }

#running frscore
run_harm_fr <- function(x, tquant, maxstep){
  def <- list(con = quote( ( (1 + 1^2) * XY/X*xy/y ) / ( 1^2 * XY/X + xy/y) ),
              cov = quote( ( (1 + 1^2) * XY/Y*xy/x ) / ( 1^2 * XY/Y + xy/x) ))
  
  asfs <- frscored_cna1(x, normalize = "truemax", outcome = "A",
                        rm.dup.factors = F, rm.const.factors = F,
                        maxstep = maxstep, ccDef = def)
  #copy-pasted from robustness paper
  out <- lapply(asfs[1], function(y) if (!is.null(y)){
    y[y$score >= quantile(y$score, tquant, na.rm = T),]})
  return(out$rean_models)
}

#### benchmark calculation functions (make sure to apply extra check for use on frscore outputs)

#csfs <- out_fr[[2]][[2]]
#GT <- GTs[[2]]



#function to calculate correctness of an output
corr_calc <- function(csfs, GT){
  if(length(csfs) == 0)
    return(0)
  if(nrow(csfs) == 0)
    return(NA)
  csfs_cr <- as.numeric(is.submodel(csfs$condition, GT))
  cr <- max(csfs_cr)
  return(cr)
}

corr_calc_csf <- function(csfs, GT){
  if(length(csfs) == 0)
    return(0)
  if(nrow(csfs) == 0)
    return(NA)
  csfs_cr <- as.numeric(lapply(csfs$condition, function(x) causal_submodel(x, GT)))
  cr <- max(csfs_cr)
  return(cr)
}


corr_calc_frscore <- function(csfs, GT){
  if(length(csfs) == 0)
    return(0)
  if(nrow(csfs) == 0)
    return(NA)
  csfs_cr <- as.numeric(is.submodel(csfs$model, GT))
  cr <- max(csfs_cr)
  return(cr)
}

corr_calc_frscore_csf <- function(csfs, GT){
  if(length(csfs) == 0)
    return(0)
  if(nrow(csfs) == 0)
    return(NA)
  # csfs_cr <- as.numeric(is.submodel(csfs$model, GT))
   csfs_cr <- as.numeric(lapply(csfs$model, function(x) causal_submodel(x, GT)))
  cr <- max(csfs_cr)
  return(cr)
}


# csfs <- out[[i]][[j]]
# GT <- GTs[[j]]
# complete_calc(csfs,GT)
# 
# asf.overlap(c("A+B <-> C", "a*b <-> C"))
#function to calculate completeness of an output
complete_calc <- function(csfs, GT){
  if(length(csfs) == 0) return(0)
  if(nrow(csfs) == 0) return(NA)  #do we want it like this??
  #get intersections with ground truth (I hope the indexing ([1]) is right, haven't encountered
  #anything with multiple maximal intersections)
  intersections <- lapply(csfs$condition, function(x) asf.overlap(c(x, GT))[1])
  #get complexity of ground truth
  GT_cm <- getComplexity(GT)
  #calculate complexity of intersection divided by complexity of ground truth
  completeness <- sapply(intersections, function(x) getComplexity(x)/GT_cm)
  #take maximum of complexity of csfs in output
  #0 is added, such that in case all intersections are empty, complexity is 0
  c_res <- max(c(unlist(completeness), 0))
  return(c_res)
}


complete_calc_frscore <- function(csfs, GT){
  if(length(csfs) == 0) return(0)
  if(nrow(csfs) == 0) return(NA)  #do we want it like this??
  #get intersections with ground truth (I hope the indexing ([1]) is right, haven't encountered
  #anything with multiple maximal intersections)
  intersections <- lapply(csfs$model, function(x) asf.overlap(c(x, GT))[1])
  #get complexity of ground truth
  GT_cm <- getComplexity(GT)
  #calculate complexity of intersection divided by complexity of ground truth
  completeness <- sapply(intersections, function(x) getComplexity(x)/GT_cm)
  #take maximum of complexity of csfs in output
  #0 is added, such that in case all intersections are empty, complexity is 0
  c_res <- max(c(unlist(completeness), 0))
  return(c_res)
}




complete_calc_csf <- function(csfs, GT){
  if(length(csfs) == 0) return(0)
  if(nrow(csfs) == 0) return(NA)  #do we want it like this??
  #get intersections with ground truth (I hope the indexing ([1]) is right, haven't encountered
  #anything with multiple maximal intersections)
  asfs_GT <- extract_asf(GT) |> unlist() |> as.list()
  asfs_csfs <- extract_asf(csfs$condition)
  asfs_csfs <- lapply(asfs_csfs, function(x) split(as.list(x),rhs(x)))
  asfs_GT <- split(asfs_GT,rhs(asfs_GT))
  # by outcome
  inter <- vector("list",length(asfs_csfs))
  for(i in 1:length(asfs_csfs)){
    aux <- vector("list",length(asfs_csfs[[i]]))
    for(j in 1:length(asfs_csfs[[i]])){
      if(any(names(asfs_GT)%in%names(asfs_csfs[[i]][j]))){
        asfs_GT_select <- asfs_GT[which(names(asfs_GT)%in%names(asfs_csfs[[i]][j]))]
        test <- expand.grid(unlist(asfs_csfs[[i]][[j]]),unlist(asfs_GT_select),stringsAsFactors = F)
        overlap <- apply(test, 1, function(x)asf.overlap(c(x[1], x[2]))[1])
        aux[[j]] <- if(!is.null(overlap)){max(getComplexity(overlap))}else{0}
      }else{aux[[j]] <- 0}
    }
    # if(any(names(asfs_GT)%in%names(asfs_csfs[[i]]))){
    # asfs_GT_select <- asfs_GT[which(names(asfs_GT)%in%names(asfs_csfs[[i]]))]
    # test <- expand.grid(unlist(asfs_csfs[i]),unlist(asfs_GT_select),stringsAsFactors = F)
    # overlap <- apply(test, 1, function(x)asf.overlap(c(x[1], x[2]))[1])
    # inter[[i]] <- if(!is.null(overlap)){max(getComplexity(overlap))}else{0}
    # }else{inter[[i]] <- 0}
    inter[[i]] <- unlist(aux) |> mean() 
    }
  #get complexity of ground truth
  # GT_cm <- getComplexity(GT)
  # #calculate complexity of intersection divided by complexity of ground truth
  # completeness <- sapply(intersections, function(x) getComplexity(x)/GT_cm)
  #take maximum of complexity of csfs in output
  #0 is added, such that in case all intersections are empty, complexity is 0
  c_res <- max(max(unlist(inter)), 0)/getComplexity(GT)
  return(c_res)
}

complete_calc_frscore_csf <- function(csfs, GT){
  if(length(csfs) == 0) return(0)
  if(nrow(csfs) == 0) return(NA)  #do we want it like this??
  #get intersections with ground truth (I hope the indexing ([1]) is right, haven't encountered
  #anything with multiple maximal intersections)
  asfs_GT <- extract_asf(GT) |> unlist() |> as.list()
  asfs_csfs <- extract_asf(csfs$model)
  asfs_csfs <- lapply(asfs_csfs, function(x) split(as.list(x),rhs(x)))
  asfs_GT <- split(asfs_GT,rhs(asfs_GT))
  # by outcome
  inter <- vector("list",length(asfs_csfs))
  for(i in 1:length(asfs_csfs)){
    aux <- vector("list",length(asfs_csfs[[i]]))
    for(j in 1:length(asfs_csfs[[i]])){
      if(any(names(asfs_GT)%in%names(asfs_csfs[[i]][j]))){
        asfs_GT_select <- asfs_GT[which(names(asfs_GT)%in%names(asfs_csfs[[i]][j]))]
        test <- expand.grid(unlist(asfs_csfs[[i]][[j]]),unlist(asfs_GT_select),stringsAsFactors = F)
        overlap <- apply(test, 1, function(x)asf.overlap(c(x[1], x[2]))[1])
        aux[[j]] <- if(!is.null(overlap)){max(getComplexity(overlap))}else{0}
      }else{aux[[j]] <- 0}
    }
    # if(any(names(asfs_GT)%in%names(asfs_csfs[[i]]))){
    # asfs_GT_select <- asfs_GT[which(names(asfs_GT)%in%names(asfs_csfs[[i]]))]
    # test <- expand.grid(unlist(asfs_csfs[i]),unlist(asfs_GT_select),stringsAsFactors = F)
    # overlap <- apply(test, 1, function(x)asf.overlap(c(x[1], x[2]))[1])
    # inter[[i]] <- if(!is.null(overlap)){max(getComplexity(overlap))}else{0}
    # }else{inter[[i]] <- 0}
    inter[[i]] <- unlist(aux) |> mean() 
    }
  #get complexity of ground truth
  # GT_cm <- getComplexity(GT)
  # #calculate complexity of intersection divided by complexity of ground truth
  # completeness <- sapply(intersections, function(x) getComplexity(x)/GT_cm)
  #take maximum of complexity of csfs in output
  #0 is added, such that in case all intersections are empty, complexity is 0
  c_res <- max(max(unlist(inter)), 0)/getComplexity(GT)
  return(c_res)
}


#function to calculate proportion of correct models in an output
corr_prop_calc <- function(csfs, GT){
  if(length(csfs) == 0)
    return(NA)
  if(nrow(csfs) == 0)
    return(NA)    #do we want it like this? choice between NA and 0
  #probably should go for NA
  csfs_cr <- as.numeric(is.submodel(csfs$condition, GT))
  cr <- mean(csfs_cr)
  return(cr)
}


corr_prop_calc_csf <- function(csfs, GT){
  if(length(csfs) == 0)
    return(NA)
  if(nrow(csfs) == 0)
    return(NA)    #do we want it like this? choice between NA and 0
  #probably should go for NA
  csfs_cr <- as.numeric(lapply(csfs$condition, function(x) causal_submodel(x, GT)))
  cr <- max(csfs_cr)
  return(cr)
}

corr_prop_calc_csf_total <- function(csfs, GT){
  if(length(csfs) == 0)
    return(NA)
  # if(nrow(csfs) == 0)
  #   return(NA)    #do we want it like this? choice between NA and 0
  #probably should go for NA
  yy <- lapply(csfs,function(x)causal_submodel(x,GT)) |> unlist()
  csfs_cr <- as.numeric(yy)
  cr <- mean(csfs_cr)
  return(cr)
}



corr_prop_calc_frscore <- function(csfs, GT){
  if(length(csfs) == 0)
    return(NA)
  if(nrow(csfs) == 0)
    return(NA)    #do we want it like this? choice between NA and 0
  #probably should go for NA
  csfs_cr <- as.numeric(is.submodel(csfs$model, GT))
  cr <- mean(csfs_cr)
  return(cr)
}


corr_prop_calc_frscore_csf <- function(csfs, GT){
  if(length(csfs) == 0)
    return(NA)
  if(nrow(csfs) == 0)
    return(NA)    #do we want it like this? choice between NA and 0
  #probably should go for NA
  csfs_cr <- as.numeric(lapply(csfs$model, function(x) causal_submodel(x, GT)))
  cr <- mean(csfs_cr)
  return(cr)
}


# 
# corr_calc_frscore_csf <- function(csfs, GT){
#   if(length(csfs) == 0)
#     return(0)
#   if(nrow(csfs) == 0)
#     return(NA)
#   # csfs_cr <- as.numeric(is.submodel(csfs$model, GT))
#    csfs_cr <- as.numeric(lapply(csfs$model, function(x) causal_submodel(x, GT)))
#   cr <- max(csfs_cr)
#   return(cr)
# }



#function to calculate ambiguity of an output
amb_calc <- function(csfs){
  if (length(csfs) == 0)
    return(NA)
  if (nrow(csfs) == 0)
    return(NA)
  return(length(csfs$condition))
}

amb_calc_frscore <- function(csfs){
  if (length(csfs) == 0)
    return(NA)
  if (nrow(csfs) == 0)
    return(NA)
  return(length(csfs$model))
}







# ------------------------------------------------------------------------------

library(stringr)

# asf1(): Extension of asf()
# Takes a cna and returns asf(x) with two additional columns:
# - min.con.msc:  minimum of consistency scores of the msc included in asf
# - asf_ok:       logical indicating whether all msc have the required consistency score
asf1 <- function(x, ccDef = x$ccDef, ...){
  stopifnot(inherits(x, "cna"))
  asfs <- asf(x)
  if (!is.null(params <- attr(ccDef, "params"))){
    .env <- params
  } else {
    .env <- parent.frame()
  }
  asfCond <- asfs$condition
  if (!identical(ccDef, x$ccDef)){
    asfs[c("consistency", "coverage")] <- 
      condTbl1(asfCond, x$configTable, ccDef = ccDef, ..., 
               .env = .env)[c("consistency", "coverage")]
  }
  extract_msc <- mapply(paste, 
                        strsplit(lhs(asfCond), "+", fixed = TRUE), 
                        rhs(asfCond), 
                        MoreArgs = list(sep = "->"))
  asfs_split <- rapply(extract_msc, 
    function(co) transform(condTbl1(co, x$configTable, ccDef = ccDef, ..., 
                                    .env = .env), 
                           ok = consistency>=x$con),
    how = "replace")
  asfs$min.con.msc <- sapply(asfs_split, function(x) min(x$consistency))
  asfs$asf_ok <- sapply(asfs_split, function(x) all(x$ok))
  asfs
}

rean_cna1 <- function (x, attempt = seq(1, 0.7, -0.1), ncsf = 20, output = c("csf", 
    "asf"), ...) 
{
    if (!inherits(x, c("configTable", "data.frame", "truthTab"))) {
    }
    cl <- match.call()
    dots <- list(...)
    if (any(c("cov", "con", "con.msc") %in% names(dots))) {
        abort("cna arguments 'con', 'cov', 'con.msc' not meaningful")
    }
    output <- match.arg(output)
    cl$attempt <- cl$asf <- cl$ncsf <- cl$csf <- cl$output <- NULL
    cl[[1]] <- as.name("cna1")
    cl$what <- if (output == "asf") {
        "a"
    }
    else {
        "c"
    }
    ccargs <- as.data.frame(expand.grid(attempt, attempt))
    colnames(ccargs) <- c("lowfirst", "lowsec")
    sols <- vector("list", length = nrow(ccargs))
    for (i in 1:length(sols)) {
        cl$con <- ccargs[i, "lowfirst"]
        cl$cov <- ccargs[i, "lowsec"]
        if (output == "csf") {
            sols[[i]] <-  asf(eval.parent(cl))
        }
        if (output == "asf") {
            sols[[i]] <- asf(eval.parent(cl))
        }
        dt <- data.frame(cnacon = rep(cl$con, nrow(sols[[i]])), 
            cnacov = rep(cl$cov, nrow(sols[[i]])))
        sols[[i]] <- cbind(sols[[i]], dt)
    }
    return(structure(sols, class = c("rean_cna", "list")))
}


frscored_cna1 <- function (x, fit.range = c(1, 0.7), granularity = 0.1, output = c("csf", 
    "asf"), scoretype = c("full", "supermodel", "submodel"), 
    normalize = c("truemax", "idealmax", "none"), verbose = FALSE, 
    maxsols = 50, test.model = NULL, print.all = FALSE, ...) 
{
    if (!inherits(x, c("configTable", "data.frame", "truthTab"))) {
        stop("invalid argument x")
    }
    cl <- match.call()
    dots <- list(...)
    if (any(c("cov", "con", "con.msc") %in% names(dots))) {
        abort("cna arguments 'con', 'cov', 'con.msc' not meaningful")
    }
    cl$fit.range <- cl$granularity <- cl$normalize <- cl$verbose <- cl$scoretype <- cl$test.model <- cl$print.all <- cl$scoretype <- cl$maxsols <- NULL
    cl[[1]] <- as.name("rean_cna1")
    if ("ncsf" %in% names(dots)) {
        cl$ncsf <- dots$ncsf
    }
    attempt <- seq(max(fit.range), min(fit.range), -granularity)
    cl$attempt <- attempt
    cl$output <- match.arg(output)
    clres <- eval.parent(cl)
    rescomb <- do.call(rbind, clres)
    rescomb <- rescomb[!is.na(rescomb[, 1]), ]
    rescombtemp <- rescomb
    rescomb <- rescomb[, -c(which(names(rescomb) %in% c("cnacon", 
        "cnacov")))]
    rescomb$condition <- gsub("\\),\\(", "\\)*\\(", rescomb$condition)
    scoretype <- match.arg(scoretype)
    normalize <- match.arg(normalize)
    if (is.null(test.model)) {
        scored <- frscore(rescomb$condition, normalize = normalize, 
            verbose = verbose, scoretype = scoretype, maxsols = maxsols)
        if (is.null(scored)) {
            warning("no solutions found in reanalysis series, perhaps consider wider fit range \n \n")
            return(NULL)
        }
    }
    else {
        if (any(sapply(rescomb$condition, function(x) cna::identical.model(x, 
            test.model)))) {
            scored <- frscore(rescomb$condition, normalize = normalize, 
                verbose = verbose, scoretype = scoretype, maxsols = maxsols)
            if (is.null(scored)) {
                warning("no solutions found in reanalysis series, perhaps consider wider fit range \n \n")
                return(NULL)
            }
        }
        else {
            abort("`test.model` not found in reanalysis series")
        }
    }
    sc <- scored[[1]]
    names(sc)[names(sc) == "model"] <- "condition"
    rescomb$condition <- as.character(rescomb$condition)
    rescombXscored <- dplyr::right_join(rescomb, sc, by = "condition") %>% 
        dplyr::filter(!is.na(.data$score))
    rescombXscored <- unique(rescombXscored)
    rescombXscored <- rescombXscored[order(rescombXscored$condition), 
        ]
    rescombXscored <- rescombXscored[order(rescombXscored$score, 
        decreasing = T), ]
    rownames(rescombXscored) <- 1:nrow(rescombXscored)
    if (!is.null(test.model)) {
        tested <- rescombXscored[sapply(rescombXscored$condition, 
            function(x) cna::identical.model(x, test.model)), 
            ]
    }
    else {
        tested <- test.model
    }
    out <- structure(list(rean_models = rescombXscored, tested = tested, 
        verbose = scored$verbose, verbout = scored$verbout, print.all = print.all, 
        fit.range = fit.range, granularity = granularity, scoretype = scoretype, 
        normal = normalize, rean.results = rescombtemp, maxsols = scored$maxsols), 
        class = c("frscored_cna", "list"))
    return(out)
}


# switch.case <- function(x){
# ifelse(useful::find.case(x), tolower(x),toupper(x))
# }


standardize <- function(GT){
  xx <- unlist(extract_asf(GT)) |> lhs() 
  yy <- unlist(extract_asf(GT)) |> rhs() 
  xx <- stdCond(xx)
  out <- mapply(function(x,y) paste0("(", x,"<->", y,")"), x=xx,y=yy) |> as.vector()
  out <- paste(out, collapse="*")
  out
}

# x <- randomAsf(6)

negate <- function(x, seed){
  left <- lhs(x)
  right <- rhs(x)
  set.seed(seed = seed)
  neg.left <- rreduce(getCond(selectCases((paste0("!","(",left,")")))), maxiter = 5000)
  out <- paste0(neg.left, "<->", right)
  out
}

# x <- randomAsf(6, compl = 2:5)
# dat1 <- selectCases(x) |> ct2df()
# dat2 <- selectCases(negate(x)) |> ct2df()
# x
# negate(x)
# out <- rhs(x)
# sum(dat1[names(dat1)%in% out])/nrow(dat1)
# sum(dat2[names(dat2)%in% out])/nrow(dat2)

condType <- function(cond, x = full.ct(cond)){
  out <- cna:::getCondType(cond, cna:::ctInfo(configTable(x)))
  as.vector(out)
}


stdCond2 <- function (x) 
{
    l <- strsplit(noblanks(x), "+", fixed = TRUE)
    if (length(l) == 0L) 
        return(character(0))
    # l <- lapply(l, unique.default)
    u <- unlist(l, use.names = FALSE, recursive = FALSE)
    out <- cna:::C_relist_Char(cna:::C_mconcat(strsplit(u, "*", fixed = TRUE), 
        "*", sorted = TRUE), lengths(l))
    cna:::C_mconcat(out, "+", sorted = TRUE)
}





asf.overlap <- function (x) # x is a vector of asfs with identical outcome 
{

     x <- noblanks(x)
     if (length(x) == 1) {
        max.intersect=lhs(x)
    }
    if (!all(condType(x) %in% "stdAtomic")) {
        stop("x must be a character vector of asf only")
    }
   if (length(unique(rhs(x)))>1) {
        stop("All asf in x must have the same outcome")
    }
    # isolate lhs and standardize syntax
    xx <- stdCond(lhs(x))
    dis <- str_split(xx,"\\+")
    max.overlap.dis <- min(unlist(lapply(dis,length)))
    conj <- lapply(dis, function(y) str_split(y, "\\*"))
    xx.grid <- expand.grid(conj,stringsAsFactors = F)
    xx.grid <- as.data.frame(xx.grid)
    attr(xx.grid, "out.attrs") <- NULL
    # str(xx.grid[1,])
    maxlen <- max(sapply(xx.grid,length))
    xx.inter <- lapply(seq(maxlen),function(i) Reduce(intersect,lapply(xx.grid,"[[",i)))
    xx.grid$overlap <- xx.inter
    xx.grid$length <- apply(xx.grid,1, function(y) length(y$overlap))
    xx.select <- xx.grid[which(xx.grid$length>0),]
    xx.select <- xx.select[order(xx.select$length),]

      # build overlaps
    b <- xx.select$overlap 
    b <- lapply(b, function(x) paste(x, collapse="*"))
    b <- unlist(b)
    if(length(b)==0){max.intersect <- NULL
    }else{
    if(length(b)>1){
    tt <- xx.select[,-which(names(xx.select)%in%c("overlap","length"))]
    n <- min(length(b), max.overlap.dis)
    xx.overlap <- vector("list",n)
    for(i in 1:n){
    xx.overlap[[i]] <- combn(1:length(xx.select$overlap), i ,simplify=F)
      }
  

    xx.overlap <- unlist(xx.overlap,recursive = F)
    xx.overlap <- xx.overlap[unlist(lapply(xx.overlap, function(x)!any(apply(tt[x,], 2, function(x)duplicated(x)))))]
    if(length(xx.overlap)>0){
      overlaps <- lapply(xx.overlap,function(x) paste(b[unlist(x)], collapse="+"))
      com <- max(nchar(unlist(overlaps)))
      overlaps <- c(b[nchar(b)==com],overlaps)
      }else{overlaps <- b}
    overlaps <- unlist(overlaps)
    }else{overlaps <- b}
    max.overlaps <- overlaps[which(nchar(overlaps)==max(nchar(overlaps)))] 
    
    max.intersect <- unique(stdCond2(max.overlaps))
    }
    return(max.intersect) 
    
}

empty.output <- structure(list(structure(list(outcome = structure(character(0), class = c("outcomeString", 
"character")), condition = structure(character(0), class = c("condString", 
"character")), consistency = numeric(0), coverage = numeric(0), 
    complexity = numeric(0), inus = logical(0), cnacon = numeric(0), 
    cnacov = numeric(0)), class = "data.frame", row.names = integer(0)), 
    structure(list(outcome = structure(character(0), class = c("outcomeString", 
    "character")), condition = structure(character(0), class = c("condString", 
    "character")), consistency = numeric(0), coverage = numeric(0), 
        complexity = numeric(0), inus = logical(0), cnacon = numeric(0), 
        cnacov = numeric(0)), class = "data.frame", row.names = integer(0)), 
    structure(list(outcome = structure(character(0), class = c("outcomeString", 
    "character")), condition = structure(character(0), class = c("condString", 
    "character")), consistency = numeric(0), coverage = numeric(0), 
        complexity = numeric(0), inus = logical(0), cnacon = numeric(0), 
        cnacov = numeric(0)), class = "data.frame", row.names = integer(0)), 
    structure(list(outcome = structure(character(0), class = c("outcomeString", 
    "character")), condition = structure(character(0), class = c("condString", 
    "character")), consistency = numeric(0), coverage = numeric(0), 
        complexity = numeric(0), inus = logical(0), cnacon = numeric(0), 
        cnacov = numeric(0)), class = "data.frame", row.names = integer(0)), 
    structure(list(outcome = structure(character(0), class = c("outcomeString", 
    "character")), condition = structure(character(0), class = c("condString", 
    "character")), consistency = numeric(0), coverage = numeric(0), 
        complexity = numeric(0), inus = logical(0), cnacon = numeric(0), 
        cnacov = numeric(0)), class = "data.frame", row.names = integer(0)), 
    structure(list(outcome = structure(character(0), class = c("outcomeString", 
    "character")), condition = structure(character(0), class = c("condString", 
    "character")), consistency = numeric(0), coverage = numeric(0), 
        complexity = numeric(0), inus = logical(0), cnacon = numeric(0), 
        cnacov = numeric(0)), class = "data.frame", row.names = integer(0)), 
    structure(list(outcome = structure(character(0), class = c("outcomeString", 
    "character")), condition = structure(character(0), class = c("condString", 
    "character")), consistency = numeric(0), coverage = numeric(0), 
        complexity = numeric(0), inus = logical(0), cnacon = numeric(0), 
        cnacov = numeric(0)), class = "data.frame", row.names = integer(0)), 
    structure(list(outcome = structure(character(0), class = c("outcomeString", 
    "character")), condition = structure(character(0), class = c("condString", 
    "character")), consistency = numeric(0), coverage = numeric(0), 
        complexity = numeric(0), inus = logical(0), cnacon = numeric(0), 
        cnacov = numeric(0)), class = "data.frame", row.names = integer(0)), 
    structure(list(outcome = structure(character(0), class = c("outcomeString", 
    "character")), condition = structure(character(0), class = c("condString", 
    "character")), consistency = numeric(0), coverage = numeric(0), 
        complexity = numeric(0), inus = logical(0), cnacon = numeric(0), 
        cnacov = numeric(0)), class = "data.frame", row.names = integer(0)), 
    structure(list(outcome = structure(character(0), class = c("outcomeString", 
    "character")), condition = structure(character(0), class = c("condString", 
    "character")), consistency = numeric(0), coverage = numeric(0), 
        complexity = numeric(0), inus = logical(0), cnacon = numeric(0), 
        cnacov = numeric(0)), class = "data.frame", row.names = integer(0)), 
    structure(list(outcome = structure(character(0), class = c("outcomeString", 
    "character")), condition = structure(character(0), class = c("condString", 
    "character")), consistency = numeric(0), coverage = numeric(0), 
        complexity = numeric(0), inus = logical(0), cnacon = numeric(0), 
        cnacov = numeric(0)), class = "data.frame", row.names = integer(0)), 
    structure(list(outcome = structure(character(0), class = c("outcomeString", 
    "character")), condition = structure(character(0), class = c("condString", 
    "character")), consistency = numeric(0), coverage = numeric(0), 
        complexity = numeric(0), inus = logical(0), cnacon = numeric(0), 
        cnacov = numeric(0)), class = "data.frame", row.names = integer(0)), 
    structure(list(outcome = structure(character(0), class = c("outcomeString", 
    "character")), condition = structure(character(0), class = c("condString", 
    "character")), consistency = numeric(0), coverage = numeric(0), 
        complexity = numeric(0), inus = logical(0), cnacon = numeric(0), 
        cnacov = numeric(0)), class = "data.frame", row.names = integer(0)), 
    structure(list(outcome = structure(character(0), class = c("outcomeString", 
    "character")), condition = structure(character(0), class = c("condString", 
    "character")), consistency = numeric(0), coverage = numeric(0), 
        complexity = numeric(0), inus = logical(0), cnacon = numeric(0), 
        cnacov = numeric(0)), class = "data.frame", row.names = integer(0)), 
    structure(list(outcome = structure(character(0), class = c("outcomeString", 
    "character")), condition = structure(character(0), class = c("condString", 
    "character")), consistency = numeric(0), coverage = numeric(0), 
        complexity = numeric(0), inus = logical(0), cnacon = numeric(0), 
        cnacov = numeric(0)), class = "data.frame", row.names = integer(0)), 
    structure(list(outcome = structure(character(0), class = c("outcomeString", 
    "character")), condition = structure(character(0), class = c("condString", 
    "character")), consistency = numeric(0), coverage = numeric(0), 
        complexity = numeric(0), inus = logical(0), cnacon = numeric(0), 
        cnacov = numeric(0)), class = "data.frame", row.names = integer(0))), class = c("rean_cna", 
"list"), GT = structure("e*C*F+b*C*G+C*B*g<->A", class = c("stdAtomic", 
"character")))
