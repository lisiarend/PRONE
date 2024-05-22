# This software is released unde The GNU General Public License:
# http://www.gnu.org/licenses/gpl.html

# EigenMS normalization
# Ref: "Normalization of peak intensities in bottom-up MS-based proteomics using
#       singular value decomposition" Karpievitch YV, Taverner T, Adkins JN,
#       Callister SJ, Anderson GA, Smith RD, Dabney AR. Bioinformatics 2009
# and
# Ref: "Metabolomics data normalization with EigenMS"
#       Karpievitch YK, Nikolic SB, Wilson R, Sharman JE, Edwards LM
#       Submitted to PLoS ONE
#
# Here we allow multiple facotrs to be preserved in ANOVA model before identifying bais.
# This requires a bit of thought on the side of the researchers, only few factors should be 'preserved'.
# For example 'treatment group' is important to preserve but 'age' may or may not be important to preserve.
# Here we do not utilize peptides with 1+ grp missing completely, this is a separate problem
# addressed by Wang et al, 2011
#
# Written by Yuliya Karpievitch, Tom Taverner, and Shelley Herbrich
# email: yuliya.k@gmail.com
#
# Version has been split into 2 main functions: eig_norm1 finds significant bias trends,
# then the user can decide if he/she wants to use that number.
# For Matabolomics data we suggest useing 20% of the number of samples examined
# as the number of bias trends to eliminated. By setting
# ints_eig1$h.c = ceil(.2 * num_samples)
# if the number of automatically identified bias trends are close to that number we suggest
# using the estimated number of bias trends.
#
# eig_norm2 function normalizes the data
# Note that rescaling has been abandoned in the latest version due to discussions abot the fact
# that systematic bias is what has been inadvertently added, thus we need to remove the bias
# but nto add any additional variation.
# Metabolomics in particular has a lot of variation that still remains that we concluded
# there is not need for rescaling.
#
# EigenMS estimates and preserves fixed effects.
# Contributions to using Random effects are welcome.

# HOWTO RUN:
# source('~/EigenMS/EigenMS.R')
# dd = # read in the data
# grps = # read in group information file
# logInts = # subset the dd matrix to only hte portion that contains intensities
#           # replace 0's with NAs, do a log-transformation to approximate Normality.
# prot.info = cbind(data.frame(rownames(dd)), data.frame(rownames(dd))) # peptideIDs, no PR IDs here, just duplicate column
#             # in case od protein IDs, those are not importnat for normalization and will be ignored.
# ints_eig1 = eig_norm1(m=logInts, treatment=grps, prot.info=prot.info)
# ints_norm = eig_norm2(rv=ints_eig1)
#
# eig_norm1 = function(m, treatment, prot.info, write_to_file='')
#       m - m x n matrix of log-transformed intensities, num peptides x num samples
#       treatment - either a single factor or a data.frame of factors (actual R factors, else code will fail)
#                   eg:  bl = gl(3,14,168) # 3 factors, repeat every 14 times, total of 168 samples
#                        it = gl(2,7,168)  # 2 factors, repeat every 7 samples, 168 total
#                        Pi = gl(2,42,168) # 2 factors, repeat twice: 42 of PI+, PI-, PI+, PI-
#                       grpFactors = data.frame(bl,it,Pi)# factors we would like to preserve in EigenMS
#       prot.info - 2 column data frame with peptide and protein IDs in that order.
#                   for metabolites both columns should contain metabolite IDs.
#       write_to_file -  if a string is passed in, 'complete' peptides will be written to that file name
#                       Some peptides could be eliminated due to too many missing values (if not able to do ANOVA)

# SUPPLEMENTARY FUNCTION USED BY EIGENMS
# plot top 3 eigentrends with a line at 0.
# in cse of a single treatment group u and v matrices are switched..
# plot.eigentrends = function(svdr, title1) {
#   v = svdr$v
#   d = svdr$d
#   ss = d ^ 2
#   Tk = signif(ss / sum(ss) * 100, 2)
#
#   titles = paste("Trend ", 1:3, " (", Tk[1:3], "%)", sep = "")
#   do.text = function(j)
#     mtext(titles[j],
#           cex = 0.7,
#           padj = -0.7,
#           adj = 1)
#   range.y = range(as.numeric(v[, 1:3]), na.rm = T)
#
#   toplot1_1 = as.numeric(v[, 1])
#   toplot1_2 = as.numeric(v[, 2])
#   toplot1_3 = as.numeric(v[, 3])
#
#   plot(
#     c(1:length(toplot1_1)),
#     toplot1_1,
#     type = 'b',
#     ann = F,
#     ylim = range.y
#   )
#   do.text(1)
#   abline(h = 0, lty = 3)
#   title(
#     title1,
#     cex.main = 1.2,
#     font.main = 1,
#     col.main = "purple",
#     ylab = NULL
#   )
#   plot(
#     c(1:length(toplot1_2)),
#     toplot1_2,
#     type = 'b',
#     ann = F,
#     ylim = range.y
#   )
#   do.text(2)
#   abline(h = 0, lty = 3)
#   plot(
#     c(1:length(toplot1_3)),
#     toplot1_3,
#     type = 'b',
#     ann = F,
#     ylim = range.y
#   )
#   do.text(3)
#   abline(h = 0, lty = 3)
#   return(Tk)
# }
#
# # plot any 3 eigentrends with a line at 0. Strting at the trend number passed in
# # pos1 parameter provide the starting index for teh 3 trends
# # in cse of a single treatment group u and v matrices are switched..
# plot.eigentrends.start = function(svdr, title1, pos1 = 1) {
#   # No check for valid range of pos1 is performed!!!
#   v = svdr$v
#   d = svdr$d
#   ss = d ^ 2
#   Tk = signif(ss / sum(ss) * 100, 2)
#   #  pe = signif(d/sum(d, na.rm=T)*100, 2)
#   titles = paste("Trend ", pos1:(pos1 + 3), " (", Tk[pos1:(pos1 + 3)], "%)", sep = "")
#   do.text = function(j)
#     mtext(titles[j],
#           cex = 0.7,
#           padj = -0.7,
#           adj = 1)
#   range.y = range(as.numeric(v[, pos1:(pos1 + 3)]), na.rm = T)
#
#   toplot1_1 = as.numeric(v[, pos1])
#   toplot1_2 = as.numeric(v[, (pos1 + 1)])
#   toplot1_3 = as.numeric(v[, (pos1 + 2)])
#
#   plot(
#     c(1:length(toplot1_1)),
#     toplot1_1,
#     type = 'b',
#     ann = F,
#     ylim = range.y
#   )
#   do.text(1)
#   abline(h = 0, lty = 3)
#   title(
#     title1,
#     cex.main = 1.2,
#     font.main = 1,
#     col.main = "purple",
#     ylab = NULL
#   )
#   plot(
#     c(1:length(toplot1_2)),
#     toplot1_2,
#     type = 'b',
#     ann = F,
#     ylim = range.y
#   )
#   do.text(2)
#   abline(h = 0, lty = 3)
#   plot(
#     c(1:length(toplot1_3)),
#     toplot1_3,
#     type = 'b',
#     ann = F,
#     ylim = range.y
#   )
#   do.text(3)
#   abline(h = 0, lty = 3)
#   return(Tk)
# }

# not used in EigenMS
make.formula.string <- function(factors, do.interactions=FALSE){
  fs <- "1"
  if(length(factors)){
    fs <- paste(factors, collapse=" + ")
    if(do.interactions && length(factors) > 1)
      fs <- paste(unlist(lapply(as.data.frame(t(gtools::combinations(length(factors), 2, factors)), stringsAsFactors=FALSE), paste, collapse="*")), collapse = " + ")
  }
  return(fs)
}


# make a string formula to use in 'lm' call when computing grp differences to preserve
makeLMFormula <- function(eff, var_name='') {
  # eff - effects used in contrasts
  # var_name - for singe factor use var-name that is passed in as variable names, otherwise it has no colnmae
  #           only used for a single factor
  if(is.factor(eff)){
    ndims <- 1
    cols1 <- var_name # ftemp in EigenMS
  }
  else
  {
    ndims <- dim(eff)[2]
    cols1 <- colnames(eff)
  }
  lhs <- cols1[1]
  lm.fm <- NULL
  # check if can have a list if only have 1 factor...

  params <- paste('contrasts=list(', cols1[1], '=contr.sum', sep=)

  if (ndims > 1) { # removed ndims[2] here, now ndims holds only 1 dimention...
    for (ii in 2:length(cols1))
    {
      lhs <- paste(lhs, "+", cols1[ii])  # bl="contr.sum",
	    params <- paste(params, ',', cols1[ii], '=contr.sum', sep='')
	  }
  }
  params <- paste(params,")")
  lm.formula <- stats::as.formula(paste('~', lhs))
  lm.fm$lm.formula <- lm.formula
  lm.fm$lm.params <- params
  return(lm.fm)
}


# First portion of EigenMS: identify bias trends
eig_norm1 <- function(m, treatment, prot.info, write_to_file=''){
# Identify significant eigentrends, allow the user to adjust the number (with causion! if desired)
# before normalizing with eig_norm2
#
# Input:
#   m: An m x n (peptides x samples) matrix of expression data, log-transformed!
#      peptide and protein identifiers come from the get.ProtInfo()
#   treatment:  either a single factor indicating the treatment group of each sample i.e. [1 1 1 1 2 2 2 2...]
#               or a frame of factors:  treatment= data.frame(cbind(data.frame(Group), data.frame(Time))
#   prot.info: 2+ colum data frame, pepID, prID columns IN THAT ORDER.
#              IMPORTANT: pepIDs must be unique identifiers and will be used as Row Names
#              If normalizing non-proteomics data, create a column such as: paste('ID_',seq(1:num_rows), sep='')
#              Same can be dome for ProtIDs, these are not used for normalization but are kept for future analyses
#   write_to_file='' - if a string is passed in, 'complete' peptides (peptides with NO missing observations)
#              will be written to that file name
#
# Output: list of:
#   m, treatment, prot.info, grp - initial parameters returned for futre reference
#   my.svd - matrices produced by SVD
#   pres - matrix of peptides that can be normalized, i.e. have enough observations for ANOVA,
#   n.treatment - number of factors passed in
#   n.u.treatment - number of unique treatment facotr combinations, eg:
#                   Factor A: a a a a c c c c
#                   Factor B: 1 1 2 2 1 1 2 2
#                   then:  n.treatment = 2; n.u.treatment = 4
#   h.c - bias trends
#   present - names/IDs of peptides on pres
#   complete - complete peptides, no missing values, these were used to compute SVD
#   toplot1 - trends automatically produced, if one wanted to plot at later time.
#   Tk - scores for each bias trend
#   ncompl - number of complete peptides with no missing observations
  print("Data dimentions: ")
  print(dim(m))
  # check if treatment is a 'factor' vs data.frame', i.e. single vs multiple factors
  if(methods::is(treatment, "factor")) { # TRUE if one factor
     n.treatment <- 1 # length(treatment)
     n.u.treatment <- length(unique(treatment))[1]
  } else { # data.frame
    n.treatment <- dim(treatment)[2]
    n.u.treatment <- dim(unique(treatment))[1] # all possible tretment combinations
  }
  # convert m to a matrix from data.frame
  m <- as.matrix(m) # no loss of information

  # filter out min.missing, here just counting missing values
  # if 1+ treatment completely missing, cannot do ANOVA, thus cannot preserve grp diff.
  # IMPORTANT: we create a composite grp = number of unique combinations of all groups, only for
  # 'nested' groups for single layer group is left as it is
  grpFactors <- treatment # temporary var, leftover from old times...

  nGrpFactors <- n.treatment # length(colnames(treatment)) # not good: dim(grpFactors)
  if(nGrpFactors > 1) { # got nested factors
    ugrps <- unique(grpFactors)
    udims <- dim(ugrps)
    grp <- NULL
    for(ii in seq_len(udims[1])) {
      pos <- grpFactors[,1] == ugrps[ii,1] # set to initial value
      for(jj in seq(2,udims[2])) {
	    pos <- pos & grpFactors[,jj] == ugrps[ii,jj]
      }
	  grp[pos] <- rep(ii, sum(pos))
    }
	grp <- as.factor(grp)
  } else {
    grp <- treatment
  }
  nobs <- array(NA, c(nrow(m), length(unique(grp)))) # noobs = number of observations

  print('Treatmenet groups:')
  print(grp)

  for(ii in seq_len(nrow(m))) {
    for(jj in seq_len(length(unique(grp)))) {
      nobs[ii,jj] <- sum(!is.na(m[ii, grp==unique(grp)[jj]])) # total number of groups num(g1) * num(g2) * ...
    }
  }
  # now 'remove' peptides with missing groups
  present.min <- apply(nobs, 1, min) # number present in each group
  ii <- present.min == 0   # 1+ obs present in ALL of the groups
  nmiss <- sum(present.min == 0) # not used, one value of how many peptides have 1+ grp missing completely
  pmiss <- rbind(m[ii,]) # these have 1+ grp missing !!!!
  # rownames must be UNIQUE, if have possible duplicates: use 'ii' ?
  rownames(pmiss) <- prot.info[ii,1]  # set rownames,

  # create matrix for peptides with enough observations for ANOVA
  # 'present' are names of the peptides (pepID) and 'pres' are abundances
  # NOTE: ! negates the proteins, so we get ones that have 1+ obs in each group
  present <- prot.info[which(!prot.info[,1] %in% rownames(pmiss)), ] # rownames OK
  # pres = m[which(!rownames(m) %in% rownames(pmiss)), ]
  pres <- m[which(!prot.info[,1] %in% rownames(pmiss)), ] # is this OK?
  rownames(pres) <- prot.info[which(!prot.info[,1] %in% rownames(pmiss)),1]

  print('Selecting complete peptides')
  # Should issue an error message if we have NO complete peptides.
  # select only 'complete' peptides, no missing values
  nobs <- array(NA, nrow(pres)) # reassign noobs to dims of 'present'
  numiter <- nrow(pres)
  for (ii in seq_len(numiter)) {
    # if(ii %% 100 == 0) { print(ii) }
    nobs[ii] <- sum(!is.na(pres[ii,]))
  }

  iii <- nobs == ncol(pres)
  complete <- rbind(pres[iii,])

  #  write out a file of complete peptides if file name is passed in
  if(write_to_file != '') {
    utils::write.table(complete, file = write_to_file, append = FALSE,
             quote = FALSE, sep = "\t",
             eol = "\n", na = "NaN", dec = ".", row.names = TRUE,
             col.names = TRUE, qmethod = c("escape", "double"))
  }

  # compute bias with 'complete' matrix and residuals from 'present'
  # calculate eigenpeptides for 'complete' data only
  # if have only 1 group, we do not need to preserve group differernces, everything is the same group, ex: QC samples
  # contrasts will fail if have only 1 group, thus have else
  if(n.u.treatment > 1) {
    print('Got 2+ treatment grps')
    # check to see if we have multiple factors
	  grpdim <- dim(treatment)

	  lm.fm <- makeLMFormula(treatment, 'TREAT') # using general function that can accomodate for 1+ number of factors
    TREAT <- treatment
    TREAT <- data.frame(treatment) # temp var to work if we got only 1 treatment vector.
    if(methods::is(treatment,"factor")) {
      colnames(TREAT) <- "TREAT"
     } else {
      colnames(TREAT) <- colnames(treatment)
    }
    attach(TREAT)

    mod.c <- stats::model.matrix(lm.fm$lm.formula, data=TREAT, eval(parse(text=lm.fm$lm.params)))
	  Y.c <- as.matrix(complete)
	  options(warn = -1)

    # use lm() to get residuals
    formula1 <- paste('t(Y.c)~', as.character(lm.fm$lm.formula)[2], sep = '')
    TREAT <- treatment
    fit_lmAll <- stats::lm(eval(parse(text=formula1)))
    R.c <- stats::residuals(fit_lmAll)  # Oct 2 messing with residuals...
  } else {  # 1 group only, set residuals to original matrix
    print('Got 1 treatment grp')
	  mod.c <- as.numeric(t(treatment))
    R.c <- t(as.matrix(complete))  # needs to be transposed to match the matrix returned from lm
    TREAT <- treatment
  }

  print('Computing SVD, estimating Eigentrends...') # let user know what is going on
  # residuals are centered around 0, here center samples not peptides/metabolites
  # centering is basic normalization

  R.c_center <- scale(R.c, center = TRUE, scale = FALSE)  # t(scale(t(R.c), center = TRUE, scale = FALSE))
  my.svd <- svd(R.c_center)  # can use wrapper below to chek if SVD has a problem...
  temp <- my.svd$u
  my.svd$u <- my.svd$v
  my.svd$v <- temp

  #identify number of eigenvalues that account for a significant amount of residual variation
  numcompletepep <- dim(complete)[1] # save to return to the user as part of the return list
  # this is important info for publications
  # tell users how many peptides/metabolites the trends are based on
  # can also be determined by doing dim(return_value_fromEIg_norm1$pres)

  print(paste('Number of treatments: ', n.u.treatment))
  h.c <- sva.id(complete, treatment, n.u.treatment, lm.fm=lm.fm)$n.sv
  print(paste("Number of significant eigenpeptides/trends", h.c) )

  # show RAW trends
  # center each peptide around zero (subtract its mean across samples)
  complete_center <- scale(t(complete), center = TRUE, scale = FALSE)
  print('Preparing to plot...')

  n.u.treatment
  toplot1 <- svd(complete_center) # scales above
  temp <- toplot1$u
  toplot1$u <- toplot1$v
  toplot1$v <- temp

  graphics::par(mfcol=c(3,2))
  graphics::par(mar = c(2,2,2,2))
  #plot.eigentrends(toplot1, "Raw Data")
  #plot.eigentrends(my.svd, "Residual Data")

  d <- my.svd$d;  ss <- d^2;
  Tk <- signif(ss/sum(ss)* 100, 2)

  retval <- list(m=m, treatment=treatment, my.svd=my.svd,
            pres=pres, n.treatment=n.treatment, n.u.treatment=n.u.treatment,
            h.c=h.c, present=present, prot.info=prot.info,
            complete=complete, toplot1=toplot1, Tk=Tk, ncompl=numcompletepep,
            grp=grp)
  return(retval)
}


# Second portion of EigenMS: remove effects of bias
# split into 2 functions to allo the user to examine the bias trends and be able to change the number
# by resetting h.c returned by eig_norm1(return_value_from_eig_norm1)
eig_norm2 <- function(rv) {
  # UNPUT:
  #   rv - return value from the eig_norm1
  #   if user wants to change the number of bias trends that will be eliminated h.c in rv should
  #   be updates to the desired number
  #
  # OUTPUT:
  #   normalized - matrix of normalized abundances with 2 columns of protein and peptdie names
  #   norm_m - matrix of normalized abundances, no extra columns
  #   eigentrends - found in raw data, bias trendsup to h.c
  #   rescrange - rescaling range for the addition of the while noise to avoid overfitting
  #   norm.svd - trends in normalized data, if one wanted to plot at later time.
  #   exPeps - excluded peptides - excluded due to exception in fitting a linear model

  m <- rv$pres # yuliya: use pres matrix, as we cannot deal with m anyways, need to narrow it down to 'complete' peptides
  treatment <- rv$treatment
  my.svd <- rv$my.svd
  pres <- rv$pres
  n.treatment <- rv$n.treatment
  n.u.treatment <- rv$n.u.treatment
  numFact <- dim(rv$treatment)[2]
  print(paste('Unique number of treatment combinations:', n.u.treatment) )
  h.c <- rv$h.c
  present <- rv$present
  toplot1 <- rv$toplot1
  # vector of indicators of peptides that threw exeptions
  exPeps <- vector(mode = "numeric", length = nrow(pres))

  print("Normalizing...")
  treatment <- data.frame(treatment) # does this need to be done?
  if(n.u.treatment > 1) {
    lm.fm <- makeLMFormula(treatment, 'ftemp')
    mtmp <- stats::model.matrix(lm.fm$lm.formula, data=treatment, eval(parse(text=lm.fm$lm.params)))  #contrasts=list(bl="contr.sum", it="contr.sum",Pi="contr.sum", tp="contr.sum"))
  } else {  # have 1 treatment group
    mtmp <- treatment # as.numeric(t(treatment))
  }
  # above needed to know how many values will get back for some matrices
  # create some variables:
  betahat <- matrix(NA,nrow=dim(mtmp)[2],ncol=nrow(pres))
  newR <- array(NA, c(nrow(pres), ncol(pres))) #, n.treatment))
  norm_m <- array(NA, c(nrow(pres), ncol(pres))) # , n.treatment))
  numsamp <- dim(pres)[2]
  numpep <- dim(pres)[1]
  betahat_n <- matrix(NA,nrow=dim(mtmp)[2],ncol=nrow(pres))
  rm(mtmp)

  V0 <- my.svd$v[,seq_len(h.c),drop=FALSE]   # residual eigenpeptides

  if(n.u.treatment == 1) { # got 1 treatment group
    for (ii in seq_len(nrow(pres))) {
      if(ii%%250 == 0) { print(paste('Processing peptide ',ii))  }
      pep <- pres[ii, ]
      pos <- !is.na(pep)
      peptemp <- as.matrix(pep[pos]) # take only the observed values
      resm <- rep(NA, numsamp)
      resm[pos] <- as.numeric(pep[pos])
      bias <- array(NA, numsamp)
      bias[pos] <- resm[pos] %*% V0[pos,] %*% t(V0[pos,])
      norm_m[ii, ] <- as.numeric(pep - bias)
    }

  } else { # got 2+ treatment groups
    for (ii in seq_len(nrow(pres))) {
      if(ii %% 100 == 0) { print(paste('Processing peptide ',ii))  }
      pep <- pres[ii, ]
      pos <- !is.na(pep)
      peptemp <- as.matrix(pep[pos]) # take only the observed values, may not be needed in R? but this works
      ftemp <- treatment[pos,]
      ftemp <- data.frame(ftemp)
      #### use try, not entirely sure if need for modt, need it for solve lm?!
      options(warn = -1)
      lm.fm <- makeLMFormula(ftemp, 'ftemp') # using general function that can accomodate for 1+ number of factors
      modt <- try(stats::model.matrix(lm.fm$lm.formula, data=ftemp, eval(parse(text=lm.fm$lm.params))), silent=TRUE)
      options(warn = 0)

      if(!inherits(modt, "try-error")) { # do nothing if could not make model matrix
        options(warn = -1)
        # if we are able to solve this, we are able to estimate bias
        bhat <-  try(solve(t(modt) %*% modt) %*% t(modt) %*% peptemp)
        options(warn = 0)
        if(!inherits(bhat, "try-error")) {
          betahat[,ii] <- bhat
          ceffects <- modt %*% bhat  # these are the group effects, from estimated coefficients betahat

          resm <- rep(NA, numsamp) # really a vector only, not m
          resm[pos] <- as.numeric(pep[pos] - ceffects)
          bias <- array(NA, numsamp)
          bias[pos] <- resm[pos] %*% V0[pos,] %*% t(V0[pos,])
          norm_m[ii, ] <- as.numeric(pep - bias)

          # yuliya:  but newR should be computed on Normalized data
          resm_n <- rep(NA, numsamp)
          bhat_n <-  solve(t(modt) %*% modt) %*% t(modt) %*% norm_m[ii, pos]
          betahat_n[,ii] <- bhat_n
          ceffects_n <- modt %*% bhat_n
          resm_n[pos] <- norm_m[ii,pos] - ceffects
          newR[ii, ] <- resm_n
        } else {
          print(paste('got exception 2 at peptide:', ii, 'should not get here...'))
          exPeps[ii] <- 2 # should not get 2 here ever...
        }
      } else {
        print(paste('got exception at peptide:', ii))
        exPeps[ii] <- 1 # keep track of peptides that threw exeptions, check why...
      }
    }
  } # end else - got 2+ treatment groups

  #####################################################################################
  # rescaling has been eliminated form the code after discussion that bias
  # adds variation and we remove it, so no need to rescale after as we removed what was introduced
  y_rescaled <- norm_m # for 1 group normalization only, we do not rescale
  # add column names to y-rescaled, now X1, X2,...
  colnames(y_rescaled) <- colnames(pres) # these have same number of cols
  rownames(y_rescaled) <- rownames(pres)
  y_resc <- data.frame(present, y_rescaled)
  rownames(y_resc) <- rownames(pres)  # rownames(rv$normalized)
  final <- y_resc # row names are assumed to be UNIQUE, peptide IDs are unique

  # rows with all observations present
  complete_all <- y_rescaled[rowSums(is.na(y_rescaled))==0,,drop=FALSE]

  #  x11() # make R open new figure window
  graphics::par(mfcol=c(3,2))
  graphics::par(mar = c(2,2,2,2))
  # center each peptide around zero (subtract its mean across samples)
  # note: we are not changing matrix itself, only centerig what we pass to svd
  complete_all_center <- t(scale(t(complete_all), center = TRUE, scale = FALSE))
  toplot3 <- svd(complete_all_center)
  #plot.eigentrends(toplot1, "Raw Data")
  #plot.eigentrends(toplot3, "Normalized Data")

  print("Done with normalization!!!")
  colnames(V0) <-  paste("Trend", seq_len(ncol(V0)), sep="_")

  maxrange <- NULL # no rescaling # data.matrix(maxrange)
  return(list(normalized=final, norm_m=y_rescaled, eigentrends=V0, rescrange=maxrange,
              norm.svd=toplot3, exPeps=exPeps))
} # end function eig_norm2



### EigenMS helper functions, a few more...
# Tom had Sig set to 0.1, I and Storey's paper has 0.05
# sva.id = function(dat, mod, n.u.treatment, B=500, sv.sig=0.05, seed=NULL) {

# yuliya Sept 1, 2014: changed parameter mod to be treatment, using lm() to obtain residuals
sva.id <- function(dat, treatment, n.u.treatment, lm.fm, B=500, sv.sig=0.05) {
# Executes Surrogate Variable Analysis
# Input:
#   dat: A m peptides/genes by n samples matrix of expression data
#   mod: A model matrix for the terms included in the analysis
#   n.u.treatment - 0 or 1, if we are normalizing data with NO groups or some groups, QC vs samples
#   B: The number of null iterations to perform
#   sv.sig: The significance cutoff for the surrogate variables
#   seed: A seed value for reproducible results
# Output
#    n.sv: Number of significant surrogate variables.
#    id: An indicator of the significant surrogate variables
#    B: number of permutation to do
#    sv.sig: significance level for surrogate variables
  print("Number of complete peptides (and samples) used in SVD")
  print(dim(dat))

  warn <- NULL
  n <- ncol(dat)
  m <- nrow(dat)

  # ncomp = length(as.numeric(n.u.treatment))
  ncomp <- n.u.treatment # JULY 2013: as.numeric(n.u.treatment)
  print(paste("Number of treatment groups (in svd.id): ", ncomp))
  # should be true for either case and can be used later

  if(ncomp > 1) { #
    formula1 <- paste('t(dat)~', as.character(lm.fm$lm.formula)[2], sep = '')
    fit_lmAll <- stats::lm(eval(parse(text=formula1)))
    res <- t(stats::residuals(fit_lmAll))
  } else {
    res <- dat
  }
  # centering was not done before...
  # center each peptide around zero (subtract its mean across samples)
  # note: we are not changing matrix itself, only centerig what we pass to svd
  res_center <- t(scale(t(res), center = TRUE, scale = FALSE))

  uu <- svd(t(res_center)) # NEED a WRAPPER for t(). the diag is min(n, m)
  temp <- uu$u
  uu$u <- uu$v
  uu$v <- temp


  # yuliya: Sept 2014: can I get around without using H??
  #  ndf = min(n, m) - ceiling(sum(diag(H)))
  #  dstat = uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
  #  dstat0 = matrix(0,nrow=B,ncol=ndf)
  #  s0 = diag(uu$d) # no need for diag here, should investigate why this is a vector already...
  s0 <- uu$d
  s0 <- s0^2
  dstat <- s0/sum(s0)  # this did not have 'diag' in it in Tom's code...
  ndf <- length(dstat) # sticking to Tom's variable name
# print(paste('length(dstat) = ndf = ', ndf))
  dstat0 <- matrix(0,nrow=B,ncol=ndf) # num samples (?)

print("Starting Bootstrap.....")
# this is the Bootstrap procedure that determines the number of significant eigertrends...
  for(ii in seq_len(B)){
    if(ii %% 50 == 0) { print(paste('Iteration ', ii)) }
    res0 <- t(apply(res, 1, sample, replace=FALSE)) # regression
    # yuliya: not sure if this is needed at all
    # not needed for 1 group normalizaiton
    ##### res0 = res0 - t(H %*% t(res0))
    # yuliya: Sept 3, 2014: REMOVED above line. Do not think this needs to be done..
    # center each peptide around zero (subtract its mean across samples)
    # note: we are not changing matrix itself, only centerig what we pass to svd
    res0_center <- t(scale(t(res0), center = TRUE, scale = FALSE))
    uu0 <- svd(res0_center)
  	temp <- uu0$u  # why did tom do this??
    uu0$u <- uu0$v
    uu0$v <- temp

    ss0 <- uu0$d  # no need for diag....
    ss0 <- ss0^2
    dstat0[ii,] <- ss0 / sum(ss0) # Tk0 in Matlab
  }

# yuliya: check p-values here, Tom had mean value...
  psv <- rep(1,n)
  for(ii in seq_len(ndf)){
	  # psv[ii] = mean(dstat0[,ii] >= dstat[ii])
    # should this be compared to a MEAN?  Should this be dstat0[ii,] ?
    posGreater <- dstat0[,ii] > dstat[ii]
    psv[ii] <- sum(posGreater) / B
  }

  # p-values for peptides have to be in monotonically increasing order,
  # set equal to previous one if not the case
  for(ii in 2:ndf){
    if(psv[(ii-1)] > psv[ii]) {
	    # psv[ii] = max(psv[(ii-1)],psv[ii])
      psv[ii] <- psv[(ii-1)]
    }
  }
  nsv <- sum(psv <= sv.sig)
  # tom - this should be at least 1
  # nsv = min(sum(psv <= sv.sig), 1, na.rm=T)
  return(list(n.sv = nsv,p.sv=psv))
}
# end sva.id


mmul <- function(A, B){
  # multiply square matrix by rectangle with NA's (???)
  X <- A
  Y <- B
  X[is.na(A)] <- Y[is.na(B)] <- 0
  R <- X %*% Y
  R[is.na(B)] <- NA
  R
}
# end mmul



#####    Not used in EigenMS, but a solution to imputing peptides with 1+ group missing completely.
#######################################################################
imp_elim_proteins <- function(nosib.miss, treatment){
# Impute 1 peptide proteins with 1 grp completely missing
#
# INPUT: proteins - 1 peptide proteins with 1 grp missing completely (NaNs)
#        grps - grouping information for observations, can be 2+ grps here,
#               but not sure what to do with more then 2 groups if 2+ gs are
#               missing completely. So only 1 grp is completely missing.
# OUTPUT: prs - imputed proteins, same size as proteins parameter
#nosib.miss <= nosib.miss
#treatment <= treatment
  tlvls <- unique(treatment)
  proteins <- nosib.miss
  ng <- length(unique(treatment))
  for (i in seq_len(nrow(proteins))) {
  	pr <- as.vector(proteins[i,])
    #% find number of missing values in ALL groups
    miss <- array(NA, c(1, ng))
    for (j in seq_len(ng)) miss[j] <- sum(is.na(pr[treatment==unique(treatment)[j]]))
#    pos = miss==0
    pos <- miss==min(miss)  # Tom 092510
    present_groups <- pr[treatment==unique(treatment)[pos]]

    # compute mean and stdev from one of the present groups, make sure no NaNs are used
	pospres <- !is.na(present_groups)
	presvals <- present_groups[pospres]
	pepmean <- mean(presvals)
	pepstd <-  stats::sd(presvals)
	if(is.na(pepstd)) next;

  #% imputing only COMPLETELY missing peptides here
	for (j in seq_len(ng)) {
		if   (!pos[j]) { #% imute only the ones not at pos complete
			imppos <- is.na(pr[treatment==tlvls[j]])  #% should be all in this group, but double check
			imppepmean <- pepmean - 6* pepstd
			imppepstd <- pepstd
			tt <- imppepmean - 3 * imppepstd
			kk <- 0
      while (tt < 0 && kk < 10){  # added kk counter - tom 092510
				offset <- .25
				gap <- imppepstd * offset
				imppepstd <- imppepstd * (1- offset)
				imppepmean <- imppepmean + 3 * gap
				tt <- imppepmean - 3 * imppepstd
				kk <- kk + 1
			}
			imp_tmp <- stats::rnorm(length(imppos), imppepmean, pepstd)
			pr[treatment==tlvls[j]] <- imp_tmp
		}
	}
  proteins[i,] <- pr
  }

  # tom - this routine gives some nearly-blank rows
  # to avoid singularities, i'm going to scan this to see which protein rows
  # have all blanks in one row and remove them
  #  proteins <= proteins
  xx <- (!is.na(proteins)) %*% stats::model.matrix(~treatment-1)
  notblank.idx <- rep(TRUE, nrow(xx))
  for(jj in seq_len(ncol(xx))) notblank.idx <- notblank.idx & xx[,jj]
#  proteins = proteins[blank.idx,,drop=FALSE]

  return(list(proteins=proteins, notblank.idx=notblank.idx))
}

#######################################################################
#######################################################################
my.Psi <- function(x, my.pi){
# calculates Psi
exp(log(1-my.pi)  + stats::dnorm(x, 0, 1, log=TRUE) - log(my.pi + (1 - my.pi) * stats::pnorm(x, 0, 1) ))
}
# end my.Psi

my.Psi.dash <- function(x, my.pi){
# calculates the derivative of Psi
-my.Psi(x, my.pi) * (x + my.Psi(x, my.pi))
}
# end my.Psi.dash

phi <- function(x){stats::dnorm(x)}

rnorm.trunc <- function (n, mu, sigma, lo=-Inf, hi=Inf){
# Calculates truncated noraml
  p.lo <- stats::pnorm (lo, mu, sigma)
  p.hi <- stats::pnorm (hi, mu, sigma)
  u <- stats::runif (n, p.lo, p.hi)
  return (stats::qnorm (u, mu, sigma))
}
# end rnorm.truncf
