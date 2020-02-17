cv.sdwd = function(x, y, lambda=NULL, 
  pred.loss=c("misclass", "loss"), nfolds=5, foldid, ...) {
  ####################################################################
  ## data setup
  if (missing(pred.loss)) 
    pred.loss = "misclass" else pred.loss = match.arg(pred.loss)
  if (!match(pred.loss, c("misclass", "loss"), FALSE)) {
    warning("Only 'misclass' and 'loss' available for 
            DWD; 'misclass' used.")
    pred.loss = "misclass"
  }
  typenames = c(misclass = "mis-classification error",
                loss     = "DWD loss")
  y = c(-1, 1)[as.factor(drop(y))]
  x = as.matrix(x)
  if (!all(y %in% c(-1, 1))) 
    stop("y should be a factor with two levels.")
  x.row = as.integer(NROW(x))
  if (length(y) != x.row) 
    stop("x and y have different number of observations.")  
  ####################################################################
  sdwd.object = sdwd(x, y, lambda=lambda, ...)
  lambda = sdwd.object$lambda
  nz = sapply(coef(sdwd.object, type="nonzero"), length)
  if (missing(foldid)) 
    foldid = sample(rep(seq(nfolds), 
      length=x.row)) else nfolds = max(foldid)
  if (nfolds < 3) 
    stop("nfolds must be bigger than 3; nfolds=5 recommended.")
  outlist = as.list(seq(nfolds))
  ## fit the model nfold times and save them
  for (i in seq(nfolds)) {
    which = foldid == i
    outlist[[i]] = sdwd(x=x[!which, , drop=FALSE], 
                        y=y[!which], lambda=lambda, ...)
  }
  ## select the lambda according to predmat
  fun = paste("cv", class(sdwd.object)[[2]], sep=".")
  cvstuff = do.call(fun, 
    list(outlist, lambda, x, y, foldid, pred.loss))
  cvm = cvstuff$cvm
  cvsd = cvstuff$cvsd
  cvname = cvstuff$name
  out = list(lambda=lambda, cvm=cvm, cvsd=cvsd, 
              cvupper=cvm+cvsd, cvlower=cvm - cvsd, nzero=nz, 
              name=cvname, sdwd.fit=sdwd.object)
  obj = c(out, as.list(getmin(lambda, cvm, cvsd)))
  class(obj) = "cv.sdwd"
  obj
} 

cv.sdwdNET = function(outlist, lambda, x, y, foldid, pred.loss) {
  ### Turn y into c(-1,1)
  y = as.factor(y)
  y = c(-1, 1)[as.numeric(y)]
  nfolds = max(foldid)
  predmat = matrix(NA, length(y), length(lambda))
  nlams = double(nfolds)
  for (i in seq(nfolds)) {
    which = foldid == i
    fitobj = outlist[[i]]
    preds = predict(fitobj, x[which, , drop=FALSE], type="link")
    nlami = length(outlist[[i]]$lambda)
    predmat[which, seq(nlami)] = preds
    nlams[i] = nlami
  }
  typenames = c(misclass = "mis-classification error",
              loss     = "DWD loss")
  cvraw = switch(pred.loss, 
    loss     = dwdloss(y * predmat),
    misclass = (y != ifelse(predmat > 0, 1, -1)))
  if (length(y)/nfolds >= 3) {
    cvob = cvcompute(cvraw, foldid, nlams)
    cvraw = cvob$cvraw
    cvn = cvob$N
  } else cvn = length(y) - colSums(is.na(predmat))    
  cvm = colMeans(cvraw, na.rm=TRUE)
  cvsd = sqrt(colMeans(scale(cvraw, cvm, FALSE)^2, 
    na.rm=TRUE)/(cvn - 1))
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
} 

predict.cv.sdwd = function(object, newx, s=c("lambda.1se", 
    "lambda.min"), ...) {
  if (is.numeric(s)) 
    lambda = s else if (is.character(s)) {
      s = match.arg(s)
      lambda = object[[s]]
    } else stop("Invalid form for s")
  predict(object$sdwd.fit, newx, s=lambda, ...)
} 

plot.cv.sdwd = function(x, sign.lambda=1, ...) {
  cvobj = x
  xlab = "log(Lambda)"
  if (sign.lambda < 0) 
    xlab = paste0("-", xlab)
  plot.args = list(x=sign.lambda * log(cvobj$lambda), 
    y=cvobj$cvm, type="n", xlab=xlab, ylab=cvobj$name,
    ylim=range(cvobj$cvupper, cvobj$cvlo))
  new.args = list(...)
  if (length(new.args)) 
    plot.args[names(new.args)] = new.args
  do.call("plot", plot.args)
  error.bars(sign.lambda * log(cvobj$lambda), cvobj$cvupper, 
    cvobj$cvlo, width=0.01, col="darkgrey")
  points(sign.lambda * log(cvobj$lambda), cvobj$cvm, 
    pch=20, col="red")
  axis(side=3, at=sign.lambda * log(cvobj$lambda), 
    labels=paste(cvobj$nz), tick=FALSE, line=0)
  abline(v=sign.lambda * log(cvobj$lambda.min), lty=3)
  abline(v=sign.lambda * log(cvobj$lambda.1se), lty=3)
  invisible()
} 

coef.cv.sdwd = function(object, 
    s=c("lambda.1se", "lambda.min"), ...) {
  if (is.numeric(s)) 
    lambda = s else if (is.character(s)) {
      s = match.arg(s)
      lambda = object[[s]]
    } else stop("Invalid form for s.")
  coef(object$sdwd.fit, s=lambda, ...)
} 
