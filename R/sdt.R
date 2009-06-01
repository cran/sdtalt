`sdt` <-
function (hits,fas,misses,cr,flat=0,pmeans=FALSE,meas="all",wk=.5,runboot=FALSE,trim=0,R=2000,confl=.95,newst=NULL,bound=FALSE,...)
 {stopifnot (length(hits)==length(fas),length(fas)==length(misses),length(misses)==length(cr))
  statsfile <- {}
  hits <- hits+flat; misses <- misses+flat; cr <- cr + flat; fas <- fas + flat
  HR <- hits/(hits+misses)
  FAR <- (fas)/(fas+cr)
  if (any(meas=="all" | meas=="HR"))
    statsfile <- cbind(statsfile,HR)
  if (any(meas=="all" | meas=="FAR"))
    statsfile <- cbind(statsfile,FAR)
  if (any(meas=="all" | meas=="d")){
    d <- qnorm(HR) - qnorm(FAR)
    if (bound)
     {d[is.infinite(d)&d>0] <- max(d[is.finite(d)])
      d[is.infinite(d)&d<0] <- min(d[is.finite(d)])
     }
    statsfile <- cbind(statsfile,d)}
  if (any(meas=="all" | meas=="csdt")){
  csdt <- -.5*(qnorm(HR) + qnorm(FAR))
    if (bound)
     {csdt[is.infinite(csdt)&csdt>0] <- max(csdt[is.finite(csdt)])
      csdt[is.infinite(csdt)&csdt<0] <- min(csdt[is.finite(csdt)])
     }
    statsfile <- cbind(statsfile,csdt)}
  if (any(meas=="all" | meas=="A")){
  A <- .5 + (sign(HR-FAR)*((HR-FAR)^2 + abs(HR-FAR))/(4*pmax(HR,FAR)-4*HR*FAR))
    statsfile <- cbind(statsfile,A)}
  if (any(meas=="all" | meas=="B")){
  B <- sign(HR-FAR)*(HR*(1-HR)-FAR*(1-FAR))/(HR*(1-HR)+FAR*(1-FAR))
    statsfile <- cbind(statsfile,B)}
  if (any(meas=="all" | meas=="lnbeta")){
    lnbeta <- (qnorm(FAR)^2 - qnorm(HR)^2)/2
    if (bound)
     {lnbeta[is.infinite(lnbeta)&lnbeta>0] <- max(lnbeta[is.finite(lnbeta)])
      lnbeta[is.infinite(lnbeta)&lnbeta<0] <- min(lnbeta[is.finite(lnbeta)])
     }
    statsfile <- cbind(statsfile,lnbeta)}
  if (any(meas=="all" | meas=="beta")){
    beta <- exp((qnorm(FAR)^2 - qnorm(HR)^2)/2)
    if (bound)
     {beta[is.infinite(beta)&beta>0] <- max(beta[is.finite(beta)])
      beta[is.infinite(beta)&beta<0] <- min(beta[is.finite(beta)])
     }
    statsfile <- cbind(statsfile,beta)}
  if (any(meas=="all" | meas=="OR")){
    OR <- (hits*cr)/(misses*fas)
    if (bound)
     {OR[is.infinite(OR)&OR>0] <- max(OR[is.finite(OR)])
      OR[is.infinite(OR)&OR<0] <- min(OR[is.finite(OR)])
     }
    statsfile <- cbind(statsfile,OR)}
  if (any(meas=="all" | meas=="lnOR")){
    lnOR <- log((hits*cr)/(misses*fas))
    if (bound)
     {lnOR[is.infinite(lnOR)&lnOR>0] <- max(lnOR[is.finite(lnOR)])
      lnOR[is.infinite(lnOR)&lnOR<0] <- min(lnOR[is.finite(lnOR)])
     }
    statsfile <- cbind(statsfile,lnOR)}
    temp1 <- hits*cr-misses*fas
    temp2 <- (hits+fas)*(fas+cr)*wk + (misses+cr)*(hits+misses)*(1-wk)
  if (any(meas=="all" | meas=="kappa")){
    kappa <- temp1/temp2
    statsfile <- cbind(statsfile,kappa)}
  if (any(meas=="all" | meas=="phi")){
    phi <- temp1/((hits+fas)*(fas+cr)*(misses+cr)*(hits+misses))^.5
    statsfile <- cbind(statsfile,phi)}
  if (any(meas=="all" | meas=="Q")){
    Q <- temp1/(hits*cr+misses*fas)
    statsfile <- cbind(statsfile,Q)}
  if (any(meas=="all" | meas=="eta")){
    eta <- sqrt((fas*misses)/(hits*cr))
    if (bound)
     {eta[is.infinite(eta)&eta>0] <- max(eta[is.finite(eta)])
      eta[is.infinite(eta)&eta<0] <- min(eta[is.finite(eta)])
     }
    statsfile <- cbind(statsfile,eta)}
  if (any(meas=="all" | meas=="PC")){
    PC <- (hits+cr)/(hits+misses+cr+fas)
    statsfile <- cbind(statsfile,PC)}
  numbns <- length(newst); counter <- numbns
  while (counter > 0) {
    func <- newst[[numbns - counter + 1]]
    newstat <- func(hits,fas,misses,cr)
    if (bound)
     {newstat[is.infinite(newstat)&newstat>0] <- max(newstat[is.finite(newstat)])
      newstat[is.infinite(newstat)&newstat<0] <- min(newstat[is.finite(newstat)])
     }
    statsfile <- cbind(statsfile,newstat)
    counter <- counter - 1}
  if (pmeans) {
    meantr <- function(x) mean(x,tr=trim)
    print("The means:",quote=FALSE)
    print(apply(statsfile,2,meantr),digits=3)
    print("The standard deviations:",quote=FALSE)
    print(sd(statsfile),digits=3)
    if (runboot){
      bootatt <- search()
      if (length(bootatt[bootatt=="package:boot"])<1) {
        bootinlib <- installed.packages()
        if (length(bootinlib[bootinlib[,1]=="boot"])<2){
          print("Need to install boot (Canty, 2007). Choose a mirror", quote=FALSE)
    install.packages("boot")}
    library("boot")}
      for (j in 1:dim(statsfile)[2]){
        runs <- boot(statsfile[,j],function(x,i) mean(x[i]),R)
        print(paste("**** BCa intervals for",colnames(statsfile)[j],"****"),quote=FALSE)
        bootout <- boot.ci(runs,conf=confl,type="bca")$bca
        colnames(bootout) <- c("Conf"," "," ","lb","ub")
        print(bootout[1,c(1,4,5)],digits=4)}}}
return(as.data.frame(statsfile))}

