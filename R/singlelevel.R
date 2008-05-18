`singlelevel` <-
function(isold,sold,covs={},lk="logit",int=FALSE,modify=TRUE){
  stopifnot(length(isold)==length(sold))
  stopifnot(length(unique(sold))==2)
  lcovs <- max(length(covs),dim(covs)[1])
  if (is.null(covs)==FALSE) stopifnot(length(sold)==lcovs)
  if (max(sold)>1 || min(sold)<0){
   print("Your values for sold were re-scaled to 1 and 0")
   sold <- as.numeric(sold == max(sold))}
  if (modify)
    isold <- as.numeric(isold == max(isold)) - .5
  mod1 <- glm(sold ~ isold,family=binomial(link=lk))
  if (lk=="probit"){
   print("Overall estimated d' ",quote=FALSE)
   print(summary(mod1)$coef[2,])
   }
  if (lk=="logit"){
   print("Overall estimated lnOR ",quote=FALSE)
   print(summary(mod1)$coef[2,])
   }
  if (is.null(covs)) return(mod1)
  ncovs <- max(1, dim(covs)[2])
  for (i in 1:ncovs){
    if (ncovs == 1) cov <- covs
    if (ncovs > 1) cov <- covs[,i]
    if (is.factor(cov) || is.ordered(cov)){
       mod2 <- glm(sold ~ -1 + cov + isold:cov, family=binomial(link=lk))
      if (lk == "probit")
        print("The estimated d' values for the groups are:",quote=FALSE)
      if (lk == "logit")
        print("The estimated lnOR values for the groups are:",quote=FALSE)
      vals <- length(unique(cov))
      print(summary(mod2)$coef[(vals+1):(2*vals),])
      m1 <- glm(sold ~ cov + isold, family=binomial(link=lk))
      m2 <- update(m1, .~. + isold:cov)
      print(anova(m1,m2,test = "Chisq"))
    }
    if (is.numeric(cov)){
      varc <- (cov-mean(cov))/sd(cov)
      mod2 <- glm(sold ~ isold*varc, family=binomial(link=lk))
      if (lk == "probit")
         print("The estimated shift in d' for a standard deviation shift is:",quote=FALSE)
      if (lk == "logit")
         print("The estimated shift in lnOR for a standard deviation shift is:",quote=FALSE)
      print(summary(mod2)$coef[4,])
      m1 <- glm(sold ~ isold+varc, family=binomial(link=lk))
      print(anova(m1,mod2,test = "Chisq"))
    }
    }
if (ncovs == 1)
 fmint <- (as.formula(paste("sold ~ isold * covs")))
if (ncovs > 1 && int==FALSE)
 fmint <- (as.formula(paste("sold ~ isold * (",paste(dimnames(covs)[[2]],collapse="+"),paste(")"))))
if (ncovs > 1 && int){
 indeps <- cbind(isold,covs)
 fmint <- (as.formula(paste("sold ~ ", paste(dimnames(indeps)[[2]],collapse="*"))))}
final <- glm(fmint,family=binomial(link=lk))
return(final)}

