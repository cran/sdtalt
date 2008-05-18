`mlmsdt` <-
function(subno,isold,sold,covs=NULL,lk="logit",vardiff=TRUE,modify=TRUE,int=FALSE,item=NULL,...)
{stopifnot(length(isold)==length(sold),length(subno)==length(sold))
 stopifnot(length(unique(sold))==2)
 lcovs <- max(length(covs),dim(covs)[1])
 if (is.null(covs)==FALSE) stopifnot(length(sold)==lcovs)
 if (max(sold)>1 || min(sold)<0){
   print(paste("Your values for", sold,"were re-scaled to 1 and 0"))
   sold <- as.numeric(sold == max(sold))}
 if (modify)
    isold <- as.numeric(isold == max(isold)) - .5
 if (sd(subno) == 0){
  print("  ",quote=FALSE)
  print("You appear to have only a single participant, ",quote=FALSE)
  print("so the function singlelevel will be used.",quote=FALSE)
  print("  ",quote=FALSE)
  return(singlelevel(isold,sold,covs,lk=lk,...))}
 # The following used when function downloaded with source function
 #l4att <- search()
 #if (length(l4att[l4att=="package:lme4"])<1) {
 #  l4inlib <- installed.packages()
 #   if (length(l4inlib[l4inlib[,1]=="lme4"])<2){
 #   print("Need to install lme4 (Bates, 2008). Choose a mirror", quote=FALSE)
 #   install.packages("lme4")}
 #  library("lme4")}
 if (vardiff == FALSE && identical(item,NULL))
   mod1 <- lmer(sold ~ isold + (1|subno),family=binomial(link=lk),method="Laplace")
 if (vardiff && identical(item,NULL))
   mod1 <- lmer(sold ~ isold + (isold|subno),family=binomial(link=lk),method="Laplace")
 if (vardiff == FALSE && identical(item,NULL)==FALSE)
   mod1 <- lmer(sold ~ isold + (1|subno) + (1|item),family=binomial(link=lk),method="Laplace")
 if (vardiff  && identical(item,NULL)==FALSE)
   mod1 <- lmer(sold ~ isold + (isold|subno) + (1|item),family=binomial(link=lk),method="Laplace")
 if (lk=="probit") print("Overall estimated d' ",quote=FALSE)
 if (lk=="logit")print("Overall estimated lnOR ",quote=FALSE)
 print(format(fixef(mod1)[2],digits=4),quote=FALSE)
 print(summary(mod1)@coefs)
 if (is.null(covs)) return(mod1)
 ncovs <- max(1,dim(covs)[2])
 for (i in 1:ncovs){
   if (ncovs == 1) cov <- covs
   if (ncovs > 1) cov <- covs[,i]
   if (is.factor(cov) || is.ordered(cov)){
     if (vardiff == FALSE && identical(item,NULL)){
      mod2 <- lmer(sold~-1+cov+isold:cov+(1|subno),family=binomial(link=lk),method="Laplace")
      m1 <- lmer(sold~cov+isold+(1|subno),family=binomial(link=lk),method="Laplace")
      m2 <- lmer(sold~cov*isold+(1|subno),family=binomial(link=lk),method="Laplace")}
     if (vardiff == FALSE && identical(item,NULL)==FALSE){
      mod2 <- lmer(sold~-1+cov+isold:cov+(1|item)+(1|subno),family=binomial(link=lk),method="Laplace")
      m1 <- lmer(sold~cov+isold+(1|item)+(1|subno),family=binomial(link=lk),method="Laplace")
      m2 <- lmer(sold~cov*isold+(1|item)+(1|subno),family=binomial(link=lk),method="Laplace")}
     if (vardiff == TRUE && identical(item,NULL)){
      mod2 <- lmer(sold~-1+cov+isold:cov+(isold|subno),family=binomial(link=lk),method="Laplace")
      m1 <- lmer(sold~cov+isold+(isold|subno),family=binomial(link=lk),method="Laplace")
      m2 <- lmer(sold~cov*isold+(isold|subno),family=binomial(link=lk),method="Laplace")}
     if (vardiff == TRUE && identical(item,NULL)==FALSE){
      mod2 <- lmer(sold~-1+cov+isold:cov+(1|item)+(isold|subno),family=binomial(link=lk),method="Laplace")
      m1 <- lmer(sold~cov+isold+(1|item)+(isold|subno),family=binomial(link=lk),method="Laplace")
      m2 <- lmer(sold~cov*isold+(1|item)+(isold|subno),family=binomial(link=lk),method="Laplace")}
    if (lk == "probit") print("The estimated d' values for the groups are:",quote=FALSE)
    if (lk == "logit") print("The estimated lnOR values for the groups are:",quote=FALSE)
    vals <- length(unique(cov))
    print(summary(mod2)@coefs)
    print(anova(m1,m2,test = "Chisq"))
    }
   if (is.numeric(cov)){
    varc <- (cov-mean(cov))/sd(cov)
     if (vardiff == FALSE && identical(item,NULL)){
      mod2 <- lmer(sold~varc+isold:varc+(1|subno),family=binomial(link=lk),method="Laplace")
      m1 <- lmer(sold~varc+isold+(1|subno),family=binomial(link=lk),method="Laplace")
      m2 <- lmer(sold~varc*isold+(1|subno),family=binomial(link=lk),method="Laplace")}
     if (vardiff == FALSE && identical(item,NULL)==FALSE){
      mod2 <- lmer(sold~varc+isold:varc+(1|item)+(1|subno),family=binomial(link=lk),method="Laplace")
      m1 <- lmer(sold~varc+isold+(1|item)+(1|subno),family=binomial(link=lk),method="Laplace")
      m2 <- lmer(sold~varc*isold+(1|item)+(1|subno),family=binomial(link=lk),method="Laplace")}
     if (vardiff == TRUE && identical(item,NULL)){
      mod2 <- lmer(sold~varc+isold:varc+(isold|subno),family=binomial(link=lk),method="Laplace")
      m1 <- lmer(sold~varc+isold+(isold|subno),family=binomial(link=lk),method="Laplace")
      m2 <- lmer(sold~varc*isold+(isold|subno),family=binomial(link=lk),method="Laplace")}
     if (vardiff == TRUE && identical(item,NULL)==FALSE){
      mod2 <- lmer(sold~varc+isold:varc+(1|item)+(isold|subno),family=binomial(link=lk),method="Laplace")
      m1 <- lmer(sold~varc+isold+(1|item)+(isold|subno),family=binomial(link=lk),method="Laplace")
      m2 <- lmer(sold~varc*isold+(1|item)+(isold|subno),family=binomial(link=lk),method="Laplace")}
   if (lk == "probit")
      print("The estimated shift in d' for a standard deviation shift is:",quote=FALSE)
    if (lk == "logit")
      print("The estimated shift in lnOR for a standard deviation shift is:",quote=FALSE)
    print(fixef(mod2))
    print(anova(m1,m2,test = "Chisq"))}}
if (ncovs == 1 && vardiff && identical(item,NULL))
 fmint <- as.formula(paste("sold~isold*covs +(isold|subno)"))
if (ncovs == 1 && vardiff==FALSE && identical(item,NULL))
 fmint <- as.formula(paste("sold~isold*covs +(1|subno)"))
if (ncovs == 1 && vardiff && identical(item,NULL)==FALSE)
 fmint <- as.formula(paste("sold~isold*covs +(1|item)+(isold|subno)"))
if (ncovs == 1 && vardiff==FALSE && identical(item,NULL)==FALSE)
 fmint <- as.formula(paste("sold~isold*covs +(1|item)+(1|subno)"))
if (ncovs > 1){
 cv <- {}
 for (i in 1:ncovs) cv <- cbind(cv,covs[,i])
 dimnames(cv) <- dimnames(covs)
 cs <- paste("(")
 for (i in 1:ncovs){
   cs <- paste(cs,dimnames(covs)[[2]][i])
   if (i < ncovs && int) cs <- paste(cs,"*")
   if (i < ncovs && int==FALSE) cs <- paste(cs,"+")
   if (i == ncovs) cs <- paste(cs,")")}
 if (vardiff == FALSE && identical(item,NULL))
   fmint <- as.formula(paste("sold ~ (1|subno) + isold*",cs))
 if (vardiff && identical(item,NULL))
   fmint <- as.formula(paste("sold ~ (isold|subno) + isold*",cs))
 if (vardiff == FALSE && identical(item,NULL)==FALSE)
   fmint <- as.formula(paste("sold ~ (1|subno)+(1|item) + isold*",cs))
 if (vardiff && identical(item,NULL)==FALSE)
   fmint <- as.formula(paste("sold ~ (isold|subno) +(1|item)+ isold*",cs))
   }
final <- lmer(fmint,family=binomial(link=lk),method="Laplace")
return(final)}

