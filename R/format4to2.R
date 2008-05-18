`format4to2` <-
function(x, cnames = c("subno","isold","saysold"),code=.5,...)
{
sold <- {}
iold <- {}
sub <- {}
if (dim(x)[2] == 4)
  x <- cbind(1:dim(x)[1],x)
for (i in 1:dim(x)[1]){
sold <- c(sold,rep(1,sum(x[i,2],x[i,3])),rep(0,sum(x[i,4],x[i,5])))
iold <- c(iold,rep(1,x[i,2]),rep(0,x[i,3]),rep(1,x[i,4]),rep(0,x[i,5]))
sub <- c(sub,rep(i,sum(c(as.numeric(x[i,2:5])))))  }
iold <- iold + code - 1
newvars <- cbind(sub,iold,sold)
dimnames(newvars)[[2]] <- cnames
return(as.data.frame(newvars))}

