`format2to4` <-
function(subno,isold,sold,
    cnames = c("subno","hits","fas","misses","crs"),...)
{
if (dim(table(isold,sold))[1]!= 2 || dim(table(isold,sold))[2]!= 2){
print("Either there is a problem with the format entered or " )
print("a row or column marginal is zero for the whole set.")
print("Since the function will not know the meaning of value")
print("it will return NULL.")
return(NULL)
}
stopifnot((length(subno)== length(isold))&(length(isold)==length(sold)))
data <- matrix(nrow=length(unique(subno)),ncol=5)
hit <- (isold==max(isold))&(sold==max(sold))
fa <- (isold==min(isold))&(sold==max(sold))
cr <- (isold==min(isold))&(sold==min(sold))
miss <- (isold==max(isold))&(sold==min(sold))
hits <- tapply(hit,subno,sum)
fas <- tapply(fa,subno,sum)
crs <- tapply(cr,subno,sum)
misses <- tapply(miss,subno,sum)
data <- cbind(unique(subno),hits,fas,misses,crs)
colnames(data) <- cnames
rownames(data) <- {}
return(data)}

