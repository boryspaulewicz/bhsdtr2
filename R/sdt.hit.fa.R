##' Utility function
##' 
##' @param dprim d'
##' @param crit main decision criterion centered at the midpoint between the evidence distributions
##' @export
sdt.hit.fa = function(dprim, crit){
    if(length(dprim) > 1){
        fun = function(x, y)cbind(hit = x, fa = y)
    }else{
        fun = function(x, y)c(hit = x, fa = y)
    }
    fun(pnorm(dprim / 2 - crit),
        pnorm(-(dprim / 2 + crit)))
}
