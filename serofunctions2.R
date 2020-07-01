log_fr <- function(r,pos,neg,u,v){
	t <- u+(1-u-v)*r
    return(pos*log2(t) + neg*log2(1-t))
}



sample_post_r_log <- function(pos,neg,se,sp,size){

	size <- 2*size # make the first half burn-in

	u <- 1-sp
	v <- 1-se

	rm <- (pos/(pos+neg)-u)/(1-u-v)

	if(rm>0 && rm<1){
		log_M <- max(log_fr(0,pos,neg,u,v), log_fr(1,pos,neg,u,v), log_fr(rm,pos,neg,u,v), na.rm=TRUE)

	} else {
		log_M <- max(log_fr(0,pos,neg,u,v), log_fr(1,pos,neg,u,v), na.rm=TRUE)
	}

	# uniform accept/reject algorithm

	samples <- c()
	n_accepted <- 0
	max_failures <- 1e-6
	failures <- 0
	i <- 0

	while(n_accepted < size){
		if(i%%1000 == 0){
			prog <- 100*n_accepted/size
		}

		r_proposal <- runif(1)
		logrv <- log2(runif(1))

		if(logrv < (log_fr(r_proposal,pos,neg,u,v)-log_M)){
			samples <- c(samples, r_proposal)
			n_accepted <- n_accepted + 1
			failures <- 0
		} else {
			failures <- failures + 1
		}
		if(failures == max_failures){
			return(-1)
		}
		i <- 1 + 1
	}

	return(tail(samples,size/2)) # return the last half of the accepted values

}