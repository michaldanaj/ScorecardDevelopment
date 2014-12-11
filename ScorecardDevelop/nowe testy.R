# TODO: Add comment
# 
# Author: Piotr
###############################################################################


library(MDBinom)


n<-1000;
x<-rnorm(n);
eps<-rnorm(n)


y<-as.numeric(1/(1+exp(x+eps))>0.8)

dd<-data.frame(x,y)

model<-glm(y~x, family='binomial')


reg_nieparam(x,y)
points(x,predict(model, type='response'), col='blue')


mapowanie<-numeric_var_treatment(x,y)
mapowanie2<-numeric_var_treatment(x,y, interactive=TRUE)

przypisz2(x, mapowanie, fitted=mapowanie$label)




reg_nieparam_bagging<-function (score, default, petli=20, buckets = 100, wytnij = 0, span = 0.7, 
		degree = 2, plot = TRUE, target = "br", new = TRUE, col_points = "black", 
		col_line = "darkblue", index = FALSE, ...) 
{
	dane <- data.frame(score, default)
	if (wytnij > 0) {
		do_usuniecia <- usun_konce(dane$score, prob = wytnij)
		if (length(do_usuniecia) > 0) 
			dane <- dane[-do_usuniecia, ]
	}
	bucket <- buckety_br(dane$score, dane$default, buckets, method = "eq_count")
	
	wynik.matrix<-matrix(nrow(bucket)*petli, ncol = petli, nrow = nrow(bucket))
	
	model_calosc <- locfit(default ~ lp(score, nn = span), family = "binomial", 
			link = "logit", data = dane)
	
	for (i in 1:petli){
		proba.idx<-sample(1:nrow(dane),size=nrow(dane), replace=TRUE)
		if (length(unique(default)) == 2) 
			l <- locfit(default ~ lp(score, nn = span), family = "binomial", 
					link = "logit", data = dane[proba.idx,])
		else l <- locfit(default ~ lp(score, nn = span), data = dane[proba.idx,])
		wynik.matrix[,i] <- predict(l, newdata = bucket$srodek)
	}
	b2<-rowMeans(wynik.matrix)
	if (target == "br") 
		bucket2 <- cbind(bucket, fitted = b2)
	else bucket2 <- cbind(bucket, fitted = log(b2/(1 - b2)))
	skala <- sqrt(bucket$n_obs/(length(score)/buckets))
	x <- bucket2$srodek
	if (index) 
		x <- bucket$nr
	if (plot) {
		if (new == TRUE) 
			plot(x, with(bucket2, get(target)), col = col_points, 
					cex = skala, ...)
		else points(x, with(bucket2, get(target)), cex = skala, 
					col = col_points, ...)
		matlines(x,wynik.matrix, new=TRUE,  col='black', lty='dotted')
		lines(x, bucket2$fitted, col = col_line, lwd=3,...)
		lines(x,predict(model_calosc, newdata=x), col='red')
	}
	bucket2
}



###############################################################################
########                                                               ########
###############################################################################


target<-c(1,1,0,0)
predicted<-c(0.9,0.9,0.1,0.1)
		
model_dev(target, predicted)
wytnij<-dane[,c('target_pwa','gprs_liczba_m1')]
model<-glm(target_pwa~gprs_liczba_m1, data=wytnij, family='binomial')

y2<-as.numeric(1/(1+exp(x*x+eps))>0.8)
		
model<-locfit(y2~x, family='binomial', deg=2)
summary(model)
wezly<-lfknots(model)
wezly<-lfknots(model)[,1]
wezly.y<-predict(model, wezly)
plot(model)
points(wezly, wezly.y, col='red')
rug(wezly, col='red')

reg_nieparam(x,y2, buckets=20)
points(model, col='red')

lfknots(model,what=c('x','coef','f1','nlx','nlx1','se','infl','infla','lik','h','deg'))

plot(lfeval(model),txt=TRUE)





options(show.error.messages = FALSE)
z=log("a")
print(.Last.value)
options(show.error.messages = TRUE)

## alternatively,
wyn<-try(log('c'), TRUE)

## run a simulation, keep only the results that worked.
set.seed(123)
x <- stats::rnorm(50)
doit <- function(x)
{
	x <- sample(x, replace=TRUE)
	if(length(unique(x)) > 30) mean(x)
	else stop("too few unique points")
}
## alternative 1
res <- lapply(1:100, function(i) try(doit(x), TRUE))
## alternative 2
## Not run: res <- vector("list", 100)
for(i in 1:100) res[[i]] <- try(doit(x), TRUE)
## End(Not run)

unlist(res[sapply(res, function(x) !inherits(x, "try-error"))])




--------------------------------------------------------------------------------
		





		
		
		###############################################################################
		########                     równoleg³e                             ########
		###############################################################################
		funkcja<-function(x, target){
	
	
	try(univariate_anal_stats1b(x, target,
					locfit=TRUE,
					NA_subst = numeric_var_treatment.params$NA_substit, 
					special_val = numeric_var_treatment.params$NA_substit), TRUE)
}


foreach(zmienna = zmienne_num_names)%dopar%{
	x=dane2[,zmienna]
	y=dane2$target_pwa
	funkcja(x, y)
}

library(foreach)
foreach(i=1:3) %do% sqrt(i)


library(snow)
library(doSNOW)
library(randomForest)
library(snowfall)

#cl<-makeCluster(c("localhost","localhost","localhost","localhost","localhost","localhost","localhost"), type = "SOCK")
sfInit(parallel=TRUE,cpus=4)
cl = sfGetCluster()
getDoParWorkers()

registerDoSNOW(cl)
xx <- matrix(runif(500), 100)
y <- gl(2, 50)

las<-randomForest(xx,y,ntree=250)
sfExport("y")
sfExport("x")
sfExport("las")
sfLibrary(randomForest)
parLapply(cl,seq(1,60),fun=las)

system.time(foreach(ntree=rep(250, 60), .combine=combine, .packages='randomForest') %dopar%
				randomForest(x, y, ntree=ntree))
stopCluster(cl)
sfStop()

install.packages(pkgs='d:/michal/snowfall_1.84.zip')






abc<-data.frame('xysjfasdklfjasklfjjhalfhjladfhasdjkfajfkafjkd'=3)

