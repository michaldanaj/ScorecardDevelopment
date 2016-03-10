# TODO: Add comment
# 
# Author: michal.danaj
###############################################################################






pokaz_gini<-function(model, subset=rep(TRUE, nrow(woe_wynik)), target=newdata$target, newdata=woe_wynik){
	
	pred<-predict(model, type='response', newdata=newdata)
	#reg_nieparam(pred_train,sample0_int$target[sample0_int$train==TRUE])
	
	zostaw_train<-subset & !is.na(pred) & newdata$train==TRUE
	zostaw_test<-subset & !is.na(pred) & newdata$train!=TRUE
	
	#infromacje o poprawnoœci kalibracji na zbiorze testowym
	#informacje_kal(pred[zostaw_test], target[zostaw_test], pred[zostaw_test])
	
	#reg_nieparam(pred[zostaw_train], sample0_int$target[zostaw_train])
	
	#reg_nieparam(pred_test[zostaw],sample0_int$target[sample0_int$train!=TRUE][zostaw])
	rbind(
			AR(-pred[zostaw_train],target[zostaw_train])[[1]],
			AR(-pred[zostaw_test],target[zostaw_test], plot=TRUE)[[1]]
	)
}

saveState<-function(filename="modele_int.RDat", dir="D:\\Michal\\IA-0219 Model CP dla PB\\Testy R"){
	file=paste(dir,filename, sep='\\')
	save(modele, modele_kom, wyniki_modeli, file=file)
}

loadState<-function(filename="modele_int.RDat", dir="D:\\Michal\\IA-0219 Model CP dla PB\\Testy R"){
	file=paste(dir,filename, sep='\\')
	load(file)
}


## data=woe_wynik
## model=model40




## x=13
## wyniki_int[[x]][[1]][,1:6]
## genSQLCase(wyniki_int[[x]][[1]], names(wyniki_int[x]), "numeric")
## 
## 
## wynik<-wyniki_reczna_dyskr
## dane<-sample0_int
## 
## 
## genSQL(wyniki_reczna_dyskr, sample0_int, model17)






calcWeightForApriori<-function(sample_apriori, target_apriori){
	
	p<-sample_apriori #p - pocz¹tkowe prawdopodobieñstwo
	pp<-target_apriori
			
	#wyliczam wagê dla 0
	p*(1-pp)/((1-p)*pp)

}

calcOffsetForApriori<-function(sample_apriori, target_apriori){
	
	p<-sample_apriori #p - pocz¹tkowe prawdopodobieñstwo
	pp<-target_apriori
	
	#wyliczam wagê dla 0
	log(1/calcWeightForApriori(p,pp))
	
}


#' Funkcja multiplikuje obserwacje tak, aby osi¹gn¹æ zadane a-priori
#'
#' Funkcja multiplikuje obserwacje dla \code{target==0} tak, aby osi¹gn¹æ zadane a-priori 
#' \code{pp}.
#' Wynikiem jest proba ze zmultiplikowanymi obserwacjami ze zbioru \code{proba}
#' @param proba \code{data.frame} z wejœciow¹ prób¹
#' @param target wektor z wartoœciami \code{c(0,1)}
#' @param pp docelowe a-priori
transformToApriori<-function(proba, target, pp){
	
	p<-mean(target) #p - pocz¹tkowe prawdopodobieñstwo
	
	#wyliczam wagê dla 0
	waga0<-calcWeightForApriori(p,pp)
	print (waga0, digits=20)
	
	#oddzielam czêœæ ca³kowit¹ od u³amkowej, ¿eby zmultiplikowaæ obserwacje "0"
	#tyle razy ile jest w czêœci ca³kowitej + jeden raz z prawdopodobieñstwem równym
	#czêœæ u³amkowej
	
	waga0_calkowita<-trunc(waga0)
	waga0_ulamkowa<-waga0-waga0_calkowita
	
	
	#losujê, ile razy zmultiplikowaæ dany rekord
	losuj<-waga0_calkowita+rbinom(nrow(proba), 1, waga0_ulamkowa)
	
	#tam gdzie target==1, generujê tylko jedn¹ obserwacjê
	losuj[target==1]<-1
	
	#nr wierszy próby wejœciowej, do wygenerowania próby wyjœciowej
	idx<-rep(1:nrow(proba), times=losuj)		
	
	proba[idx,]
}



#' Norma L1
L1Norm<-function(x){
	mean(abs(x-mean(x)))
}




get_password <- function() {
	cat("Has³o, buraku!!!: ")
	system("stty -echo")
	a <- readline()
	system("stty echo")
	cat("\n")
	invisible(a)
}

#get_password()