# TODO: Add comment
# 
# Author: michal.danaj
###############################################################################






pokaz_gini<-function(model, subset=rep(TRUE, nrow(woe_wynik)), target=newdata$target, newdata=woe_wynik){
	
	pred<-predict(model, type='response', newdata=newdata)
	#reg_nieparam(pred_train,sample0_int$target[sample0_int$train==TRUE])
	
	zostaw_train<-subset & !is.na(pred) & newdata$train==TRUE
	zostaw_test<-subset & !is.na(pred) & newdata$train!=TRUE
	
	#infromacje o poprawno�ci kalibracji na zbiorze testowym
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
	
	p<-sample_apriori #p - pocz�tkowe prawdopodobie�stwo
	pp<-target_apriori
			
	#wyliczam wag� dla 0
	p*(1-pp)/((1-p)*pp)

}

calcOffsetForApriori<-function(sample_apriori, target_apriori){
	
	p<-sample_apriori #p - pocz�tkowe prawdopodobie�stwo
	pp<-target_apriori
	
	#wyliczam wag� dla 0
	log(1/calcWeightForApriori(p,pp))
	
}


#' Funkcja multiplikuje obserwacje tak, aby osi�gn�� zadane a-priori
#'
#' Funkcja multiplikuje obserwacje dla \code{target==0} tak, aby osi�gn�� zadane a-priori 
#' \code{pp}.
#' Wynikiem jest proba ze zmultiplikowanymi obserwacjami ze zbioru \code{proba}
#' @param proba \code{data.frame} z wej�ciow� pr�b�
#' @param target wektor z warto�ciami \code{c(0,1)}
#' @param pp docelowe a-priori
transformToApriori<-function(proba, target, pp){
	
	p<-mean(target) #p - pocz�tkowe prawdopodobie�stwo
	
	#wyliczam wag� dla 0
	waga0<-calcWeightForApriori(p,pp)
	print (waga0, digits=20)
	
	#oddzielam cz�� ca�kowit� od u�amkowej, �eby zmultiplikowa� obserwacje "0"
	#tyle razy ile jest w cz�ci ca�kowitej + jeden raz z prawdopodobie�stwem r�wnym
	#cz�� u�amkowej
	
	waga0_calkowita<-trunc(waga0)
	waga0_ulamkowa<-waga0-waga0_calkowita
	
	
	#losuj�, ile razy zmultiplikowa� dany rekord
	losuj<-waga0_calkowita+rbinom(nrow(proba), 1, waga0_ulamkowa)
	
	#tam gdzie target==1, generuj� tylko jedn� obserwacj�
	losuj[target==1]<-1
	
	#nr wierszy pr�by wej�ciowej, do wygenerowania pr�by wyj�ciowej
	idx<-rep(1:nrow(proba), times=losuj)		
	
	proba[idx,]
}



#' Norma L1
L1Norm<-function(x){
	mean(abs(x-mean(x)))
}




get_password <- function() {
	cat("Has�o, buraku!!!: ")
	system("stty -echo")
	a <- readline()
	system("stty echo")
	cat("\n")
	invisible(a)
}

#get_password()