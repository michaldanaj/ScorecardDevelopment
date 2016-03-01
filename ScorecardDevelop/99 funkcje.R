# TODO: Add comment
# 
# Author: michal.danaj
###############################################################################


dopasowanie_do_zmiennej<-function( x, bucket, estim, subset=NULL,...){
	if (!is.null(subset)){
		x<-x[subset]
		estim<-estim[subset]
	}
	if (any(is.na(estim))){
		ile_na<-sum(is.na(estim))
		warning(sprintf("W 'estim' bylo %s brakow danych. Zostaly usuniete.", ile_na))
		x<-x[!is.na(estim)]
		estim<-estim[!is.na(estim)]
	}
	bucket_new<-bucket
	bucket$fitted<-rownames(bucket)
	pred<-przypisz2(x, bucket)
	estim_bucket<-tapply(estim, pred, mean)
	bucket_new$model<-estim_bucket[rownames(bucket_new)]
	plot(bucket_new$nr[-nrow(bucket_new)], bucket_new$br[-nrow(bucket_new)],...)
	points(bucket_new$nr, bucket_new$model, col="green")
	bucket_new
}




korelacje_zmiennych<-function(model, data){
	zmienne<-names(coef(model))[-1]
	edit(cor(data[,zmienne]))
}


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

make_model_formula<-function(target, vars){
	
	#je�li w nazwach zmiennych jest nazwa targetu, to j� usuwam
	czy_jest_target<-target==vars
	vars<-vars[!czy_jest_target]
	
	#suma zmiennych
	suma_zmiennych<-paste(vars, collapse = ' + ')
	
	#wynik
	formula(paste(target, suma_zmiennych, sep='~'))
}


## data=woe_wynik
## model=model40

#licz� korelacje zmiennych z modelu ze zmiennymi do dodania
step_bez_kor<-function(data, model){
	
	#zmienne w modelu
	zmienne_model<-names(coef(model)[-1])
	
	#zmienne z danych z budowy modelu
	#zmienne_budowa<-nazwy_zmiennych
	zmienne_budowa<-names(data)
	
	#korelacja mi�czy nimi
	korel<-as.data.frame(cor(data[,zmienne_budowa]))
	korel_zm_model<-korel[zmienne_model,]
	
	#gdzie akceptowalna korelacja
	korelacje_max<-apply(abs(korel_zm_model), 2, max)
	czy_przekracza<-as.data.frame(abs(korel_zm_model)>0.75)
	czy_przekracza<-sapply(czy_przekracza, any)
	zmienne<-names(czy_przekracza[czy_przekracza==FALSE & !is.na(czy_przekracza)])
	
	#usuwam target
	zmienne<-zmienne[zmienne!='target']
	
	#robi� stepa
	#form<-make_model_formula('target',zmienne_budowa)
	form<-make_model_formula('target', c(".",zmienne))
	dodana_zmienne<-add1(model, scope= form, test='Chisq')	
	kolejnosc_aic<-order(dodana_zmienne$AIC)
			
	cbind(dodana_zmienne[kolejnosc_aic,], cor_max=korelacje_max[rownames(dodana_zmienne[kolejnosc_aic,])])
}



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

genSQL<-function(wyniki, dane, model){
	
	#nazwy zmiennych wynikowych z case'�w
	varNames<-paste("var", 1:length(wyniki), sep="")
	sql<-c(
		genSQLCases(wyniki, dane, varNames),
		",",
		genSQLScore(model, varNames)
	)
	write.table(sql, file='clipboard', quote=FALSE, row.names=FALSE, col.names=FALSE)
	sql
}

## 
## x=3
## genSQLCase(bucket=wyniki_reczna_dyskr[[x]][[1]], 
##         nazwaZmiennej=names(wyniki_reczna_dyskr[x]), 
##         typZmiennej="numeric", 
##         nazwaZmiennejOut='xyz')

genSQLCase<-function(bucket, nazwaZmiennej, typZmiennej, nazwaZmiennejOut){

	#je�li jest wiersz z podsumowaniem, to go usuwam
	bucket<-bucket[bucket$label != 'TOTAL',] 
	 
	#robi� zaokr�glenie
	bucket$woe<-round(bucket$woe, 7)
	
	kodzik_dyskr<-character()
	kodzik_ciagle<-character()
	kodzik_null<-character()
	
	#Najpierw obs�uguj� warto�ci dyskretne
	gdzie_dyskretne<-is.na(bucket$srodek)
	if (any(gdzie_dyskretne)){

		ciapek=''
		if (typZmiennej=='character')
			ciapek="'"

		dyskretne <- bucket[gdzie_dyskretne,]
		kodzik<-paste("WHEN ", nazwaZmiennej, ' = ', ciapek, dyskretne$discret, ciapek, ' THEN ', dyskretne$woe, sep='')
				
		#podmieniam dla warto�ci null
		#TODO Obs�u�y� to lepiej! Najpierw wygenerowa� woe dla null, w zale�no�ci czy jawne czy nie
		#a dopiero p�niej wygenerowa� kod w zale�no�ci od typu zmiennej
		gdzie_null<-dyskretne$discret==numeric_var_treatment.params$NA_substit | dyskretne$discret==''
		kodzik[gdzie_null]<-paste("WHEN", nazwaZmiennej, 'IS NULL', 'THEN', dyskretne$woe[gdzie_null])
		kodzik_dyskr<-kodzik

	}
	
	#teraz zmienne ci�g�e	
	gdzie_ciagle<- !is.na(bucket$srodek)
	
	if (any(gdzie_ciagle)){
		
		ciagle <- bucket[gdzie_ciagle,]
		kodzik<-paste("WHEN", nazwaZmiennej, '>=', ciagle$od,'AND', nazwaZmiennej, '<', ciagle$do,'THEN', ciagle$woe)
		
		#podmieniam pierwszy wiersz, aby by� od minus niesko�czono�ci
		kodzik[1]<-paste("WHEN", nazwaZmiennej, '<', ciagle$do[1],'THEN', ciagle$woe[1])
		
		#podmieniam ostatni wiersz, aby by� od minus niesko�czono�ci		
		kodzik[nrow(ciagle)]<-paste("WHEN", nazwaZmiennej, '>=', ciagle$od[nrow(ciagle)],'THEN', ciagle$woe[nrow(ciagle)])
		
		kodzik_ciagle<-kodzik
		
	}
	
	#jeszcze else na koniec
	kodzik_else<-"     ELSE NULL"
				
	#generuj� ca�y kod, dodaj�c warunek else
	kodzik_all<-c('CASE',
			paste("    ",c(kodzik_dyskr, kodzik_ciagle, kodzik_null)),
			kodzik_else,
		paste('END as ', nazwaZmiennejOut,'\n', sep=''))	

	kodzik_all_cat<-paste(kodzik_all, collapse='\n')


	kodzik_all_cat
	
}

#genSQLCases(wyniki_reczna_dyskr, sample0_int, varNames)

#TODO doda� do struktur w dyskretyzacj�, oraz do metadanych, typ danych
#Na razie robi� to r�cznie. Dlatego potrzebuj� zmeinnej \code{dane}
#aby okre�li� typy danych
genSQLCases<-function(wynik, dane, varNames){
	
	cases<-list()
	
	#wyznaczam typ danych, p�niej trzeba by to zast�pi� metadanymi
	for (i in 1:length(wynik)){
		typ_danych<-typeof(dane[,names(wynik[i])])
		
		cases[[i]]<-genSQLCase(wynik[[i]][[1]], names(wynik[i]), typ_danych, nazwaZmiennejOut = varNames[i])
		names(cases)[i]<-names(wynik[i])
		
	}
	
	
	## cases<-lapply(wynik, function(x){
	##             genSQLCase(x[[1]], names(x), attributes(x[[1]])$type)
	##         }
	## )
	
	cases2<-character()
	for (i in 1:length(wynik)){
		cases2<-c(cases2, cases[[i]])
		if (i<length(wynik))
			cases2<-c(cases2,",")		
	}
	
	cat("\n\n#########   	Kod z przekszta�ceniami do wklejenia:   		############\n\n")
	cat(cases2)
	cat("\n\n#########   Koniec kodu z przekszta�ceniami do wklejenia   ############\n\n")
	
	cases2
}

#TODO doda� do struktur w dyskretyzacj�, oraz do metadanych, typ danych
#Na razie robi� to r�cznie. Dlatego potrzebuj� zmeinnej \code{dane}
#aby okre�li� typy danych
# mapping - mapowanie nazwy �r�d�owej (przed zastosowniem dyskretyzacji z \code{wynik}
# 			na nazw�, kt�ra jest p�niej w modelu. Je�li jedna i druga nazwa jest taka sama,
#			to znaczy �e nie by�o robionej dyskretyzacji i case jest dla tej zmeinnej
#			pomini�ty
genSQLCases2<-function(wynik, dane, mapping){
	
	cases<-list()
	
	#ograniczam list� wynik tylko do tych zmiennych, kt�re maj�
	#zmienion� nazw� w mapping
	rozneNazwy<-mapping$sourceVarName != mapping$modelVarName
	wynik<-wynik[mapping$sourceVarName[rozneNazwy]]
	
	#wyznaczam typ danych, p�niej trzeba by to zast�pi� metadanymi
	for (i in 1:length(wynik)){
		typ_danych<-typeof(dane[,names(wynik[i])])
		
		cases[[i]]<-genSQLCase(wynik[[i]][[1]], names(wynik[i]), typ_danych, nazwaZmiennejOut = varNames[i])
		names(cases)[i]<-names(wynik[i])
		
	}
	
	
	## cases<-lapply(wynik, function(x){
	##             genSQLCase(x[[1]], names(x), attributes(x[[1]])$type)
	##         }
	## )
	
	cases2<-character()
	for (i in 1:length(wynik)){
		cases2<-c(cases2, cases[[i]])
		if (i<length(wynik))
			cases2<-c(cases2,",")		
	}
	
	cat("\n\n#########   	Kod z przekszta�ceniami do wklejenia:   		############\n\n")
	cat(cases2)
	cat("\n\n#########   Koniec kodu z przekszta�ceniami do wklejenia   ############\n\n")
	
	cases2
}


genSQLModel<-function(model, varNames=NULL){
	
	wspolczynniki<-round(coef(model),7)
	
	if (is.null(varNames))
		varNames=names(coef(model))[-1]
		
	#dodaj� nazw� NA na pocz�tek, bo we wsp�czynnikach mam dodatkowo itercept
	varNames<-c(NA,varNames)
	
	 
	score_lin<-paste(wspolczynniki,' * ', varNames, '+')
	
	#podmieniam intercept
	score_lin[1]<-paste(wspolczynniki[1],'+')
	
	#usuwam + z ostatniego wiersza
	i<-length(score_lin)
	score_lin[i]<-paste(wspolczynniki[i],' * ', varNames[i])
	
	score_lin<-c(score_lin, 'as score_lin')
	
	score_lin<-c(score_lin,"score_pd = 1/(1+exp(-score_lin));")
	
	score_lin_cat<-paste(score_lin, collapse="\n")
	
	cat("\n\n#########   	Kod scoringowy do wklejenia:   		############\n\n")
	cat(score_lin_cat)
	cat("\n\n#########   Koniec kodu scoringowy do wklejenia   ############\n\n")
	score_lin_cat
}	







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






przypisz_z_listy<-function(bucket_list, data, vars=names(bucket_list), colname='fitted', varname_sufix=colname){
	
	data_out<-NULL
	
	for (zmienna in names(bucket_list)){
		
		if (!(zmienna %in% vars)) 
			next;
		
		####   wyliczam   woe    ######
		
		#Wyci�gam element listy
		bucket<-bucket_list[[zmienna]]
		#bucket<-bucket[rownames(bucket)!='TOTAL',]
		
		#je�li bad�w lub good�w jest 0, to przyjmuj� �e jest 0.5	
		
		fitted = bucket[,colname]
		
		####   przypisuj� woe    ######
		
		fitted_x<-przypisz2(data[,zmienna],
				bucket_list[[zmienna]], 
				fitted=fitted,
				NA_subst = numeric_var_treatment.params$NA_substit,
				interpol=FALSE)
		
		if (is.null(data_out))	{
			data_out <- data.frame(fitted_x)
			names(data_out)<-zmienna
		}
		else
			data_out[,zmienna] <- c(fitted_x)
		
	}
	
	names(data_out)<-paste(names(data_out), varname_sufix, sep="_")
	data_out
}	



get_password <- function() {
	cat("Has�o, buraku!!!: ")
	system("stty -echo")
	a <- readline()
	system("stty echo")
	cat("\n")
	return(a)
}