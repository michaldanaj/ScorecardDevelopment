#library(sqldf);
library(MDBinom);
library(reshape);
#library(kalibracja);
library(R2HTML)
library(lattice)
library(Hmisc)





#########################################

numeric_var_treatment.params<-list(
		#przy jakiej wartoœci unikalnych wartoœci zmiennej ma j¹ traktowaæ jako dyskretn¹
		discret_threshold=15,
		
		#kody ró¿nych wartoœci specjalnych. Bêd¹ one traktowane osobno, jako wartoœci dyskretne
		spcial_val=-.Machine$integer.max,
		
		#wartoœæ do zast¹pienia missing value
		NA_substit=-.Machine$integer.max,
		
		#Graniczny udzia³ wartoœci, powy¿ej której traktujemy j¹ jako wartoœæ specjaln¹ 
		#(w sposób dyskretny, wydzielon¹ z pozosta³ych).
		separate_value_thr=0.1,		
		
		#maksymalna g³êbokoœæ budowy drzewa moim algorytmem
		max_gleb=3,
		
		#minimalna liczba obserwacji - do sprawdzenia - w liœciu/w wêŸle do podzia³u
		min_bucket=200,
		
		#wartoœæ graniczna nulli. Poni¿ej robimy imputacjê, powy¿ej traktujemy je jako osobn¹ grupê
		nulle_do_imp_thr=0.0

)



###########################################
#  zmiana funkcji do dat. Dodanie okresu pó³rocznego #
###########################################

#' Zwraca ostatni dzieñ podanego okresu.
#'
#' @param dat Obiekt klasy \link{Date}, lub daj¹cy siê skonwertowaæ funkcj¹ \link{as.Date}.
#' @param unit Jednostka czasu, która ma byæ uwzglêdniona przy wyznaczaniu ostatniego dnia.
#'  Mo¿liwe wartoœci: \code{c("month","quater", "year")}, które odpowiadaj¹ ostatniemu dniu miesi¹ca, kwarta³u, roku,
#'  do którego nale¿y data \code{dat}.
#' @param ... Opcje do funkcji \link{as.Date}.
#' @author Micha³ Danaj
lastDay<-function(dat, unit=c("month","quarter", "half_year", "year"),...) {
	
	# Wstêpne przygotowania
	unit<-match.arg(unit);
	
	dat<-as.Date(dat,...);
	
	# Pobieram rok i miesi¹c
	rok<-as.numeric(format(dat, '%Y'));
	mies<-as.numeric(format(dat, '%m'));
	
	# Jeœli ostatni dzieñ kwarta³u, to wybieram ostatni miesi¹c kwarta³u
	if (unit=="quarter")
		mies<-((mies-1)%/%3+1)*3
	
	if (unit=="half_year")
		mies<-((mies-1)%/%6+1)*6
	
	if (unit=="year")
		mies<-12;
	
	# W³aœciwy algorytm
	ostatni<-mies==12;
	mies[!ostatni]<-mies[!ostatni]+1;
	mies[ostatni]<-1;
	rok[ostatni]<-  rok[ostatni]+1;
	
	# zamiana na datê
	new_date<-as.Date(ISOdate(rok, mies, 1));
	return(new_date-1);
}


###########################################
#  zmiana funkcji do rysowania kalibracji #
###########################################
plotCalibr<-
		function (score, response, estim, target = c("br", "logit"), ylab=y_name,
				xlab="score",    ...)
{
	target <- match.arg(target)
	if (!is.list(estim))
		estim <- list(estim)
	y_name <- "PD"
	if (target == "logit")
		y_name <- "logit"
	buck <- reg_nieparam(score, response, target = target, ylab = ylab,
			xlab = xlab, ...)
	grupy_skala <- data.frame(od = buck$od, do = buck$do, fitted = buck$nr)
	grupy <- przypisz(score, grupy_skala)
	kolejnosc <- order(score)
	for (i in 1:length(estim)) {
		estim_grp <- tapply(estim[[i]], grupy, mean)
		buck <- cbind(buck, estim_grp)
		names(buck) <- c(names(buck)[1:(ncol(buck) - 1)], paste(names(estim[i]),
						"br", sep = "_"))
		if (target == "logit") {
			estim[[i]] <- logit(estim[[i]])
			buck <- cbind(buck, logit(estim_grp))
			names(buck) <- c(names(buck)[1:(ncol(buck) - 1)],
					paste(names(estim[i]), "logit", sep = "_"))
		}
		lines(score[kolejnosc], estim[[i]][kolejnosc], col = i +
						1, lty = i + 1, pch = i + 19)
	}
	legend(x = "topright", col = 1:(length(estim) + 1), lty = 1:(length(estim) +
						1), legend = c("Dependence from data", names(estim)),
			bty = "n", inset = 0.05)
	return(buck)
}




NA_substit=-.Machine$integer.max;
special_val=c(-99999999, -10000001, -10000000, -9999999,  -9999998,  -9999997, -9999996);

na.subst<-function(x, val){
	x[is.na(x)]<-val;
	x
}

#' Dyskretyzuje zmienn¹ i wylicza na niej statystyki
#'
#' W przypadku, gdy liczba unikalnych wartoœci zmiennej jest <= \code{discret_treshold}
#' lub zmienna nie jest zmienn¹ numeryczn¹,
#' uznaje ¿e zmienna jest dyskretna i jedynie wylicza dla niej statystyki. W przeciwnym
#' wypadku dyskretyzuje zmienn¹ i wylicza statystyki.
#' @item discret_treshold jeœli liczba unikalnych wartoœci zmiennej jest nie wiêksza
#'        ta wartoœæ, zmienna uznana jest za dyskretn¹ i nie jest poddawana dyskretyzacji.
#' @item interactive TRUE, jeœli zmienna ma byæ dyskretyzowana interaktywnie. W
#'                   przeciwnym razie, co jest wartoœci¹ domyœln¹, dyskretyzacja
#'                   jest automatyczna.
#' @item breaks zamiast automatycznego dzielenia, mo¿na podaæ wartoœci przedzia³ów (from,to].
#' @item forceContinous wymusza potraktowanie zmiennej jako ci¹g³¹, mimo ¿e liczba
#'                      unikalnych wartoœci jest mniejsza ni¿ \code{discret_treshold}.
#' @seealso \link{\code{buckety_stat}}.

univariate_anal_stats1<-function(x,y, discret_treshold=15,
		special_val=NULL, max_gleb=3, min_bucket=200, interactive=FALSE,
		breaks=NULL, mapping=NULL, forceContinous=FALSE,...){
	
	if (length(x)!=length(y))
		stop("paramet ry 'x' i 'y' maj¹ ró¿ne d³ugoœci!");
	
	#Mimo, ¿e przygotowywya³em funkcjê do obs³ugi null-i, to rezygnujê z tego
	#ze wzglêdów bezpieczeñstwa.
	if (any(is.na(y)))
		stop ("W 'y' nie mo¿e byæ NA!");
	
	## jeœli jest to ju¿ zmienna dyskretna;
	
	## jeœli jest to zmienna dyskretna
	if ((length(unique(x))<=discret_treshold || !is.numeric(x))&&
			is.null(breaks) && !forceContinous){
		discret<-buckety_stat(x, y, total=TRUE);
		#uzupe³niam statystyki
		discret$od<-NA;
		discret$do<-NA;
		discret$srodek<-NA;
		
		nam<-rownames(discret)
		if (is.numeric(x)){
			nam[length(nam)]<-NA
			discret$discret<-as.numeric(nam)
		}
		else{
			nam[length(nam)]<-"<TOTAL>";
			discret$discret<-nam;
		}
		
		discret<-discret[,c('nr','label','discret', 'od','srodek','do','n_good','pct_good','n_bad','pct_bad','n_obs','pct_obs',
						'br','woe','logit')]
	}
	## jeœli jest to zmienna ci¹g³a
	else{
		discret<-discretization(x,y, special_val=special_val,
				max_gleb=max_gleb,min_bucket=min_bucket,breaks=breaks,
				interactive=interactive,...);
	}
	
	discret$label<-rownames(discret);
	return(discret);
}

#' Robi rozk³ad zmiennej po czasie (lub innym podziale)
#'
#' Wylicza licznoœci dla ka¿dego poziomu \code{x_discr} oraz œredni¹ wartoœæ \code{y}
#' w podziale na zadane grupy czasowe \code{czas} (lub podzia³ innego typu).
#' @item x_discr zmienna objaœniaj¹ca.
#' @item y zmienna odpowiedzi.
#' @item czas.
univariate_anal_stats2<-function(x_discr, y, czas, estim){
	
	#licznoœci po czasie i zmiennej
	total<-table(x_discr);
	total_czas<-table(czas);
	all_tbl<-cbind(table(x_discr, czas), TOTAL=total);
	all_tbl<-rbind(all_tbl, c(total_czas, length(x_discr)));
	rownames(all_tbl)[nrow(all_tbl)]<-'TOTAL';
	
	
	#rozk³ad ka¿dego bucketu po czasie (agreguj¹c po czasie zsumuje siê do 1)
	pct_all_tbl<-sweep(all_tbl, 2, STATS=c(total_czas,length(x_discr)), FUN="/");
	
	#œrednie LGD po czasie i zmiennej
	total<-tapply(y, list(x_discr), mean)
	avg_y<-tapply(y, list(x_discr,czas), mean);
	avg_y<-cbind(avg_y, TOTAL=total);
	total_czas<-c(tapply(y, list(czas), mean), mean(y));
	avg_y<-rbind(avg_y, total_czas);
	rownames(avg_y)[nrow(avg_y)]<-'TOTAL';
	return(list(obs_all_tbl = all_tbl, pct_all_tbl=pct_all_tbl, avg_t_tbl=avg_y,
					estim=tapply(estim, czas, mean)));
}

#' Wylicza AR po czasie i w podziale na zadane próby
#'
#' @item score zmienna, po której procedura bêdzie sortowaæ.
#' @item y zmienna odpowiedzi.
#' @item czas podzia³ na okresy czasowe.
#' @item proby \code{data.frame}, w którym ka¿da kolumna jest wektorem logicznym,
#'       zawieraj¹cym informacje, czy obserwacja nale¿y do danej próby.
univariate_anal_stats3<-function (score, y, czas, proby){
	razem<-data.frame(score,y,czas);
	razem$czas<-as.factor(razem$czas);
	
	if (!is.data.frame(proby)){
		proby<-as.data.frame(proby);
	}
	
	
	for (i in 1:ncol(proby)){
		proba<-razem[proby[,i], ];
		AR_calosc<-AR(proba$score, proba$y)[[1]]['AR'];
		aery<-sapply(split(proba, proba$czas), FUN=function(pr){
					
					if (nrow(pr)==0)
						return(NA);
					wyn<-AR(pr$score, pr$y)[[1]]['AR'];
					names(wyn)<-NULL;
					return(wyn)
				}
		);
		
		if (i==1){
			wyniki_AR<-data.frame(aery);
			wyniki_AR<-cbind(AR_calosc,t(wyniki_AR));
		}
		else
			wyniki_AR<-rbind(wyniki_AR, c(AR_calosc,aery));
	}
	
	rownames(wyniki_AR)<-colnames(proby);
	colnames(wyniki_AR)<- colnames(wyniki_AR);
	return(wyniki_AR);
}

#' Dyskretyzuje zmienn¹ i wylicza dla niej statystyki
#'
#' Wylicza statystyki i zwraca je w postaci listy.
#' @item x zmienna, po której procedura bêdzie sortowaæ.
#' @item y zmienna odpowiedzi.
#' @item czas podzia³ na okresy czasowe.
#' @item proby \code{data.frame}, w którym ka¿da kolumna jest wektorem logicznym,
#'       zawieraj¹cym informacje, czy obserwacja nale¿y do danej próby.
#' @item interactive TRUE, jeœli zmienna ma byæ dyskretyzowana interaktywnie. W
#'                   przeciwnym razie, co jest wartoœci¹ domyœln¹, dyskretyzacja
#'                   jest automatyczna.
#' @item min_bucket minimalna liczba obserwacji w buckecie, przy dzieleniu drzewem.
#' @item breaks zamiast automatycznego dzielenia, mo¿na podaæ wartoœci przedzia³ów (from,to].
#' @item forceContinous wymusza potraktowanie zmiennej jako ci¹g³¹, mimo ¿e liczba
#'                      unikalnych wartoœci jest mniejsza ni¿ \code{discret_treshold}.
#' @item NA_subst wartoœæ jaka ma byæ przypisana w miejsce braków danych. Dalsze analizy
#'		bêd¹ przeporwadzone w standardowy sposób. Jeœli jednak wartoœæ \code{NA_subst} zostanie
#' 		dodana do \code{special_val}, zostanie ona potraktowana jako wartoœæ dyskretna.  
#' @author Micha³ Danaj
univariate_anal_stats<-function(x,y,czas,proby=rep(TRUE, length(y)),
		interactive=FALSE, min_bucket=floor(0.05*length(x)), breaks=NULL,
		mapping=NULL, forceContinous=FALSE, 
		special_val=numeric_var_treatment.params$spcial_val, 
		NA_subst=numeric_var_treatment.params$NA_substit,
		span=0.9){
	
	# zamieniam braki danych na liczbê.
	if (is.numeric(x) & !is.null(NA_subst))
		x[is.na(x)]<-NA_subst;


	# dyskretyzujê zmienn¹ i wyliczam pierwsze statystyki
	stat1<-univariate_anal_stats1b(x,y, special_val=special_val,
			max_gleb=3,plot=FALSE, min_bucket=min_bucket,
			interactive=interactive, breaks=breaks, mapping=mapping,
			forceContinous=forceContinous, span=span);
	
	stat2<-NULL;
	stat3<-NULL;
	#Dalsze statystyki robiê pod warunkiem, ¿e jest wiêcej ni¿ jedna wartoœæ dyskretna
	#(W stat1 jest te¿ <TOTAL>, dlatego 2)
	if (nrow(stat1)>2){
		
		# przypisujê nazwê bucketu
		stat1$fitted<-stat1$label;
		x_discr<-przypisz2(x,stat1);
		
		# przypisujê BR
		stat1$fitted<-stat1$br;
		BR_discr<-przypisz2(x,stat1);
		
		# wyliczam drugie statystyki
		stat2<-univariate_anal_stats2(x_discr, y, czas, BR_discr);
		
		# wyliczam trzecie statystyki (GINI) po zadanych próbach i czasie
		
		stat3<-univariate_anal_stats3(score=-BR_discr, y, czas, proby);
	}
	return(list(dyskretyzacja=stat1, rozklady=stat2, dyskryminacja=stat3));
}


#'   Dyskretyzuje zmienn¹ ciag³¹ drzewkiem w oparciu o zmienn¹ odpowiedzi
#'
#' W parametrze \code{x} ani \code{y} nie mo¿e byæ NULLi. Mo¿na je jakoœ zakodowaæ.
#' @item x zmienna ci¹g³a.
#' @item y zmienna celu (PD, LGD).
#' @item special_val wartoœci specjalne, które powinny zostaæ uwzglêdnione jako
#'                   osobne klasy. Wartoœæ tak¹ jest równie¿ \code{NA}, automatycznie
#'                   uwzglêdniana jako osobna klasa.
#' @item interactive TRUE, jeœli zmienna ma byæ dyskretyzowana interaktywnie. W
#'                   przeciwnym razie, co jest wartoœci¹ domyœln¹, dyskretyzacja
#'                   jest automatyczna.
#' @item from zamiast automatycznego dzielenia, mo¿na podaæ wartoœci przedzia³ów (from,to].
#' @item to   zamiast automatycznego dzielenia, mo¿na podaæ wartoœci przedzia³ów (from,to].
#' @item ... inne parametry do funkcji \link{\code{drzewo}}.
#' @seealso \code{drzewo}
#' @author Micha³ Danaj
#TODO wyci¹gn¹æ parametry drzewa
discretization<-function(x, y, special_val=NULL, min_bucket=20, max_gleb=3,
		interactive=FALSE,breaks=NULL,testy=FALSE,...){
	
	if(testy==TRUE) {
		print ('========   funckcja discretizetion   ==========')
		print('nulle w x:')
		print(table(is.na(x)))
		print('nulle w y:')
		print(table(is.na(y)))
		print ('========  ==========')
	}
	
	if (length(x)!=length(y))
		stop("discretization: parametry 'x' i 'y' maj¹ ró¿ne d³ugoœci!");
	
	#Mimo, ¿e przygotowywya³em funkcjê do obs³ugi null-i, to rezygnujê z tego
	#ze wzglêdów bezpieczeñstwa.
	if (any(is.na(x)) | any(is.na(y)))
		stop ("discretization: W 'x' ani 'y' nie mo¿e byæ NA!");
	
	# Wartoœci specjalne
	special_idx<-is.na(x)|x %in% special_val;
	special_val<-unique(x[special_idx & !is.na(x)]);
	#czy s¹ NA
	sa_na<-any(is.na(x));
	
	## Obs³ugujê wartoœci ci¹g³e
	#jeœli zosta³y podane zakresy przedzi³ów, to dzielimy wg nich
	if (!is.null(breaks)){
		bucket_drzewo<-buckety_stat2(breaks, x[!special_idx], y[!special_idx], total=FALSE);
	} else {
		if (!interactive){
			bucket_drzewo<-drzewo(x[!special_idx],y[!special_idx], min_bucket=min_bucket, max_gleb=max_gleb, n_buckets=5, wytnij=0.01, ...)
		}else{
			bucket_drzewo<-interactive_tree(score=x[!special_idx],def=y[!special_idx],
					span=0.80, min_split=200, min_bucket=min_bucket,
					buckets=60, max_gleb=2)
		}
	}
	
	#przypisujê nr przedzia³ów dla wartoœci ci¹g³ych
	bucket_drzewo$fitted<-bucket_drzewo$nr;
	classing<-rep(NA,length(x));
	print('przed przypisaniem')
	classing[!special_idx]<-przypisz(x[!special_idx], bucket_drzewo);
	print('po przypisaniu')
	#print(cbind(x,special_idx, x[!special_idx], classing))
	print(bucket_drzewo)
	
	#nadajê indeksy wartoœciom specjalnym i je przypisujê
	special_map<- -(length(special_val):1);
	names(special_map)<-special_val;
	classing[special_idx]<-special_map[as.character(x[special_idx])];
		
	#i jeszcze NA
	classing[is.na(x)]<- 0;
	
	
	print(table(classing,useNA='always'))
	#liczê statystyki
	classing_stat<-buckety_stat(classing, y, total=TRUE);
	
	#zmieniam nazwy wierszy, ¿eby nie by³y numery a labele klas
	#mapping<-c(names(special_map), "<NA>", rownames(bucket_drzewo), 'TOTAL');
	mapping<-c(names(special_map), "<NA>", rownames(bucket_drzewo), 'TOTAL');
	names(mapping)<-c(special_map,"0",bucket_drzewo$nr, 'TOTAL');
	
	rownames(classing_stat)<-mapping[rownames(classing_stat)];
	
	#które wartoœci s¹ specjalne (dyskretne)
	classing_stat$discret<-rep("",nrow(classing_stat));
	if (sa_na)
		classing_stat[c("<NA>", special_val),"discret"]<-rownames(classing_stat[c("<NA>", special_val),])
	else if (length(special_val)>0)
		classing_stat[as.character(special_val),"discret"]<-rownames(classing_stat[as.character(special_val),]);
	
	classing_stat$discret[nrow(classing_stat)]<-"<TOTAL>";
	
	classing_stat$do<-classing_stat$srodek<-classing_stat$od<-NA;
	classing_stat[classing_stat$discret=="",c("od","srodek","do")]<-bucket_drzewo[,c("od","srodek","do")];
	
	classing_stat<-classing_stat[,c('nr','label','discret', 'od','srodek','do','n_good','pct_good','n_bad','pct_bad','n_obs','pct_obs',
					'br','woe','logit')];
	if(testy==TRUE) 
		print('+++++++++ koniec discret +++++++++++')
	return(classing_stat);
}

####################


drzewo<-function(score, def, freq=NULL, wytnij=0, min_split=30, min_bucket=10, max_gleb=4, n_buckets=20, plot=TRUE, testy=FALSE,...)
{
	
	#musze teraz wywolac buckety_br, bo pozniej agreguje wektory
	#do jednego scoru i wykorzystuje freq, a buckety_br musi miec
	#jeden rekord=jedna obserwacja
	#if (plot)
	#	b<-buckety_br(score , def, n_buckets, method="eq_lengt");
	if (testy==TRUE){
		print('=================   funkcja drzewo  ======================')
		print(length(score))
		print(length(def))
	}
	k<-order(score);
	score<-score[k];
	def<-def[k];
	
	if (is.null(freq))
		freq<-rep(1,length(score))
	else
	{
		stop("poki co nie obslugiwane freq inne od NULL");
		freq<-freq[k];
	}
	if (wytnij>0)
	{
		usun<-usun_konce(score, prob=wytnij);
		if (length(usun)>0){
			score<-score[-usun];
			def<-def[-usun];
			freq<-freq[-usun];
		}
	}
	
	print('---xx')
	
	print(length(def))
	print(length(freq))
	print(length(score))
	
	def_a<-tapply(def,score,sum);
	freq_a<-tapply(freq, score, sum);
	score_a<-unique(score);
	
	
	print('---')
	
	print(length(def_a))
	print(length(freq_a))
	print(length(score_a))
	
	#Zabezpieczam siê przed zminiejszeniem precyzji liczb w wektorze score podczas konwersji
	#na zmienn¹ znakow¹, wykonywan¹ podczas wykonania funkcji tapply
	if (length(def_a)!=length(score_a)){
		score_a<-as.numeric(names(def_a))
		warning("W funkcji 'drzewo' wykonywana jest konwersja score na ci¹g znaków, która
			zmniejszy³a precyzjê. Dalsze dzia³anie funkcji przebiega w oparciu o zmniejszon¹
			precyzjê liczb. Zaleca siê zmniejszenie precyzji danych wejœciowych.")
		
	}

	
	#vec_stats(score_a);
	w<-drzewo_podzial(score_a, def_a, 1, min(score), max(score), freq_a, 0, min_split, min_bucket, max_gleb);
	
	#wybieram liscie
	w<-w[is.na(w$poprawa),];
	w<-w[order(w$od),];
	
	#i robie dla nich statystyki
	breaks<-sort(unique(c(w$od, w$do)));
	
	# jeœli jest tylko jeden liœæ
	if (length(breaks)==1)
		bucket<-buckety_stat(score, def, total=FALSE)
	else
		bucket<-buckety_stat(cut(score, breaks, include.lowest=TRUE), def, total=FALSE);
	
	#uzupe³niam statystyki
	bucket$fitted<-bucket$br;
	
	bucket$od<-w$od;
	bucket$do<-w$do;
	bucket$srodek<-(w$od+w$do)/2;
	#	----- rysowanie -----
	if (plot)
	{
		
		plot(bucket$srodek, bucket$br, xlab="", ylab="", xlim=c(range(breaks)),...);
		for (i in 1:length(w$od))
			lines(c(bucket$od[i], bucket$do[i]), c(bucket$fitted[i],bucket$fitted[i]), col="blue");
	}
	
	if (testy==TRUE){
		print('+++++++++++++++++++   koniec funkcja drzewo  +++++++++++++++++++++')
	}
	
	bucket
}

drzewo_plot<-function(liscie_drzewa,...){
	#liscie<-liscie_drzewa[liscie_drzewa$discret=="",];
	liscie<-liscie_drzewa[!is.na(liscie_drzewa$od),];
	
	# Jeœli s¹ jakieœ liœcie ci¹g³e
	if (nrow(liscie)>1){
		# Wartoœci nieskoñczone œrodka zamieniam krañcowymi
		liscie$srodek[liscie$srodek==-Inf]<-min(liscie$do)
		liscie$srodek[liscie$srodek==Inf]<-max(liscie$od)
		
		breaks<-sort(unique(c(liscie$od,liscie$do)));
		#usuwam nieskoñczonoœci
		niesk<-which(is.infinite(breaks))
		if (length(niesk)>0)
			breaks<-breaks[-niesk]
		
		skala <- sqrt(liscie$n_obs/(sum(liscie$n_obs)/nrow(liscie)));
		plot(liscie$srodek, liscie$br, xlim=c(range(breaks)),cex=skala,...);
		for (i in 1:length(liscie$od))
			lines(c(liscie$od[i], liscie$do[i]), c(liscie$br[i],liscie$br[i]), col="blue");
	}
	# W przeciwnym razie
	else{
		
		#usuwam totala
		liscie<-liscie_drzewa[!(is.na(liscie_drzewa$discret) | as.character(liscie_drzewa$discret)=='<TOTAL>'),];
		
		#jeœli coœ tam jest
		if(nrow(liscie)>0){
		
			skala <- sqrt(liscie$n_obs/(sum(liscie$n_obs)/nrow(liscie)));
			
			plot(liscie$br, cex=skala, axes=FALSE,...);
			axis(1, at=1:nrow(liscie), labels=liscie$discret);
			axis(2);
			box();
		}
		else{
			#rysujê pusty wykres
			plot(0,pch='x');
		}
	}
}

drzewo_podzial<-function(score, def, nr_wezla, od, do, freq, glebokosc,
		min_split=200, min_bucket=100, max_gleb=3, testy=FALSE)
{
	if (testy==TRUE){
		print("===============   funkcja drzewo_podzial   ====================");
		print(nr_wezla);
		print(length(score))
		print(length(def))
	}
	all_obs<-sum(freq);
	print(table(freq,useNA='always'))
	print(all_obs)
	all_bad<-sum(def);
	br_akt<-all_bad/all_obs;
	gini_akt<-br_akt*(1-br_akt)*all_obs;
	wezel<-data.frame(nr_wezla, rodzic=floor(nr_wezla/2), od, do, n_obs=all_obs, n_bad=all_bad, br=br_akt, poprawa=NA, podzial=NA);
	
	wynik<-wezel;
	
	#jesli ilosc obserwacji jest wystarczajaca, aby zrobic podzial w wezle
	#i jeszcze mo¿emy dorobiæ liœcie
	if (all_obs>min_split)
	{
		cum_bad_lewo<-cumsum(def);
		cum_obs_lewo<-cumsum(freq);
		
		cum_bad_prawo<-(all_bad-cum_bad_lewo);
		cum_obs_prawo<-(all_obs-cum_obs_lewo);
		
		br_lewo<-cum_bad_lewo/cum_obs_lewo;
		br_prawo<-cum_bad_prawo/cum_obs_prawo;
		
		gini_lewo<-br_lewo*(1-br_lewo)*cum_obs_lewo;
		gini_prawo<-br_prawo*(1-br_prawo)*cum_obs_prawo;
		
		gini_roz<-gini_akt-(gini_prawo+gini_lewo);
		#print("gini");print(gini_akt);print(gini_prawo);print(gini_lewo)
		
		#zostawiam podzialy, dla ktorych spelnione sa wymogi na ilosc
		#obserwacji w wynikowych lisciach
		zostaw<-(cum_obs_lewo>min_bucket)&(cum_obs_prawo>min_bucket);
		gini_roz[!zostaw]<-NA;
		
		#nr podzialu maksymalizujacego roznice gini
		nr<-which.max(gini_roz);
		if (length(nr)>0 & glebokosc<max_gleb)
		{
			wezel$poprawa<-gini_roz[nr];
			wezel$podzial<-(score[nr]+score[nr+1])/2;
			l<-length(score);
			print ('testy podzialu')
			print(length(score)) 
			print(length(def))
			print('length(freq):')
			print(length(freq))
			print('nr:') 
			print(nr)
			print('nic')
			wl<-drzewo_podzial(score[1:nr], def[1:nr], nr_wezla*2, od, wezel$podzial, freq[1:nr], glebokosc+1,
					min_split, min_bucket, max_gleb);
			wr<-drzewo_podzial(score[(nr+1):l], def[(nr+1):l], nr_wezla*2+1, wezel$podzial, do, freq[(nr+1):l], glebokosc+1,
					min_split, min_bucket, max_gleb);
			wynik<-rbind(wezel,wl,wr);
		}
	}
	
	if (testy==TRUE){
		print("===============   koniec funkcja drzewo_podzial   ====================");
	}
	wynik
}


genRaportBody<-function(wyniki, kolejnosc, dir, plik_main, scale){
	for (i in 1:length(wyniki)){
		#for (i in 1:1){
		wynik<-wyniki[[kolejnosc[i]]];
		nazwa_zmiennej<-names(wyniki)[kolejnosc[i]];
		
		cat(sprintf('<a name="%s">', nazwa_zmiennej), file=plik_main, append=TRUE)
		HTML.title(nazwa_zmiennej);
		cat('</a>', file=plik_main, append=TRUE)
		#    windows();
		
		### dyskryminacja ###
		HTML.title("Discrimination GINI", HR=3);
		if (!is.null(wynik$rozklady$pct_all_tbl)){
			do_wykresu<-melt(wynik$dyskryminacja)
			do_wykresu<-do_wykresu[do_wykresu$X2!='AR_calosc',];
			X1_order<-ordered(do_wykresu$X1, levels=rownames(wynik$dyskryminacja));
			print(str(do_wykresu));
			print(xyplot(value ~ X2 , group = X1_order, data=do_wykresu, type='b',
							xlab="Date", ylab="GINI", main=nazwa_zmiennej));
			print('1');
			HTMLplot(Caption = "", file = plik_main, append = TRUE, GraphDirectory = dir,   GraphFileName = paste(nazwa_zmiennej, ' discrimination'), GraphSaveAs = "png", GraphBorder = 1,  Align = "center",
					Width = 400, Height = 400, WidthHTML = NULL,     HeightHTML = NULL, GraphPointSize = 12, GraphBackGround = "white",     GraphRes = 72)
			
			HTML(wynik$dyskryminacja);
		}
		
		###   rysunek PIT/TTC   ###
		HTML.title("Point in Time or Through the Cycle", HR=3);
		if (!is.null(wynik$rozklady$avg_t_tbl)){
			plot(wynik$rozklady$avg_t_tbl['TOTAL',-ncol(wynik$rozklady$avg_t_tbl)], main="PIT/TTC",
					ylab="Mean LGD", xlab="Date");
			points(wynik$rozklady$estim, col="green")
			HTMLplot(Caption = "Does changes in variable distribution follow changes of portfolio LGD?",
					file = plik_main, append = TRUE, GraphDirectory = dir,   GraphFileName = paste(nazwa_zmiennej, 'cycle'), GraphSaveAs = "png", GraphBorder = 1,  Align = "center",
					Width = 400, Height = 400, WidthHTML = NULL,     HeightHTML = NULL, GraphPointSize = 12, GraphBackGround = "white",     GraphRes = 72)
			
			HTML(    t(data.frame("Portfolio LGD" = wynik$rozklady$avg_t_tbl['TOTAL',-ncol(wynik$rozklady$avg_t_tbl)],
									"Estimated LGD" = wynik$rozklady$estim))
			);
		}
		
		###   dyskretyzacja   ###
		HTML.title("Buckets", HR=3);
#    dev.off();
		windows(1400,700);
		par(mfrow=c(1,2));
		
		#nie wiem, czemu by³ tu wymóg rysowania tylko ci¹g³ych wartoœci
		#zobaczymy, jak to bêdzie po usuniêciu tego.
		#ciagle<-nchar(wynik$dyskretyzacja$discret)==0
		#drzewo_plot(wynik$dyskretyzacja[ciagle,], xlab=nazwa_zmiennej, ylab="Mean LGD",
		#		main=paste(nazwa_zmiennej,"discretization"));
		drzewo_plot(wynik$dyskretyzacja, xlab=nazwa_zmiennej, ylab="Mean LGD",
				main=paste(nazwa_zmiennej,"discretization"));
		ile_row<-nrow(wynik$dyskretyzacja);
		b<-barplot(wynik$dyskretyzacja$pct_obs[-ile_row], names.arg=rownames(wynik$dyskretyzacja)[-ile_row], xlab=nazwa_zmiennej,
				ylab='Distribution',main="Distribution with LGD");
		par(usr=c(par()$usr[1:2], scale))
		lines(b, wynik$dyskretyzacja$br[-ile_row],type="o", col="red", lty="solid", pch="x")
		axis(4)
		HTMLplot(Caption = "Results of discretization", file = plik_main, append = TRUE, GraphDirectory = dir,   GraphFileName = paste(nazwa_zmiennej, 'tree'), GraphSaveAs = "png", GraphBorder = 1,  Align = "center",
				Width = 800, Height = 400, WidthHTML = NULL,     HeightHTML = NULL, GraphPointSize = 12, GraphBackGround = "white",     GraphRes = 72)
		
		par(mfrow=c(1,1));
		HTML(wynik$dyskretyzacja);
		
		###   rozk³ady    ###
		HTML.title("Distribution of buckets", HR=3);
		if (!is.null(wynik$rozklady$pct_all_tbl)){
			do_wykresu<-melt(wynik$rozklady$pct_all_tbl)
			do_wykresu<-do_wykresu[do_wykresu$X1!='TOTAL' & do_wykresu$X2!='TOTAL',];
			X1_order<-ordered(do_wykresu$X1, levels=rownames(wynik$rozklady$pct_all_tbl));
			
			#dev.off();
			#png(filename = paste(dir,"xxx.png",sep="/"), width = 480, height = 480)
			print(barchart(value ~ X2|X1_order , stack=TRUE, data=do_wykresu, main=paste("Distribution of", nazwa_zmiennej),
							xlab='Date', ylab='Percent in given date'))
			#plot(1:10);
			#dev.off()
			
			#
			HTMLplot(Caption = "", file = plik_main, append = TRUE, GraphDirectory = dir,   GraphFileName = paste(nazwa_zmiennej, 'distribution'),
					GraphSaveAs = "png", GraphBorder = 1,  Align = "center",
					Width = 800, Height = 400, WidthHTML = NULL,     HeightHTML = NULL, GraphPointSize = 12, GraphBackGround = "white",     GraphRes = 72)
			
			HTML(wynik$rozklady$obs_all_tbl, caption="Number of observations");
			HTML(wynik$rozklady$pct_all_tbl, caption="% share at given date");
			
			#     œredni LGD    #
			HTML.title("Mean LGD", HR=3);
			do_wykresu<-melt(wynik$rozklady$avg_t_tbl)
			do_wykresu<-do_wykresu[do_wykresu$X1!='TOTAL' & do_wykresu$X2!='TOTAL',];
			X1_order<-ordered(do_wykresu$X1, levels=rownames(wynik$rozklady$pct_all_tbl));
			print(xyplot(value ~ X2|X1_order , data=do_wykresu, type='b', xlab="Date", ylab="Mean LGD", main=nazwa_zmiennej,
							strip=strip.custom(bg='green')));
			
			HTMLplot(Caption = "", file = plik_main, append = TRUE, GraphDirectory = dir,   GraphFileName = paste(nazwa_zmiennej, 'LGD by bucket'), GraphSaveAs = "png", GraphBorder = 1,  Align = "center",
					Width = 800, Height = 400, WidthHTML = NULL,     HeightHTML = NULL, GraphPointSize = 12, GraphBackGround = "white",     GraphRes = 72)
			print(xyplot(value ~ X1_order |X2 , data=do_wykresu, type='b', xlab="Bucket", ylab="LGD", main=nazwa_zmiennej));
			HTMLplot(Caption = "", file = plik_main, append = TRUE, GraphDirectory = dir,   GraphFileName = paste(nazwa_zmiennej, 'LGD by time'), GraphSaveAs = "png", GraphBorder = 1,  Align = "center",
					Width = 800, Height = 400, WidthHTML = NULL,     HeightHTML = NULL, GraphPointSize = 12, GraphBackGround = "white",     GraphRes = 72)
			
			HTML(wynik$rozklady$avg_t_tbl, caption="Mean LGD");
			
			dev.off();
		}
		HTML("<HR><HR><HR><HR><HR>")
	}
}



###########   generujê plik z menu   #################

genRaportMenu<-function(wyniki, dir){

	plik_menu<-HTMLInitFile(dir, 'univariate_menu');
	
	GINI<-sapply(wyniki, function(wynik){
				gini<-wynik$dyskryminacja[1,'AR_calosc']
				if (is.null(gini))
					return(NA);
				return(gini);
			})
		
	
	#GINI
	kolej<-rev(order(GINI));
	HTML.title('Sortowanie po GINI');
	for (i in 1:length(wyniki)){
		nazwa_zmiennej<-names(wyniki)[kolej[i]];
		cat(sprintf('<a href="univariate_main.html#%s" target=main>%s (%f)</a></br>\n',nazwa_zmiennej,nazwa_zmiennej, round(GINI[kolej[i]],3))
				, file=plik_menu, append=TRUE);
	}
	
	
	
	
	#alfabetycznie
	kolejnosc<-order(names(wyniki));
	HTML.title('Sortowanie Alfabetyczne');
	for (i in 1:length(wyniki)){
		nazwa_zmiennej<-names(wyniki)[kolejnosc[i]];
		cat(sprintf('<a href="univariate_main.html#%s" target=main>%s</a></br>\n',nazwa_zmiennej,nazwa_zmiennej)
				, file=plik_menu, append=TRUE);
	}
	
	HTMLEndFile(plik_menu)
}

#TODO wyci¹gn¹æ scale!
#' 
#' @param wyniki - lista z wynikami dyskretyzacji itp, z funkcji \link{\code{univariate_anal_stats}} 
#' @param kolejnosc - kolejnoœæ wg której zmienne maj¹ byæ wyœwietlone. 
#' @param dir - katalog z raportem, jako pe³na bezwzglêdna œcie¿ka! Katalog musi byæ stworzony.  
#' 
#' @author Micha³ Danaj
genRaport<-function(wyniki, dir, kolejnosc=1:length(wyniki), scale=c(0,0.2)){
	
	makeCSSFile(dir)
	
	HTMLStart(dir , "univariate", HTMLframe=TRUE, Title="Univariate analysis",
			echo=TRUE);
	HTMLStop();
	
	
	plik_main<-HTMLInitFile(dir, 'univariate_main', CSSFile='R2HTML MD.css');
	
	genRaportBody(wyniki, kolejnosc, dir, plik_main, scale)
	genRaportMenu(wyniki, dir)
	
}

#zaci¹ga z plików za zapisanymi wyliczonymi zmiennycmi podane zmienne
#z okreœlonych dat.

getVariables<-function(variables, dates){
	#wczytujê nazwy zmiennych z informacj¹, w którym pliku siê znajduj¹
	nazwy_zmiennych<-read.delim("budowa_willcc_nazwy_zm2.txt", col.names=c('variable_name', 'file'), as.is=TRUE);
	#ograniczam siê do zaci¹ganych zmiennych
	nazwy_zmiennych<-nazwy_zmiennych[nazwy_zmiennych$variable_name %in% variables,];
	dane_all<-NULL;
	
	for (i in unique(nazwy_zmiennych$file)){
		#for (i in 9:9){
		#for (i in 2:2){
		print(paste("Plik nr ",i));
		
		
		baza<-sprintf('budowa_willcc_variables_%02.0f.db',i);
		sqldf(sprintf("attach '%s' as new", baza));
		
		dates_str<-paste(dates, collapse="','");
		dates_str<-paste("('",dates_str,"')", sep='');
		variable_names_str<-paste(variables, collapse="','");
		variable_names_str<-paste("('",variable_names_str,"')", sep='');
		
		sql<-sprintf("select distinct account_contract, VARIABLE_NAME, VAL, REPORTINGDATE from budowa_willcc_variables_%02.0f
						where REPORTINGDATE in %s and variable_name in %s order by account_contract, reportingdate;",i, dates_str,
				variable_names_str);
		dane_sql<-sqldf(sql, dbname=baza);
		names(dane_sql)<-tolower(names(dane_sql));
		
		print("Wczytano dane, wykonujê przekszta³cenia");
		dane_trans<-cast(dane_sql, account_contract + reportingdate~ variable_name, value='val');
		
		if (is.null(dane_all))
			dane_all<-dane_trans
		else
		{
			#test kolejnoœci
			if (any(dane_trans$account_contract!=dane_all$account_contract ||
							dane_trans$reportingdate!=dane_all$reportingdate))
				stop("Uwaga! Kolejnoœæ account_contract lub reportingdate jest niezgodna!");
			tmp<-data.frame(dane_trans[,c(-1,-2)]);
			names(tmp)<-names(dane_trans)[-c(1,2)];
			dane_all<-cbind(dane_all, tmp);
		}
		
		rm(dane_sql);
		gc();
		gc(reset=TRUE);
		gc(reset=TRUE);
		gc(reset=TRUE);
		gc(reset=TRUE);
		gc(reset=TRUE);
	}
	return(dane_all);
}

################################################################################
#                    zmienione funkcjie z pakietu binom                        #
################################################################################





#' Interaktywny podzia³ zmiennej ci¹g³ej na buckety
#'
#' Interaktywny podzia³ zmiennej ci¹g³ej na buckety. Uwaga! W danych
#' nie mo¿e byæ wartoœci NULL. Nieopisane zmienne s¹ to zmienne z
#' wykorzystywanej funkcji \link{\code{drzewo}} oraz \link{\code{reg_nieparam}}.
#' Po wykonaniu zmiany wyœwietlane s¹ statystyki po bucketach oraz statystyka AR.
#' @item score zmienna score'owa.
#' @item def zmienna odpowiedzi z zakresu [0,1]. Np. default, LGD.
#' @seealso \link{\code{drzewo}}, \link{\code{quick_AR}},  \link{\code{buckety_stat2}}.
#' @return \code{data.frame} ze statystykami.
#' @author Micha³ Danaj
interactive_tree<-function(score, def, span=0.80, min_split=200, min_bucket=100,
		buckets=60, max_gleb=2)
{
	
	if(any(is.na(score)|is.na(def)))
		stop('interactive_tree: Niedozwolone wartoœci NULL!')
	#wylicza pozycjê punktów symuluj¹cych menu i je rysuje
	punkty_menu<-function(){
		
		minx<-min(axis(1))
		maxx<-max(axis(1))
		
		miny<-min(axis(2))
		maxy<-max(axis(2))
		
		deltax=(maxx-minx)/20;
		deltay=(maxy-miny)/20;
		
		x_fun=c(maxx-deltax,maxx);
		y_fun=c(maxy,maxy);
		
		points(x_fun,y_fun, col=c("green","red"), pch=19);
		
		return(identify(x_fun,y_fun,n=1));
	}
	
	#posortuj wg score
	kolej<-order(score);
	score<-score[kolej];
	def<-def[kolej];
	
	#par(mfrow=c(1,2), xpd=F);
	reg<-reg_nieparam(score=score , default=def, buckets=buckets, wytnij=0.01, span=span, degree=2);
	drz<-drzewo(score, def, wytnij=0.01, min_split=min_split, min_bucket=min_bucket, max_gleb=max_gleb);
	reg<-reg_nieparam(score , def, buckets=buckets, wytnij=0.01, span=span, degree=2);
	for (i in 1:length(drz$od))
		lines(c(drz$od[i], drz$do[i]), c(drz$br[i],drz$br[i]), col="blue");
	
	drz_pocz<-drz;
	
	#par(xpd=NA);
	#bar<-barplot(drz$br);
	#text(bar,drz$br,paste(drz$n_default, "\n", drz$n_obs, "\n", round(drz$br*100,1)), cex=0.8, adj=c(1);
	drz<-buckety_stat2(c(drz$od,drz$do), score, def, total=FALSE);
	points(drz$median, drz$br, col="red",pch="X");
	#points(drz$srodek, drz$br, col="red",pch="X");
	#points(drz$n_bad/drz$n_obs, drz$br, col="green",pch="I");
	
	print(drz[,c('od', 'do','n_bad','n_obs','pct_obs','br','woe')]);
	kolej<-order(-drz$br);
	print(paste('AR =',with(drz, AR_quick(n_bad[kolej], n_obs[kolej]))));
	
	nr<-punkty_menu();
	
	#print(data.frame(drz$od, drz$do, drz$n_default,drz$n_obs, drz$br));
	
	while (length(nr)>0)
	{
		#polacz
		if (nr==1)
		{
			#text(x_fun[1]-deltax, y_fun[1]-deltay,"Wybierz dwa przedzialy",cex=0.7);
#  		id<-identify(drz$srodek,drz$br,n=2, pos=1, col="black")$ind;
			id<-identify(drz$median,drz$br,n=2, pos=1, col="black")$ind;
			#print(id);
			if (length(id)==2)
			{
				id<-id[order(id)];
				drz[id[1],"do"]<-drz[id[2],"do"];
				drz_old<-drz;
				drz<-buckety_stat2(c(drz$od[-id[2]], drz$do[-id[2]]), score, def, total=FALSE);
				drz$od<-drz_old$od[-id[2]];
				drz$do<-drz_old$do[-id[2]];
				#drz<-buckety_stat(score, def, drz$od[-id[2]], drz$do[-id[2]]);
				#print(data.frame(drz$od, drz$do, drz$n_default,drz$n_obs, drz$br));
				drz$fitted<-drz$br;
			}
		}
		#rozlacz
		if (nr==2)
		{
			#text(x_fun[1]-deltax, y_fun[1]-deltay,"Wybierz jeden przedzial",cex=0.7);
			id<-identify(drz$median,drz$br,n=1, pos=1, col="black")$ind;
			
			if (length(id)==1)
			{
				wybr<-score>drz[id,"od"] & score<=drz[id,"do"];
				
				def_a<-tapply(def[wybr],score[wybr],sum);
				freq_a<-tapply(rep(1,length(score[wybr])), score[wybr], sum);
				score_a<-unique(score[wybr]);
				#vec_stats(score_a);
				
				wl<-drzewo_podzial(score_a, def_a, 1, drz[id,"od"], drz[id,"do"], freq_a, 0, min_split, min_bucket, 1);
				wl<-wl[is.na(wl$podzial),];
				
				drz<-drz[-id,];
				od<-c(drz$od,wl$od);
				od<-od[order(od)];
				do<-c(drz$do, wl$do);
				do<-do[order(do)];
				#i robie dla nich statystyki
				breaks<-sort(unique(c(od, do)));
				
				# jeœli jest tylko jeden liœæ
				if (length(breaks)==1)
				{drz<-buckety_stat2(c(od, do), score, def, total=FALSE);
				}else
					drz<-buckety_stat2(c(od, do), score, def, total=FALSE);
				
				drz$fitted<-drz$br;
				
				#print(data.frame(drz$od, drz$do, drz$n_default,drz$n_obs, drz$br));
			}
		}
		reg_nieparam(score , def, buckets=buckets, wytnij=0.01, span=span, degree=2);
		for (i in 1:length(drz$od))
			lines(c(drz$od[i], drz$do[i]), c(drz$br[i],drz$br[i]), col="blue");
		
		points(drz$median, drz$br, col="red",pch="X");
		points(drz$mean, drz$br, col="green",pch="I");
		
		print(drz[,c('od', 'do','n_bad','n_obs','pct_obs','br','woe')]);
		kolej<-order(-drz$br);
		print(paste('AR =',with(drz, AR_quick(n_bad[kolej], n_obs[kolej]))));		
		nr<-punkty_menu();
	}
	
	return(drz);
}




#' £¹czy ze sob¹ buckety
#' 
#' £¹czy buckety. W przypadku bucketów z przedzia³ami zmiennej ci¹g³ej, mo¿liwe jest 
#' po³¹czenie tylko
#' przedzia³ów przlegaj¹cych do siebie. Dla nowo powsta³ych bucketów wylicza statystyki. 
#' W przypadku ³¹czenia przedzia³ów zmiennej ci¹g³ej, wiersze z tymi bucketami zostan¹ ze sob¹
#' po³¹czone i powstanie \{data.frame} z liczb¹ wierszy o jeden mniejsz¹ ni¿ w \code{bucket}.
#' Przy ³¹czeniu bucketów dyskretnych lub dyskretnego i ci¹g³ego, wiersze nie zostan¹ usuniête. 
#' Zostanie im nadany wspólny label oraz wspólny numer. Jeœli liczba bucketów do po³¹czenia jest<=1,
#' zwracany jest wejœciowy podzia³. 
#' @param x zmienna score'owa.
#' @param y zmienna odpowiedzi z zakresu [0,1]. Np. default, LGD.
#' @param buckets buckety.
#' @param nr1 numer wiersza w \{buckets} z pierwszym bucketem do po³¹czenia.
#' @param nr2 numer wiersza w \{buckets} z pierwszym bucketem do po³¹czenia.
#' @param new_label label jaki bêdzie nadany po³¹czonym bucketom. W przypadku braku podania, zostanie zatosowany domyœlny.
#' @returnType \code{data.frame}.
#' @author Micha³ Danaj
polacz_buckety<-function(x, y, buckets, row_idxs, new_label=NULL)
{
	row_idxs<-sort(unique(row_idxs));
	#sprawdzam, czy jest co ³¹czyæ
	if (length(row_idxs)<=1)
		return(buckets)
	#sprawdzam, czy nie wyszliœmy poza zakres
	if(min(row_idxs)<1 | max(row_idxs)>nrow(buckets)-1) 
		stop("Numery wierszy s¹ poza zakresem zmiennej 'buckets'");
	#sprawdzam, czy nie ma ju¿ takiego labela w innych wierszach
	if (!is.null(new_label) & any( buckets$label[-row_idxs]==new_label))
		warning("Podana wartoœæ 'new_label' znajduje siê ju¿ w 'buckets' w wierszu innym ni¿ aktualnie ³¹czone wiersze.")
	for (i in 1:(length(row_idxs)-1)){
		nr1<-row_idxs[i];
		nr2<-row_idxs[i+1];
		#jeœli ³¹czymy dane ci¹g³e
		if (all((!is.na(buckets$od) & !is.na(buckets$do))[c(nr1,nr1)])){
			if (nr2-nr1!=1)
				stop("B³êdnie podane wiersze do po³¹czenia! Przedzia³y powinny byæ do siebie przyleg³e!")
			#jeœli label nie podany, to go tworzê
			if (is.null(new_label)){
				new_label<-strsplit(buckets$label[nr1],',')[[1]][1]
				new_label<-paste(new_label,strsplit(buckets$label[nr2],',')[[1]][2],sep=',')
			}
			#jeszcze przepisujê krañce przedzia³ów i inne wartoœci nie wyliczane w buckety_stat a
			#wyliczane w buckety_stat2
			buckets$do[nr1]<-buckets$do[nr2];
			buckets$srodek[nr1]<-c(buckets$od[nr1]+buckets$do[nr1])/2;
			buckets$label[nr1]<-new_label;
			rownames(buckets)[nr1]<-new_label;
			#usuwam drugi bucket
			buckets<-buckets[-nr2,];
			
		}else{
			if (is.null(new_label))
				new_label<-paste(buckets$label[nr1],buckets$label[nr2],sep=',');
			buckets$label[c(nr1,nr2)]<-new_label;
			buckets$nr[nr2]<-buckets$nr[nr1];
		}
		buckets$fitted<-buckets$label;
		x_buckets<-przypisz2(x,buckets);
		buckets_new<-buckety_stat(x_buckets, y)
		buckets<-cbind(buckets[,c('nr','label','discret','od','srodek','do')], 
				buckets_new[buckets$label,c('n_good','pct_good','n_bad','pct_bad','n_obs','pct_obs','br','woe','logit')]);
	}	
	return(buckets);
}

#' Wylicza statystyki dla zmiennej score'owej i objaœnianej wg zadanych bucketów
#'
#' Wylicza statystyki i zwraca je w postaci listy. Do wyliczenia \code{AR} u¿ywa jako zmiennej
#' score'owej oryginalnych wartoœci BR przypisanych do bucketu. Tzn, mimo ¿e na nowych danych kolejnoœæ
#' bucketów, sortuj¹c je po br mo¿e byæ inna ni¿ oryginalnie, to stosowana jest oryginalna kolejnoœæ.
#' @tiem buckets \code{data.frame} z podzia³em zmiennej na buckety.
#' @item x zmienna, po której procedura bêdzie sortowaæ.
#' @item y zmienna odpowiedzi.
#' @item czas podzia³ na okresy czasowe.
#' @item proby \code{data.frame}, w którym ka¿da kolumna jest wektorem logicznym,
#'       zawieraj¹cym informacje, czy obserwacja nale¿y do danej próby.
#' @author Micha³ Danaj
univariate_stats_new_data<-function(buckets,x,y,czas,proby=rep(TRUE, length(y))){
	buckets$fitted<-buckets$label;
	nowe_wartosci<-przypisz2(x, buckets)
	buckets_new<-univariate_anal_stats(nowe_wartosci, y, czas=czas, proby=proby)
	stat1<-cbind(buckets[,c('nr','label','discret','od','srodek','do')], 
			buckets_new$dyskretyzacja[buckets$label,c('n_good','pct_good','n_bad','pct_bad','n_obs','pct_obs','br','woe','logit')],
			br_orig=buckets$br);

	stat2<-NULL;
	stat3<-NULL;
	#Dalsze statystyki robiê pod warunkiem, ¿e jest wiêcej ni¿ jedna wartoœæ dyskretna
	#(W stat1 jest te¿ <TOTAL>, dlatego 2)
	if (nrow(stat1)>2){
		# przypisujê nazwê bucketu
		stat1$fitted<-stat1$label;
		x_discr<-przypisz2(x,stat1);
		
		# przypisujê BR
		stat1$fitted<-stat1$br;
		BR_discr<-przypisz2(x,stat1);
		
		# wyliczam drugie statystyki
		stat2<-univariate_anal_stats2(x_discr, y, czas, BR_discr);
		
		# wyliczam trzecie statystyki (GINI) po zadanych próbach i czasie
		buckets$fitted<-buckets$br;
		stare_BR<-przypisz2(x,buckets);
		stat3<-univariate_anal_stats3(-stare_BR, y, czas, proby);
		attr(stat3,'comment')<-"AR policzony w oparciu o oryginalne wartoœci BR."
	}
	return(list(dyskretyzacja=stat1, rozklady=stat2, dyskryminacja=stat3));
}



#' Generuje tabelê z coarse classing do dokumentacji
#' 
#' Jeœli Woe lub IV wychoch¹ +-Inf, to wstawia zamiast tego NA
#' @item wyniki lista z wynikami dyskretyzacji
makeCoarseClassingTables<-function(wyniki){
	lapply(wyniki, function(z){
				x<-z$dyskretyzacja;
				last_row<-nrow(x);
				x$GB_ODDS<-x$n_good/x$n_bad;
				gb_odds_total<-x$GB_ODDS[last_row];
				
				temp<-x$GB_ODDS/gb_odds_total;
				x$GB_index<-temp;
				
				print(x$GB_ODDS)
				x$GB_index[na.subst(x$GB_ODDS>gb_odds_total,FALSE)]<-paste(round(temp[na.subst(x$GB_ODDS>gb_odds_total,FALSE)]*100), 'G', sep='')
				x$GB_index[na.subst(x$GB_ODDS<=gb_odds_total, FALSE)]<-paste(round(temp[na.subst(x$GB_ODDS<=gb_odds_total, FALSE)]*100), 'B', sep='')
				
				x$IV<-(x$pct_good-x$pct_bad)*x$woe;
				
				#zamienia nieskoñczonoœci na braki danych
				x$IV[abs(x$IV)==Inf]<-NA;
				x$woe[abs(x$woe)==Inf]<-NA;
				
				x$IV[last_row]<-sum(x$IV[-last_row], na.rm=TRUE);
				x<-x[,c('label','n_obs','pct_obs','n_good','pct_good','n_bad','pct_bad','GB_ODDS','GB_index',
								'br','woe', 'IV')]
				x$discret[x$discret=='']<-x$label[x$discret==''];
				names(x)<-c('Coarse Classes','# Applicants','% Applicants',
						'# Good','% Good','# Bad','% Bad','GB Odds','GB Index',
						'LGD','WoE', 'IV')
				x<-unique(x);
				return(x);
			})
}


	
#' Wylicza punkty score'owe na podstawie parametrów modelu
#' @param model model, którego parametry zostan¹ przekszta³cone na punkty score'owe.
#' @param from od jakiej wartoœci maj¹ siê rozpoczynaæ punkty score'owe.
#' @param to do jakiej wielkoœci maj¹ byæ wartoœci punktów socre'owych.
#' @param test jeœli TRUE, zwraca testy istotnoœci atrybutów.
getScoreCard<-function(model, from, to, test=FALSE){
	wynik<-melt(model$xlevels);
	wynik$by<-paste(wynik$L1,wynik$value,  sep='')
	wynik$coeff<- -1*coef(model)[wynik$by]
	if (test)
		wynik$test<- (summary(model)$coefficients[,4])[wynik$by]
	wynik$coeff[is.na(wynik$coeff)]<-0;
	
	rng<-tapply(wynik$coeff, wynik$L1, function(x)max(x)-min(x))
	max_rng<-max(rng);
	
	wynik$points<-round((wynik$coeff-ave(wynik$coeff,wynik$L1, FUN=min))/max_rng*(to-from)+from)
	if (test==FALSE){
		names(wynik)<-c('value','variable','x','coeff', 'points');
		wynik<-wynik[,c('variable','value','coeff', 'points')];
	}
	else{
		names(wynik)<-c('value','variable','x','coeff', 'test', 'points');
		wynik<-wynik[,c('variable','value','coeff', 'test', 'points')];		
	}
	
	return(wynik);
}

#' Wylicza punkty score'owe na podstawie parametrów modelu
#' @param model model, którego parametry zostan¹ przekszta³cone na punkty score'owe.
#' @param from od jakiej wartoœci maj¹ siê rozpoczynaæ punkty score'owe.
#' @param to do jakiej wielkoœci maj¹ byæ wartoœci punktów socre'owych.
#' @param test jeœli TRUE, zwraca testy istotnoœci atrybutów.
getScoreCard2<-function(model, from, to, test=FALSE){
	wynik<-melt(model$xlevels);
	wynik$by<-paste(wynik$L1,wynik$value,  sep='')
	wynik$coeff<- -1*coef(model)[wynik$by]
	if (test)
		wynik$test<- (summary(model)$coefficients[,4])[wynik$by]
	wynik$coeff[is.na(wynik$coeff)]<-0;
	
	sum_min_est<-sum(tapply(wynik$coeff, wynik$L1, min))
	sum_max_est<-sum(tapply(wynik$coeff, wynik$L1, max))
	
	wynik$points<-round((wynik$coeff-sum_min_est)/(sum_max_est-sum_min_est)*(to-from)+from)
	if (test==FALSE){
		names(wynik)<-c('value','variable','x','coeff', 'points');
		wynik<-wynik[,c('variable','value','coeff', 'points')];
	}
	else{
		names(wynik)<-c('value','variable','x','coeff', 'test', 'points');
		wynik<-wynik[,c('variable','value','coeff', 'test', 'points')];		
	}
	
	return(wynik);
}


#' Przypisuje score z definicji karty scoringowej.
#' @param scoreCard definicja karty scoringowej. Patrz \link{\code{getScores}}.
#' @param x \code{data.frame} z kolumnami o nazwach takich, jak w definicji karty
#' @param sufix ci¹g znaków dodany do nazw wynikowych kolumn.
#' @seealso \link{\code{getScores}}.
#' scoringowej.
assignScore<-function(scoreCard, x, sufix='_points'){
	if (!is.data.frame(x))
		stop('"x" powinno byæ typu data.frame.');
	
	czy_ok<-c('variable','value','coeff', 'points')%in%names(scoreCard);
	if(all(czy_ok)==FALSE)
		stop(paste("W scoreCard brak kolumn(y)", c('variable','value','coeff', 'points')[!czy_ok]));
	
	#if (class(scoreCard)!='scoreCard')
	#	stop('"scoreCard" powinno byæ klasy scoreCard (patrz getScores).')
	wynik<-data.frame();
	nazwy<-unique(scoreCard$variable);
	nazwy2<-character(); #nazwy, które s¹ i w x i w scoreCard
	for (i in 1:length(nazwy)){
		nazwa<-nazwy[i];
		czesc<-scoreCard[scoreCard$variable==nazwa,];
		if (nazwa %in% names(x)==FALSE){
			warning(paste("W x brak kolumny",nazwa));
			next;
		}
		nazwy2<-c(nazwy2,nazwa);
		rownames(czesc)<-czesc$value;
		if (nrow(wynik)==0)
			wynik<-data.frame(czesc[x[,nazwa],'points'])
		else
			wynik[,nazwa]<-czesc[x[,nazwa],'points'];
	}	
	
	if(any(is.na(wynik))){
		warning('W x wyst¹pi³y wartoœci nie zdefiniowane w karcie scoringowej! W te miejsca przypisano NA.')
	}
	names(wynik)<-paste(nazwy2,sufix, sep='');
	wynik[,'score']=rowSums(wynik);	
	wynik
}


interakcja<-function(zm1, zm2, target){
	dwa<-tapply(target, list(zm1,zm2), mean)
	jeden<-log(dwa/(1-dwa))
	#jeden<-jeden[-grep('-999', rownames(jeden)),];
	#jeden<-jeden[,-grep('-999', colnames(jeden))];
	mfrow_oryg<-par('mfrow')
	par(mfrow=c(2,1))
	matplot(dwa, type=rep('b', ncol(jeden)))
	matplot(jeden, type=rep('b', ncol(jeden)))
	par(mfrow=mfrow_oryg)
	dwa
}



###############################################################################
#####																	  #####
###############################################################################


makeCSSFile<-function(dir){
tekst<-'
body {
	background: #FFFFFF;
			color: #000000;
			font-family: Verdana, Arial, Helvetica, sans-serif;
	font-size: 10pt;
	font-weight: normal
}

.tablesort {
	cursor: pointer;
	behavior: url(tablesort.htc);
	-moz-binding: url(tablesort.htc);
}

H1 {
	font-family: Arial, Helvetica, sans-serif;
	font-size: 30pt;
	font-style: normal;
	font-weight: bold;
	color: #3333CC;
			background: #004080;
			text-align: center;
	margin: 10pt 2.5%
}

H2 {
	font-family: Arial, Helvetica, sans-serif;
	font-size: 17pt;
	font-style: normal;
	font-weight: bold;
	color: #FFFFFF;
			background: #0050d0;
			text-align: center
}

H2.index {
	font-family: Arial, Helvetica, sans-serif;
	font-size: 17pt;
	font-style: normal;
	font-weight: normal;
	color: #FFFFFF;
			background: #0050d0;
			text-align: center;
	margin: 10pt 5%
}

H3 {
	font-family: Arial, Helvetica, sans-serif;
	font-size: 14pt;
	font-style: normal;
	font-weight: bold;
	text-align: center;
	color: #004080
}

H4 {
	font-family: T, Helvetica, sans-serif;
	font-size: 10pt;
	font-style: normal;
	font-weight: bold;
	color: #000000;
			line-height: 16pt
}

LI {
	font-family: Verdana, Arial, Helvetica, sans-serif;
	font-size: 10pt
}

A {
	font-family: Verdana, Arial, Helvetica, sans-serif;
	font-size: 10pt;
	text-decoration: none
}

.caption {
	font-style: italic
}

.title2 {
	font-family: Arial, Helvetica, sans-serif;
	font-size: 14pt;
	font-style: normal;
	font-weight: bold;
	color: #004080
}

.equation{
	font-weight: bold;
}

.command {
	font-family=verdana, arial;
	color=red
}

.partitle {
	font-family=verdana, arial;
	font-weight: bold
}

XMP {
	font-family: Verdana, Arial, Helvetica, sans-serif;
	font-size: 10pt
}

.function {
	font-family=courier;
	color=blue;
	font-size: 10pt
}

TR {
	font-family: Arial, Helvetica, Times, Helvetica, sans-serif;
	font-size: 10pt;
	font-style: normal;
	padding: 0 0
}

TR.firstline {
	color: #FFFFFF;
			background: #000000;
			text-align=center;
	font-weight: bold
}
TR.ListBackTitle {
	color: #FFFFFF;
			background: #000000;
			text-align=left;
	font-weight: bold
}
TD {
	background=#FFFFFF;
			padding: 0 0
}
TD.ListBackMain {
	background: #E0E0E0;
			padding: 0 0
}
TD.firstcolumn {
	padding: 5 10;
	background: #C0C0C0;
			text-align=right
}
TD.cellinside {
	padding: 5 10;
	background: #FFFFFF;
			text-align=right
}
/* CORRELATION MATRIX TRAFFIC HIGHLIGHT*/
		TD.corvarname {
			background-color="#FFFFFF";
			color=black;
			height: 1.1cm;
			text-align: right;
			font-weight: bold
		}
TD.corsep {
	width: 0.5cm
}
TD.cordiag {
	background-color=#fffff;
			color=white
}
TD.cor0 {
	background-color="#FFFFFF";
	color=black;
	width: 1.1cm;
	height: 1.1cm;
	text-align: center
}
TD.cor1 {
	background-color="#E6E6E6";
	color=black;
	width: 1.1cm;
	height: 1.1cm;
	text-align: center
}
TD.cor2 {
	background-color="#CCCCCC";
	color=black;
	width: 1.1cm;
	height: 1.1cm;
	text-align: center
}
TD.cor3 {
	background-color="#B3B3B3";
	color=black;
	width: 1.1cm;
	height: 1.1cm;
	text-align: center
}
TD.cor4 {
	background-color="#999999";
	color=black;
	width: 1.1cm;
	height: 1.1cm;
	text-align: center
}
TD.cor5 {
	background-color="#808080";
	color=white;
	width: 1.1cm;
	height: 1.1cm;
	text-align: center
}
TD.cor6 {
	background-color="#666666";
	color=white;
	width: 1.1cm;
	height: 1.1cm;
	text-align: center
}
TD.cor7 {
	background-color="#4D4D4D";
	color=white;
	width: 1.1cm;
	height: 1.1cm;
	text-align: center
}
TD.cor8 {
	background-color="#333333";
	color=white;
	width: 1.1cm;
	height: 1.1cm;
	text-align: center
}
TD.cor9 {
	background-color="#1A1A1A";
	color=yellow;
	width: 1.1cm;
	height: 1.1cm;
	text-align: center
}
TD.cor10 {
	background-color="#000000";
	color=yellow;
	width: 1.1cm;
	height: 1.1cm;
	text-align: center
}'

	nazwa<-paste(dir,'/R2HTML MD.css', sep='')
	file.create(nazwa)
	cat(tekst, file=nazwa)
}

		
#' Na podstawie zadanych warunków przypisuje etykiety
#' 
#' Na podstawie zadanych warunków w \code{mapping$war} przypisuje etykietê \code{mapping$label}.
#' W \code{mapping$label} nie mo¿e byæ wartoœci pustych stringów. Pusty string stosowany jest jako
#' brak danych. Jeœli dla jakieœ elementu ¿aden warunek nie zachodzi, zwracany jest pusty string.    
#'
#' Mo¿liwe jest okreœlenie warunku \code{else} poprzez wpisanie w \code{mapping$war} stringa "else".
#'
#' @param data \code{data.frame} do którego bêdzie przypisywane mapowanie
#' @param maping   \code{data.frame} z dwoma kolumnami znakowymi (lub do przerobienia przez \code{as.character})
#'			  \code{war} oraz \code{label}
#' @author Micha³ Danaj
mapuj<-function(data, mapping){
	
	
	#jeœli nie znakowe, to przerabiam na znakowe
	if (!is.character((mapping$war)))
		mapping$war<-as.character(mapping$war)
	if (!is.character((mapping$label)))
		mapping$label<-as.character(mapping$label)	
	if (!is.data.frame(data))
		stop("Funkcja mapuj: Dane musz¹ byæ typu data.frame, z kolumnami wykorzystywanymi w warunkach mapuj¹cych!")
	
	wynik<-rep("", nrow(data))
	
	gdzie_else<-which(mapping$war=='else')
	
	#jeœli jest else, to go wydzielam
	if (length(gdzie_else>0)){
		lab_else<-mapping$label[gdzie_else]
		mapping<-mapping[-gdzie_else,]
	}
	
	text=sprintf("with(data,%s)",mapping$war)	
	for (i in 1:nrow(mapping)){	
		
		war<-eval(parse(text=text[i]))		
		wynik[war]<-mapping$label[i]
	}
	
	#jeœli jest else, to go stosujê
	if (length(gdzie_else)>0)
		wynik[wynik==""]<-lab_else
	
	#sprawdzam, czy coœ siê nie przypisa³o
	if (any(wynik==""))
		warning("mapuj: Nie wszystkie wartoœci zosta³y przypisane.")
	
	return(wynik)
}

univariate_anal_stats4<-function(dane, mapowanie, czas=lastDay(dane$reportingdate, unit = "quater"),...){
	dyskretne<-mapuj(dane,mapowanie)
	wynik<-univariate_anal_stats(dyskretne, dane$def, czas,...)
	
	#funkcja univariate_anal_stats zmienia kolejnoœæ, dlatego wracamy j¹
	kolejnosc<-match(mapowanie$label,wynik$dyskretyzacja$label)
	#dodajê totala
	kolejnosc<-c(kolejnosc,length(kolejnosc)+1)
	wynik$dyskretyzacja<-wynik$dyskretyzacja[kolejnosc,]
	wynik$rozklady$obs_all_tbl<-wynik$rozklady$obs_all_tbl[kolejnosc,]
	wynik$rozklady$pct_all_tbl<-wynik$rozklady$pct_all_tbl[kolejnosc,]
	wynik$rozklady$avg_t_tbl<-wynik$rozklady$avg_t_tbl[kolejnosc,]
	
	wynik$dyskretyzacja<-cbind(wynik$dyskretyzacja[,1:3], 
			mapping_war=c(as.character(mapowanie$war),""),
			wynik$dyskretyzacja[,4:length(wynik$dyskretyzacja)])
	wynik$dyskretyzacja$discret<-""
	wynik$dyskretyzacja$discret[nrow(wynik$dyskretyzacja)]="<TOTAL>"
	wynik
}

#TODO! Coœ nie dzia³a obs³uga,  gdy w x s¹ NA
#np. zmienna Wiek
przypisz2<-function(x, bucket, interpol=FALSE, fitted=NULL, NA_substit=-2147483647)
{
	
	if (!is.null(fitted))
		bucket$fitted<-as.vector(fitted);
	
	if (is.null(bucket$fitted))
		stop("Brak kolumny 'bucket$fitted'.")
	
	if (is.factor(bucket$fitted))
		warning("przypisz2: Uwaga! bucket$fitted jest typu factor, co prowadzi do dziwnych wyników!");
	
	#inicjujê wektor z wynikami
	if (is.numeric(bucket$fitted))
		wynik<-rep(NA, length(x))
	else
		wynik<-rep("<NA>", length(x));
	#print(bucket)
	#Jeœli buckety okreœlone s¹ przez warunki mapuj¹ce
	jest_mapowanie=FALSE
	if (!is.null(bucket$mapping_war))
		if(any(!is.na(bucket$mapping_war)))
			jest_mapowanie=TRUE
	if (jest_mapowanie){
		mmm<-data.frame(war=bucket$mapping_war, label=bucket$fitted, discret=bucket$discret)
		mmm<-mmm[mmm$discret!="<TOTAL>",]
		wynik<-mapuj(x, mmm[,1:2])
	}
	else{
		#rozdzielam na wartoœci dyskretne i ci¹g³e
		ciagle_buck<-bucket[!is.na(bucket$od),];
		#wybieram dyskretne wartoœci i usuwam totala
		dyskretne_buck<-bucket[is.na(bucket$od) & bucket$discret!="<TOTAL>",];
		
		#czy s¹ jakieœ przedzia³y
		if (nrow(ciagle_buck)>0){
			if (interpol)
			{
				#obcinam zakresy do wartosci srodkow krancowych przedzialow
				x[x<min(ciagle_buck$srodek)]<-min(ciagle_buck$srodek);
				x[x>max(ciagle_buck$srodek)]<-max(ciagle_buck$srodek);
				
				przedzial<-findInterval(x, ciagle_buck$srodek, rightmost.closed = TRUE, all.inside = TRUE);
				
				od<-ciagle_buck$srodek[przedzial];
				do<-ciagle_buck$srodek[przedzial+1];
				wynik<-(x-od)/(do-od)*bucket$fitted[przedzial+1]+(do-x)/(do-od)*ciagle_buck$fitted[przedzial];
			}
			else
			{
				granice<-sort(unique(c(ciagle_buck$od, ciagle_buck$do)));
				if (length(granice)==1)
					granice<-c(granice, granice);
				przedzial<-findInterval(x, granice, rightmost.closed = TRUE, all.inside = TRUE);
				wynik<-ciagle_buck$fitted[przedzial];
			}
		}
		

		# Jeœli bucket zdefiniowa³ obs³ugê braków danych, to wszystkie braki podmieniam
		# wartoœci¹ specjaln¹ oznaczajac¹ brak danych
		if (any(rownames(bucket)== NA_substit))
			x[is.na(x)]<-NA_substit;
		
		## I jeszcze wartoœci dyskretne lub specjalne
		# gdzie s¹
		spec_bool<- x %in% dyskretne_buck$discret;
		#jakie_sa
		spec_idx<- match(x, dyskretne_buck$discret);
		#nadpisujê
		wynik[spec_bool]<- dyskretne_buck$fitted[na.omit(spec_idx)];
		
	}
	
	
	#sprawdzam, czy s¹ jakieœ nieprzypisane wartoœci. Jeœli tak, to rzucam ostrze¿enie;
	if (is.numeric(bucket$fitted)){
		if (any(is.na(wynik)))
			warning("przypisz2: Nie wszystkie wartoœci zosta³y przypisane do bucketa. Pozosta³y braki danych.")	
	} 	else if(any(wynik=="<NA>"))
		warning("przypisz2: Nie wszystkie wartoœci zosta³y przypisane do bucketa. Pozosta³y braki danych,
						oznaczone jako <NA>'.")
	
	if (!is.numeric(wynik))
		wynik<-factor(wynik, levels=unique(bucket$fitted[bucket$discret!="<TOTAL>"]));
	
	wynik
}





#'   Dyskretyzuje zmienn¹ ciag³¹ drzewkiem w oparciu o zmienn¹ odpowiedzi
#'
#' W parametrze \code{x} ani \code{y} nie mo¿e byæ NULLi. Mo¿na je jakoœ zakodowaæ.
#' @item x zmienna ci¹g³a.
#' @item y zmienna celu (PD, LGD).
#' @item special_val wartoœci specjalne, które powinny zostaæ uwzglêdnione jako
#'                   osobne klasy. Wartoœæ tak¹ jest równie¿ \code{NA}, automatycznie
#'                   uwzglêdniana jako osobna klasa.
#' @item interactive TRUE, jeœli zmienna ma byæ dyskretyzowana interaktywnie. W
#'                   przeciwnym razie, co jest wartoœci¹ domyœln¹, dyskretyzacja
#'                   jest automatyczna.
#' @item from zamiast automatycznego dzielenia, mo¿na podaæ wartoœci przedzia³ów (from,to].
#' @item to   zamiast automatycznego dzielenia, mo¿na podaæ wartoœci przedzia³ów (from,to].
#' @item ... inne parametry do funkcji \link{\code{drzewo}}.
#' @seealso \code{drzewo}
#' @author Micha³ Danaj
#TODO wyci¹gn¹æ parametry drzewa
numeric_var_treatment<-function(x, y, 
		special_val=numeric_var_treatment.params$spcial_val,
		min_bucket=numeric_var_treatment.params$min_bucket, 
		max_gleb=numeric_var_treatment.params$max_gleb,
		interactive=FALSE, locfit=FALSE, breaks=NULL,span=0.9, ...){
	if (length(x)!=length(y))
		stop("discretization: parametry 'x' i 'y' maj¹ ró¿ne d³ugoœci!");
	
	#Mimo, ¿e przygotowywya³em funkcjê do obs³ugi null-i, to rezygnujê z tego
	#ze wzglêdów bezpieczeñstwa.
	if (any(is.na(x)) | any(is.na(y)))
		stop ("discretization: W 'x' ani 'y' nie mo¿e byæ NA!");
	
	# Wartoœci specjalne
	special_idx<-is.na(x)|x %in% special_val;
	special_val<-unique(x[special_idx & !is.na(x)]);
	#czy s¹ NA
	sa_na<-any(is.na(x));
	
	## Obs³ugujê wartoœci ci¹g³e
	#jeœli zosta³y podane zakresy przedzi³ów, to dzielimy wg nich
	if (!is.null(breaks)){
		bucket_drzewo<-buckety_stat2(breaks, x[!special_idx], y[!special_idx], total=FALSE);
	} else {
		if (locfit){
			bucket_drzewo<-try(reg_nieparam(x[!special_idx],y[!special_idx], span=span, wytnij=0.01), TRUE)
			
			#jeœli wyliczy³o siê z b³êdem, to próbujê jeszcze na dwa sposoby...
			if (class(bucket_drzewo)=="try-error"){
				bucket_drzewo<-try(reg_nieparam(x[!special_idx],y[!special_idx], span=span, wytnij=0.01, buckets=50), TRUE)
				if (class(bucket_drzewo)=="try-error"){
					bucket_drzewo<-try(reg_nieparam(x[!special_idx],y[!special_idx], span=span, wytnij=0.01, buckets=30), TRUE)
					
					#jeœli wci¹¿ siê nie powiod³o, to zwracamy wartoœæ b³êdu
					if (class(bucket_drzewo)=="try-error")
						return(bucket_drzewo)
				}	
				print("  2  ")
			}
			
			#zmieniam nazwê, bo póŸniej j¹ wykorzystujê (historycznie)			
			bucket_drzewo$predicted=bucket_drzewo$fitted
		}
		else if (!interactive){
			bucket_drzewo<-drzewo(x[!special_idx],y[!special_idx], min_bucket=min_bucket, max_gleb=max_gleb, n_buckets=5, wytnij=0, ...)
			bucket_drzewo$predicted<-bucket_drzewo$br
		}else{
			bucket_drzewo<-interactive_tree(score=x[!special_idx],def=y[!special_idx],
					span=span, min_split=200, min_bucket=min_bucket,
					buckets=60, max_gleb=2)
			bucket_drzewo$predicted<-bucket_drzewo$br
		}
	}
	
	#przypisujê nr przedzia³ów dla wartoœci ci¹g³ych
	bucket_drzewo$fitted<-bucket_drzewo$nr;
	classing<-rep(NA,length(x));
	classing[!special_idx]<-przypisz(x[!special_idx], bucket_drzewo);
	
	#nadajê indeksy wartoœciom specjalnym i je przypisujê
	special_map<- -(length(special_val):1);
	names(special_map)<-special_val;
	classing[special_idx]<-special_map[as.character(x[special_idx])];
	
	#i jeszcze NA
	classing[is.na(x)]<- 0;
	
	#liczê statystyki
	classing_stat<-buckety_stat(classing, y, total=TRUE);
	
	#zmieniam nazwy wierszy, ¿eby nie by³y numery a labele klas
	#mapping<-c(names(special_map), "<NA>", rownames(bucket_drzewo), 'TOTAL');
	mapping<-c(names(special_map), "<NA>", rownames(bucket_drzewo), 'TOTAL');
	names(mapping)<-c(special_map,"0",bucket_drzewo$nr, 'TOTAL');
	
	rownames(classing_stat)<-mapping[rownames(classing_stat)];
	
	#które wartoœci s¹ specjalne (dyskretne)
	classing_stat$discret<-rep("",nrow(classing_stat));
	if (sa_na)
		classing_stat[c("<NA>", special_val),"discret"]<-rownames(classing_stat[c("<NA>", special_val),])
	else if (length(special_val)>0)
		classing_stat[as.character(special_val),"discret"]<-rownames(classing_stat[as.character(special_val),]);
	
	classing_stat$discret[nrow(classing_stat)]<-"<TOTAL>";
	
	classing_stat$do<-classing_stat$srodek<-classing_stat$od<-NA;
	classing_stat[classing_stat$discret=="",c("od","srodek","do")]<-bucket_drzewo[,c("od","srodek","do")];
	
	#dodajê predykcjê. W przypadku drzewka by³ to br, w przypadku locfit by³ to wynik regressji
	classing_stat$predicted<-NA
	classing_stat[classing_stat$discret=="",]$predicted<-bucket_drzewo$predicted
	#dodajê predykcjê do wartoœci specjalnych. Bêdzie to br
	classing_stat$predicted[classing_stat$discret!="" & classing_stat$discret!="<TOTAL>"]<-classing_stat$br[classing_stat$discret!="" & classing_stat$discret!="<TOTAL>"]
	
	classing_stat<-classing_stat[,c('nr','label','discret', 'od','srodek','do','n_good','pct_good','n_bad','pct_bad','n_obs','pct_obs',
					'br','predicted','woe','logit')];
	return(classing_stat);
}







#' Dyskretyzuje zmienn¹ i wylicza na niej statystyki
#'
#' W przypadku, gdy liczba unikalnych wartoœci zmiennej jest <= \code{discret_treshold}
#' lub zmienna nie jest zmienn¹ numeryczn¹,
#' uznaje ¿e zmienna jest dyskretna i jedynie wylicza dla niej statystyki. W przeciwnym
#' wypadku dyskretyzuje zmienn¹ i wylicza statystyki.
#' @item discret_treshold jeœli liczba unikalnych wartoœci zmiennej jest nie wiêksza
#'        ta wartoœæ, zmienna uznana jest za dyskretn¹ i nie jest poddawana dyskretyzacji.
#' @item interactive TRUE, jeœli zmienna ma byæ dyskretyzowana interaktywnie. W
#'                   przeciwnym razie, co jest wartoœci¹ domyœln¹, dyskretyzacja
#'                   jest automatyczna.
#' @item breaks zamiast automatycznego dzielenia, mo¿na podaæ wartoœci przedzia³ów (from,to].
#' @item forceContinous wymusza potraktowanie zmiennej jako ci¹g³¹, mimo ¿e liczba
#'                      unikalnych wartoœci jest mniejsza ni¿ \code{discret_treshold}.
#' @seealso \link{\code{buckety_stat}}.

univariate_anal_stats1b<-function(x,y, 
		locfit=FALSE, 
		discret_treshold=15,
		special_val=numeric_var_treatment.params$spcial_val, 
		max_gleb=3, 
		min_bucket=200, 
		interactive=FALSE,
		breaks=NULL, 
		mapping=NULL, 
		forceContinous=FALSE,
		span=0.9,...){
	
	if (length(x)!=length(y))
		stop("paramet ry 'x' i 'y' maj¹ ró¿ne d³ugoœci!");
	
	#Mimo, ¿e przygotowywya³em funkcjê do obs³ugi null-i, to rezygnujê z tego
	#ze wzglêdów bezpieczeñstwa.
	if (any(is.na(y)))
		stop ("W 'y' nie mo¿e byæ NA!");
	
	
	## jeœli s¹ jakieœ nulle w x, to odpowiednio siê nimi zajmujê
	nulle<-is.na(x)
	if (any(nulle)){
		
		#jeœli nulli jest mniej ni¿ za³o¿ona czêœæ populacji, to imputujê je. W przeciwnym razie przypisujê jako osobn¹ grupê.
		ile_nulli<-prop.table(table(nulle))
		if (ile_nulli["TRUE"]<numeric_var_treatment.params$nulle_do_imp_thr)
			#TODO zobaczyæ, czy y ma dwie wartoœci i jest to 0 i 1
			x<-missing_bin_target(x, y)
		else
			x[nulle]<-numeric_var_treatment.params$NA_substit;
		
	}
	
	## patrzê, czy nie ma skupisk w jakichœ wartoœciach. Jeœli tak, to bêdê je traktowaæ jako wartoœci specjalne
	## skupiska_freq<-prop.table(table(x))
	## skupiska <- skupiska_freq>numeric_var_treatment.params$separate_value_thr;
	## special_val<-unique(c(special_val,names(skupiska_freq)[skupiska]))
	
	
	## jeœli jest to zmienna dyskretna lub mapowanie
	if (!is.null(mapping)||((length(unique(x))<=discret_treshold || !is.numeric(x))&&
			is.null(breaks) && !forceContinous)){
		
		if (!is.null(mapping))
			x<-mapuj(x, mapping)	
		
		discret<-buckety_stat(x, y, total=TRUE);

		
		## uzupe³niam statystyki ##
		
		# ci¹g³e
		discret$od<-NA;
		discret$do<-NA;
		discret$srodek<-NA;
		
		#dyskretne
		nam<-rownames(discret)
		if (is.numeric(x)){
			nam[length(nam)]<-NA
			discret$discret<-as.numeric(nam)
		}
		else{
			nam[length(nam)]<-"<TOTAL>";
			discret$discret<-nam;
		}
	
		#mapowanie
		if(!is.null(mapping)) 
			discret$mapping_war<-mapping$war
		else
			discret$mapping_war<-NA
		
		discret<-discret[,c('nr','label','discret', 'mapping_war', 'od','srodek','do','n_good','pct_good','n_bad','pct_bad','n_obs','pct_obs',
						'br','woe','logit')]
		discret$predicted<-discret$br
	}
	## jeœli jest to zmienna ci¹g³a
	else{
		discret<-numeric_var_treatment(x,y, special_val=special_val,
				max_gleb=max_gleb,min_bucket=min_bucket,breaks=breaks,
				interactive=interactive, locfit=locfit, span=span, ...);
	}
	
	discret$label<-rownames(discret);
	return(discret);
}



model_dev<-function(predicted, target){
	-2*sum(target*log(predicted)+(1-target)*log(1-predicted))
}



#przypisuje woe na podstawie br z bucketa. Korzysta z funkcji \code{woe}, która 
#w razie niewyst¹pienia w buckecie wartoœci z klasy 0 lub 1, przyjmuje wartoœæ 0.5
#
#bucket_list - lista z opisami dyskretyzacji. Nazwy elementów listy powinny byæ zgodne z nazwami zmiennych w \code{data}
#data - \code{data.frame} z oryginalnymi zmiennymi
#var - wektor z nazwami zmiennych do ograniczenia


przypisz_woe<-function(bucket_list, data, vars=names(bucket_list), varname_sufix='woe'){
	
	data_out<-NULL
	
	for (zmienna in names(bucket_list)){
		
		if (!(zmienna %in% vars)) 
			next;
		
		####   wyliczam   woe    ######
		
		#Wyci¹gam element listy
		bucket<-bucket_list[[zmienna]]
		
		#jeœli badów lub goodów jest 0, to przyjmujê ¿e jest 0.5	
		
		pct_good = (pmax(bucket$n_good,0.5))/(max(bucket['TOTAL','n_good'],0.5))
		pct_bad = (pmax(bucket$n_bad,0.5))/(max(bucket['TOTAL','n_bad'],0.5))
		woe = log(pct_good/pct_bad)
		
		
		####   przypisujê woe    ######
		
		woe<-przypisz2(data[,zmienna],
				bucket_list[[zmienna]], 
				fitted=woe,
				NA_subst = numeric_var_treatment.params$NA_substit,
				interpol=FALSE)
		
		if (is.null(data_out))	{
			data_out <- data.frame(woe)
			names(data_out)<-zmienna
		}
		else
			data_out[,zmienna] <- c(woe)
		
	}
	
	names(data_out)<-paste(names(data_out), varname_sufix, sep="_")
	data_out
}	

#do skasowania!!!
woe<-function(br){
	br[br==0]<-0.5/length(br)
	br[br==1]<-(length(br)-0.5)/length(br)
	log(br/(1-br))
}





#generuje kod do zmiany roli zmiennych
#gen_code - jeœli TRUE, wynikiem funkcji jest kod do zmiany wartoœci. W przeciwnym razie, zwracany
#			jest data.frame ze zmienionymi rolami
editVariablesRole<-function(zmienne_rola, pattern=NULL, gen_code=TRUE){
	
	rola<-c("rejected", "explanatory", "target", "keep");
	
	if (is.null(pattern))
		wynik<-edit(zmienne_rola)
	else{
		indeksy<-grep(pattern, rownames(zmienne_rola), ignore.case = TRUE)
		wynik<-zmienne_rola
		wynik[indeksy,]<-edit(zmienne_rola[indeksy,])
		
		if (any(!(wynik$rola%in%rola))){
			blad<-"Niepoprawne wartoœci ról!"
			#stop(rownames(zmienne_rola)[!wynik$rola%in%rola])
			stop(blad)
		}
	}
	
	#sprawdzam, których zmiennych role siê zmieni³y
	zmiany<- zmienne_rola$rola!=wynik$rola
	nazwy_zmienionych<-rownames(zmienne_rola)[zmiany]
	nowe_wartosci<-wynik$rola[zmiany]
	
	#dla tych zmiennych generujê kod do zmiany ich ról
	nazwy_przecinek<-paste("c('", paste(nazwy_zmienionych, collapse = "','"), "')", sep="")
	wartosci_przecinek<-paste("c('", paste(nowe_wartosci, collapse = "','"), "')", sep="")
	
	if (gen_code==TRUE)
		paste(deparse(substitute(zmienne_rola)),"[",nazwy_przecinek,",'rola']<-", wartosci_przecinek, sep="")
	else
		wynik
}



################################################################################
########		TEstowe funckje z wagami   #####################################
################################################################################

## ## x<-rnorm(10000)
## ## y<-as.numeric(x+2*runif(1000)>0.5)
## ## 
## ## waga<-rep(1,length(x))
## ## waga[y==1]<-0.1
## ## #weights=waga
## ## 
## ## 
## ## wyn1<-reg_nieparam(x,y,30)
## ## windows(10,10)
## ## wyn2<-wtd.reg_nieparam(x,y,10, weights=waga)
## ## 
## ## score=x
## ## default=y
## ## 

wtd.reg_nieparam<-function (score, default, buckets = 100, subset=NULL, wytnij = 0, span = 0.7,
		degree = 2, plot = TRUE, target = "br", new = TRUE, col_points = "black",
		col_line = "darkblue", index = FALSE, weights=NULL, estim=NULL, ...)
{
	dane <- data.frame(score, default)
	
	if (is.null(weights)){
		weights=rep(1,length(score))
		weighst_locfit=weights	
	}
	else{
		
		# TODO zmieniæ to!
		# na cele tego zadania przyhardkorzy³em i obchodzê b³¹d locfita w ten sposób,
		# ¿e ustalam wagê 1 dla obserwacji z target=1
		
		waga_1<-weights[which(default==1)[1]]
		weights_locfit<-weights/waga_1
	}
	
	#jeœli okreœlono subset, to ograniczam dane na których pracujemy 
	if (!is.null(subset)){
		dane<-dane[subset,]
		weights<-weigths[subset]
		estim<-estim[subset]
	}
	
	if (wytnij > 0){
		do_usuniecia<-usun_konce(dane$score, prob = wytnij, weights=weights);
		if (length(do_usuniecia)>0)
			dane <- dane[-do_usuniecia,]
	}
	bucket <- buckety_br(x=dane$score, y=dane$default, n=buckets, method = "eq_count", weights=weights)
	

	
	#jeœli s¹ dwie wartoœci y, to uznaje to jest to zmienna binarna, i stosuje regresjê logistyczn¹
	if (length(unique(default)) == 2)
		l <- locfit(default ~ lp(score, nn = span), family = "binomial",
				link = "logit", data = dane, weights=weights_locfit)
	#w przeciwnym razie robi regresje liniow¹
	else 
		l <- locfit(default ~ lp(score, nn = span), data = dane, weights=weights)
	
	b2 <- predict(l, newdata = bucket$srodek)
	if (target == "br")
		bucket2 <- cbind(bucket, fitted = b2)
	else bucket2 <- cbind(bucket, fitted = log(b2/(1 - b2)))
	
	#liczê wielkoœæ b¹belka
	skala <- sqrt(bucket$n_obs/(sum(weights)/buckets))
	
	#liczê wartoœci wyestymowane
	#estim_aggr<-buckety_stat(b2, default, )
	
	#rysowanie
	x <- bucket2$srodek
	if (index)
		x <- bucket$nr
	if (plot) {
		if (new == TRUE)
			plot(x, with(bucket2, get(target)), col = col_points,
					cex = skala, ...)
		else points(x, with(bucket2, get(target)), cex = skala,
					col = col_points, ...)
		lines(x, bucket2$fitted, col = col_line, ...)
	}
	bucket2
}


usun_konce<-function (score, prob = 0.01, weights=NULL)
{
	
	if (is.null(weights))
		weights<-rep(1,length(score))
	
	po_score<-tapply(weights, score, sum)
	s <- cumsum(po_score/sum(weights))
	
	new_min <- as.numeric(names(which.max(s[s <= prob])))
	
	if (length(new_min) == 0)
		new_min <- -.Machine$double.xmax
	
	#od którego elementu powinienm wycinaæ wartoœci przekraczaj¹ce 1-prob
	temp<-which(s >= 1 - prob)[1]+1;
	if (is.na(temp))
		new_max <-  .Machine$double.xmax
	else
		new_max <- as.numeric(names(s[temp]));
	
	return(which(score <= new_min | new_max <= score))
}

## buckety_br(x, y, 10, method = "eq_count", weights=waga)
## wtd.quantile(x, prob=0:n/n, type='quantile', weights=waga)

buckety_br<-function(x, y, n, method=c("eq_length", "eq_count"), one.value.action=c("none","combine"),
		weights=NULL)
{                                         
#	TODO - domyœlna wartoœæ method
	method<-match.arg(method);
	
	if (is.null(weights))
		weights<-rep(1,length(x))
	
	#dolny kraniec bucketa (wartosc x)
	od=c(1:n);
	#gorny kranieb bucketa (wartosc x)
	do=c(1:n);
	
	if (method=="eq_length")
		granice<-seq(min(x),max(x),length.out=n+1)
	else{
		#w przypadku, gdy jedna wartoœæ jest dla wielu kwantyli, zdarzaj¹ siê problemy numeryczne
		#¿e wartoœæ teoretycznie jest taka sama, ale ró¿ni siê na 15-tym miejscu po przecinku.
		#tak na szybko, brute force obejœcie: sortujê
		granice<-sort(as.vector(unique(wtd.quantile(x, prob=0:n/n, type='quantile', weights=weights))));
	}
	#oznaczam, do ktorego przedzialu nalezy dana obserwacja
	przedzial<-findInterval(x, granice, rightmost.closed = TRUE, all.inside = TRUE);
	
	#i licze potrzebne rzeczy
	od<-granice[1:(length(granice)-1)];
	do<-granice[2:length(granice)];		
	
	srodek<-as.vector(tapply(x, przedzial, median));

	n_obs<- as.vector(tapply(weights, przedzial, sum));
	mean<-as.vector(tapply(x*weights, przedzial, FUN=sum))/n_obs;
	n_default<- as.vector(tapply(y*weights, przedzial, sum));
	br<- n_default/n_obs;
	logit<-log(br/(1-br));
	probit<-qnorm(br);
	
	#z tego co pamiêtam, mamy dziêki temu wyrzuciæ puste przedzia³y
	zostaw<-sort(unique(przedzial));
	
	as.data.frame(cbind(nr=zostaw,  od=od[zostaw], do=do[zostaw], 
					srodek, mean,  n_default, n_obs, br, logit, probit, var=br*(1-br)/n_obs))
}


to_excel<-function(table, row.names=FALSE){
	write.table(table, file='clipboard', sep='\t', row.names=row.names, dec=',')
}