### Libraries ###

library(MALDIquant)
library(MALDIquantForeign)
library(dplyr)
library(purrr)
library(tidyr)
library(binda)
library(broom)

### Custom Functions ###

Binary.analisis<-function(feature.matrix,Peaks,spot,
                          sensi,other,title) {
    Ytrain1<-sensi
    categorias<-table(Ytrain1)
    total.de.NA<-sum(is.na(intensityMatrix(Peaks))) 
    promedio.de.NA<-mean(apply(is.na(intensityMatrix(Peaks)), 1, sum)) 
    mz.full <- as.double(colnames(feature.matrix))
    mz <- round( mz.full )
    picos.ducplicados<-any(duplicated(mz))
    thr <- optimizeThreshold(feature.matrix, Ytrain1)
    Xtrain.b <- dichotomize(feature.matrix, thr)
    colnames(Xtrain.b) <-mz
    rownames(Xtrain.b) <- paste(spot,other,sensi, sep=".")
    Chequeo.dicotomico<-is.binaryMatrix(Xtrain.b)
    Xtrain.b.naive <- ifelse(is.na(intensityMatrix(Peaks)), as.integer(0), as.integer(1))
    rownames(Xtrain.b.naive) <-paste(spot,other,sensi, sep=".")
    colnames(Xtrain.b.naive) <-mz
    Chequeo.binario<-is.binaryMatrix(Xtrain.b.naive)
    Feature.dicho3<-Xtrain.b
    Feature.bina3<-Xtrain.b.naive
    save(Feature.dicho3, Feature.bina3,
         file=paste("Matrix", title, "rda", sep="."))
    print(categorias)
    print(paste("Total NA",total.de.NA,sep=" "))
    print(paste("Average NA",promedio.de.NA,sep=" "))
    print(paste("Any duplicate peak",picos.ducplicados,sep=" "))
    print(paste("Dicho is binary",Chequeo.dicotomico,sep=" "))
    print(paste("Bina is binary",Chequeo.binario,sep=" "))
}

Peaks_extraction_BDA<-function(Matrix, Clasif, Num.peaks, 
                               Name.Class1,Name.Class2, title) {
    br <- binda.ranking(Matrix, Clasif, verbose=FALSE)
    br2 <- br[1:Num.peaks, c(1:4)]
    selPeaks23<-br2[,2]
    Most.dif2<-tidy(selPeaks23)
    BDA.Rank<-Most.dif2
    names(BDA.Rank)<-c("Picos","BDA.score") 
    selPeaks232<-br2[,3]
    Most.dif22<-tidy(selPeaks232)
    names(Most.dif22)<-c("Picos",Name.Class1)
    selPeaks233<-br2[,4]
    Most.dif23<-tidy(selPeaks233)
    names(Most.dif23)<-c("Picos",Name.Class2)
        BDA.Rank.2<-BDA.Rank%>%
            left_join(Most.dif22,by="Picos")%>%
            left_join(Most.dif23,by="Picos")
        BDA.Rank<-gather(BDA.Rank.2,Class,t.score,c(3:4))
    Picos.significativos<-BDA.Rank
    save(Picos.significativos,
         file=paste("Peak", title, "rda", sep="."))
}


### Read Bruker data ###

so2<-importBrukerFlex("Path_to_spectra_folder", 
                      verbose=FALSE)

# Labeling spectra
Newman.table<-read.csv2("Path_to_csv_data")
Newman.table<-Newman.table%>%
    mutate(spot.a.1=paste(ID,Bacteria,Medio.Cultivo,sep="."))

names<-Newman.table$spot.a.1
repeticion<-Newman.table$n.esp
Names.Sa.S.Pocillo<- rep(names,repeticion)
Names.arr<-array(Names.Sa.S.Pocillo)
for (h in 1:length(so2)) {
    metaData(so2[[h]])$spot<-Names.arr[[h]]
}
# Spectra Pre-proccessing
so5<-so2
spot.factor <- factor(sapply(so5,function(x)metaData(x)$spot))
so5 <- trim(so5)
so5 <- transformIntensity(so5, method="sqrt")
so5 <- smoothIntensity(so5, method="SavitzkyGolay",
                       halfWindowSize=30)
so5 <- removeBaseline(so5, method="SNIP", iterations=100)
so5 <- calibrateIntensity(so5, method="TIC")
so5 <- alignSpectra(so5, halfWindowSize=30,SNR=5,
                    tolerance=0.2, warpingMethod="quadratic")
avgSpectra.Newman.HCCA <-averageMassSpectra(so5, labels=spot.factor, method="sum")
peaks <- detectPeaks(avgSpectra.Newman.HCCA, SNR=5, method="MAD", halfWindowSize=30)
peaks <- binPeaks(peaks,tolerance=0.2)
species.Ave<-factor(Newman.table$Growth.medium) 
spot.factor.Avera<-factor(Newman.table$spot.a.1) 
peaks <- filterPeaks(peaks, minFrequency=c(rep(1/2,2)),
                     labels = species.Ave,
                     mergeWhitelists=TRUE)
featureMatrix <- intensityMatrix(peaks, avgSpectra.Newman.HCCA)

# Binarization

Binary.analisis(feature.matrix= featureMatrix,
                Peaks=peaks,
                spot= spot.factor.Avera,
                sensi=species.Ave,
                other="",
                title="NewmanData")
load("Matrix.NewmanData.rda")

# Discriminant-peaks selection

Peaks_extraction_BDA(Feature.dicho3, species.Ave, 30 ,"Sal","Sin sal","30.BDA")
load("Peak.30.BDA.rda")





