library("qtl")
library("ASMap")
library("xlsx")
Joppa10Map = read.cross("csvr","Y:/Yuan/integrated map/Joppa10Ae564-corrected maps after publication",
           "GBS_90K_Joppa_YL_2019.csv",genotypes=c("A","H","B"),na.strings = "-"
           ,estimate.map = F, F.gen = 7,crosstype = "riself")
identicalMarker = findDupMarkers(Joppa10Map,exact.only = F,adjacent.only = F)

duplicateMarker = unlist(identicalMarker)

Joppa10MapNoDup = drop.markers(Joppa10Map,duplicateMarker)

Joppa10MapNoDup = drop.markers(Joppa10MapNoDup,
                               names(ntyped(Joppa10MapNoDup,"mar")[ntyped(Joppa10MapNoDup,"mar") < 100]))

summary(Joppa10MapNoDup)

#plotMissing(Joppa10MapNoDup)

mapConstruction = mstmap(Joppa10MapNoDup,id = "Genotype", bychr = TRUE,
                         suffix = "numeric", anchor = FALSE, dist.fun = "kosambi",
                         objective.fun = "COUNT", p.value = 1e-18, noMap.dist = 30,
                         noMap.size = 0, miss.thresh = 1, mvest.bc = FALSE,
                         detectBadData = FALSE, return.imputed = FALSE,
                         trace = FALSE)
mapConstruction$geno$UN.1$map

setwd("C:/Users/liu.yuan/Desktop/temp")


for(i in 1:length(mapConstruction$geno)){
  if(length(mapConstruction$geno[[i]]$map )> 15){
    write.xlsx(mapConstruction$geno[[i]]$map,"temp_map.xlsx",sheetName = names(mapConstruction$geno)[i],
               col.names = T,row.names = T,append = T,showNA = T)    
  }
}

##################################
##DP527
DP527GBS90K = read.cross("csvr","Y:/Yuan/integrated map/DP527",
                         "DP527_GBS_90K.csv",genotypes=c("A","H","B"),na.strings = "-"
                         ,estimate.map = F, F.gen = 7,crosstype = "riself")
identicalMarker = findDupMarkers(DP527GBS90K,exact.only = F,adjacent.only = F)
duplicateMarker = unlist(identicalMarker)
DP527GBS90KNoDup = drop.markers(DP527GBS90K,duplicateMarker)
DP527GBS90KNoDup = drop.markers(DP527GBS90KNoDup,
                                names(ntyped(DP527GBS90KNoDup,"mar")[ntyped(DP527GBS90KNoDup,"mar") < 110]))

DP527MapConstruction = mstmap(DP527GBS90KNoDup,id = "Genotype", bychr = TRUE,
                              suffix = "numeric", anchor = FALSE, dist.fun = "kosambi",
                              objective.fun = "COUNT", p.value = 1e-18, noMap.dist = 30,
                              noMap.size = 0, miss.thresh = 1, mvest.bc = FALSE,
                              detectBadData = FALSE, return.imputed = FALSE,
                              trace = FALSE)

setwd("C:/Users/liu.yuan/Desktop/temp")

for(i in 1:length(DP527MapConstruction$geno)){
  if(length(DP527MapConstruction$geno[[i]]$map )> 15){
    write.xlsx(DP527MapConstruction$geno[[i]]$map,"temp_map.xlsx",sheetName = names(DP527MapConstruction$geno)[i],
               col.names = T,row.names = T,append = T,showNA = T)    
  }
}

lapply(X = identicalMarker[1:10], FUN = function(x) {
  write.table(x, append = T, file = "test.csv", ncolumns = length(x),sep = ",")
})

for(i in 1:length(identicalMarker)){
  write(c(names(identicalMarker[i]),identicalMarker[[i]]), file = "test.csv"
        ,append = T, ncolumns = length(identicalMarker[i])+1,sep = ",")
}


#####################################
##BP025
BP025GBS9K = read.cross("csvr","Y:/Yuan/integrated map/Ben_PI41025",
                        "BP025_9K_GBS_test.csv",genotypes=c("A","H","B"),na.strings = "-"
                        ,estimate.map = F, F.gen = 7,crosstype = "riself")
identicalMarker = findDupMarkers(BP025GBS9K,exact.only = F,adjacent.only = F)
duplicateMarker = unlist(identicalMarker)
BP025GBS9KNoDup = drop.markers(BP025GBS9K,duplicateMarker)

BP025MapConstruction = mstmap(BP025GBS9KNoDup,id = "Genotype", bychr = TRUE,
                              suffix = "numeric", anchor = FALSE, dist.fun = "kosambi",
                              objective.fun = "COUNT", p.value = 1e-15, noMap.dist = 30,
                              noMap.size = 0, miss.thresh = 1, mvest.bc = FALSE,
                              detectBadData = FALSE, return.imputed = FALSE,
                              trace = FALSE)

setwd("C:/Users/liu.yuan/Desktop/temp")

for(i in 1:length(BP025MapConstruction$geno)){
  if(length(BP025MapConstruction$geno[[i]]$map )> 15){
    write.xlsx(BP025MapConstruction$geno[[i]]$map,"temp_map_BP025.xlsx",sheetName = names(BP025MapConstruction$geno)[i],
               col.names = T,row.names = T,append = T,showNA = T)    
  }
}
