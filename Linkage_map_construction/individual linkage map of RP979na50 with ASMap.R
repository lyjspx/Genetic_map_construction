#Set working directory under Mac
setwd("/Users/xuehui.li/Desktop/XL/research/crops/durum_wheat/p2_Integrated genetic linkage map in tetraploid wheat/linkage map/RP979")

#load packages "qtl"
RP979=read.cross("csvr", "/Users/xuehui.li/Desktop/XL/research/crops/durum_wheat/p2_Integrated genetic linkage map in tetraploid wheat/linkage map/RP979", "RP979_rotated_Rqtl.csv", na.strings=c("-"), estimate.map = FALSE, crosstype = "riself")
RP979=read.cross("csvr", "/Users/xuehui.li/Desktop/XL/research/crops/durum_wheat/p2_Integrated genetic linkage map in tetraploid wheat/linkage map/RP979", "RP979_rotated_Rqtl_6624mapped.csv", na.strings=c("-"), estimate.map = FALSE, crosstype = "riself")

str(RP979)
nind(RP979)
nchr(RP979)
totmar(RP979)
nmar(RP979)
plot(RP979)
plotMap(RP979)
plotMissing(RP979)

#remove samples with a lot of missing marker values
par(mfrow=c(1,1), las=1)
plot(ntyped(RP979), xlab="Sample", ylab="No. typed markers", main="No. typed markers by samples") 
RP979 <- subset(RP979, ind=(ntyped(RP979)>2800))
nind(RP979)

#remove the markers with lot of missing values
plot(ntyped(RP979, "mar"), xlab="Marker", ylab="No. typed samples",main="No. typed samples by markers")
nt.bymar <- ntyped(RP979, "mar")
todrop <- names(nt.bymar[nt.bymar < 91])
RP979na50 <- drop.markers(RP979, todrop)

nind(RP979na50)
totmar(RP979na50)
nmar(RP979na50)

#identify duplicate individuals
cg <- comparegeno(RP979na50)
str(cg)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes") 
rug(cg[lower.tri(cg)])
wh <- which(cg > 0.95, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
genotype=pull.geno(RP979na50)
str(genotype)
dim(genotype)
table(genotype[70,], genotype[81,])

gc=genClones(RP979na50, tol=0.9)
gc$cgd
RP979na50=fixClones(RP979na50, gc$cgd, consensus = TRUE)

#remove extremly distorted markers
##estimate segregation distortion with Bonferroni correction
gt <- geno.table(RP979na50)
str(gt)
dim(gt)
gt[1:5,]

todrop1 <- rownames(gt[gt$P.value < 1e-10,])
length(todrop1)
todrop2=rownames(gt[gt$P.value < 0.05/totmar(RP979na50),])
length(todrop2)

##estimate minor allele frequency
MAF=rep(NA, 6735)
for (i in 1:6735)
{
MAF[i]=min(gt$AA[i], gt$BB[i])/sum(gt$AA[i]+gt$BB[i])
}

MAF[1:5]
gt=data.frame(gt, MAF)
str(gt)
todrop.MAF=rownames(gt[gt$MAF<0.3, ])
length(todrop.MAF)

RP979na50 <- drop.markers(RP979na50, todrop.MAF)

nind(RP979na50)
totmar(RP979na50)
nmar(RP979na50)

sg1=statGen(RP979na50, stat.type = c("xo", "dxo", "miss"))
profileGen(RP979na50, stat.type = c("xo", "dxo", "miss"))
sm1=statMark(RP979na50,stat.type = c("marker", "interval"), map.function ="kosambi")
profileMark(RP979na50, stat.type = "marker", use.dist = TRUE,map.function = "kosambi")
profileMark(RP979na50, stat.type = "seg.dist", use.dist = TRUE,map.function = "kosambi")


#load packages "asmap" and construct linkage map
RP979na50.map.pE15=mstmap(RP979na50, id="Genotype", pop.type = "RIL7", dist.fun = "kosambi", suffix = "numeric", anchor = FALSE, objective.fun = "COUNT", p.value = 1e-15, noMap.dist = 15,
                          noMap.size = 0, miss.thresh = 1, mvest.bc = FALSE, detectBadData =FALSE, return.imputed = FALSE, trace = FALSE)

RP979na50.map.pE18=mstmap(RP979na50, id="Genotype", pop.type = "RIL7", dist.fun = "kosambi", suffix = "numeric", anchor = FALSE, objective.fun = "COUNT", p.value = 1e-18, noMap.dist = 15,
                          noMap.size = 0, miss.thresh = 1, mvest.bc = FALSE, detectBadData =FALSE, return.imputed = FALSE, trace = FALSE)

RP979na50.map.pE21=mstmap(RP979na50, id="Genotype", pop.type = "RIL7", dist.fun = "kosambi", suffix = "numeric", anchor = FALSE, objective.fun = "COUNT", p.value = 1e-21, noMap.dist = 15,
                     noMap.size = 0, miss.thresh = 1, mvest.bc = FALSE, detectBadData =FALSE, return.imputed = FALSE, trace = FALSE)

RP979na50.map.pE24=mstmap(RP979na50, id="Genotype", pop.type = "RIL7", dist.fun = "kosambi", suffix = "numeric", anchor = FALSE, objective.fun = "COUNT", p.value = 1e-24, noMap.dist = 15,
                     noMap.size = 0, miss.thresh = 1, mvest.bc = FALSE, detectBadData =FALSE, return.imputed = FALSE, trace = FALSE)

str(RP979na50.map)
nmar(RP979na50.map.pE15)
plotMap(RP979na50.map.pE15)
nmar(RP979na50.map.pE18)
plotMap(RP979na50.map.pE18)
nmar(RP979na50.map.pE21)
plotMap(RP979na50.map.pE21)
nmar(RP979na50.map.pE24)
plotMap(RP979na50.map.pE24)

pull.map(RP979na50.map.pE15)[[1]]
pull.map(RP979na50.map.pE15)[[2]]
pull.map(RP979na50.map.pE15)[[3]]
pull.map(RP979na50.map.pE15)[[4]]
pull.map(RP979na50.map.pE15)[[5]]
pull.map(RP979na50.map.pE15)[[6]]

pull.map(RP979na50.map.pE18)[[1]]
pull.map(RP979na50.map.pE18)[[2]]
pull.map(RP979na50.map.pE18)[[3]]
pull.map(RP979na50.map.pE18)[[4]]
pull.map(RP979na50.map.pE18)[[5]]
pull.map(RP979na50.map.pE18)[[6]]
pull.map(RP979na50.map.pE18)[[7]]
pull.map(RP979na50.map.pE18)[[8]]
pull.map(RP979na50.map.pE18)[[9]]
pull.map(RP979na50.map.pE18)[[10]]
pull.map(RP979na50.map.pE18)[[11]]

pull.map(RP979na50.map.pE21)[[1]]
pull.map(RP979na50.map.pE21)[[2]]
pull.map(RP979na50.map.pE21)[[3]]
pull.map(RP979na50.map.pE21)[[4]]
pull.map(RP979na50.map.pE21)[[5]]
pull.map(RP979na50.map.pE21)[[6]]
pull.map(RP979na50.map.pE21)[[7]]
pull.map(RP979na50.map.pE21)[[8]]
pull.map(RP979na50.map.pE21)[[9]]
pull.map(RP979na50.map.pE21)[[10]]
pull.map(RP979na50.map.pE21)[[11]]
pull.map(RP979na50.map.pE21)[[12]]
pull.map(RP979na50.map.pE21)[[13]]
pull.map(RP979na50.map.pE21)[[14]]

names(RP979na50.map.pE18$geno)=c("1A", "un1", "un2", "1B1", "5A", "3A_7B_7A_3B", "4A1", "4A_2B_2A", "5B_4B", "6A_6B", "1B2")
names(pull.map(RP979na50.map.pE18))

RP979na50.map.pE18.br=breakCross(RP979na50.map.pE18, split=list('3A_7B_7A_3B'=c("3A_3545059","7B_652390","7A_180999206"), '4A_2B_2A'=c("7A_159655817","2B_2692666"), '5B_4B'="5B_53893", '6A_6B'="6A_207552873"),
                                 suffix =list('3A_7B_7A_3B' = c("3A", "7B", "7A", "3B"), '4A_2B_2A'=c("4A2", "2B", "2A"), '5B_4B'=c("5B", "4B"), '6A_6B'=c("6A", "6B")))

names(pull.map(RP979na50.map.pE18.br))
nmar(RP979na50.map.pE18.br)

RP979na50.map.pE18.br_mer=mergeCross(RP979na50.map.pE18.br, merge = list('1B' = c("1B1", "1B2"), '4A'=c("4A1", "4A2")))

names(pull.map(RP979na50.map.pE18.br_mer))
nmar(RP979na50.map.pE18.br_mer)
pull.map(RP979na50.map.pE18.br_mer)[[7]]

RP979na50.map.pE18.bm=mstmap(RP979na50.map.pE18.br_mer, dist.fun = "kosambi", suffix = "numeric", bychr=TRUE, anchor = FALSE, trace = FALSE, p.value=1e-6)
nmar(RP979na50.map.pE18.bm)
plot(RP979na50.map.pE18.bm)
plotMap(RP979na50.map.pE18.bm)
str(pull.map(RP979na50.map.pE18.bm))

RP979na50.map.pE18.bm=breakCross(RP979na50.map.pE18.bm, split=list('4A'="4A_483620"), suffix =list('4A' = c("4A", "un3")))
RP979na50.map.pE18.bm=mstmap(RP979na50.map.pE18.bm, dist.fun = "kosambi", suffix = "numeric", bychr=TRUE, anchor = FALSE, trace = FALSE, p.value=1e-6)
nmar(RP979na50.map.pE18.bm)
plot(RP979na50.map.pE18.bm)
plotMap(RP979na50.map.pE18.bm)
str(pull.map(RP979na50.map.pE18.bm))
pull.map(RP979na50.map.pE18.bm)[[17]]

sg2=statGen(RP979na50.map.pE18.bm, stat.type = c("xo", "dxo", "miss"))
profileGen(RP979na50.map.pE18.bm, stat.type = c("xo", "dxo", "miss"))
profileGen(RP979na50.map.pE18.bm, stat.type = "xo")
profileGen(RP979na50.map.pE18.bm, stat.type = c("dxo"))
profileGen(RP979na50.map.pE18.bm, stat.type = c("miss"))
sm2=statMark(RP979na50.map.pE18.bm, stat.type = c("marker", "interval"))
profileMark(RP979na50.map.pE18.bm, stat.type = "marker", use.dist = TRUE)
profileMark(RP979na50.map.pE18.bm, stat.type = c("seg.dist"), use.dist = TRUE)


RP979na50.GeneticMap=c(pull.map(RP979na50.map.pE18.bm)[[1]], pull.map(RP979na50.map.pE18.bm)[[2]],pull.map(RP979na50.map.pE18.bm)[[3]], pull.map(RP979na50.map.pE18.bm)[[4]],
pull.map(RP979na50.map.pE18.bm)[[5]], pull.map(RP979na50.map.pE18.bm)[[6]],pull.map(RP979na50.map.pE18.bm)[[7]], pull.map(RP979na50.map.pE18.bm)[[8]],
pull.map(RP979na50.map.pE18.bm)[[9]], pull.map(RP979na50.map.pE18.bm)[[10]],pull.map(RP979na50.map.pE18.bm)[[11]], pull.map(RP979na50.map.pE18.bm)[[12]],
pull.map(RP979na50.map.pE18.bm)[[13]], pull.map(RP979na50.map.pE18.bm)[[14]])

write.table(RP979na50.GeneticMap, "RP979na50.GeneticMap_ASMap.txt")

#Synteny between genetic linkage map and physical map
RP979_GvsP=read.delim("RP979na50.GeneticMap.6625markers.txt", header=T, na.strings=c("NA"))
RP979_GvsP=read.delim("RP979na50.GeneticMap.6624markers.txt", header=T, na.strings=c("NA"))

dim(RP979_GvsP)
RP979_GvsP[1:5,]

LG=RP979_GvsP[,2]
chr=RP979_GvsP[,5]
order.genet=RP979_GvsP[,4]
order.phy=RP979_GvsP[,7]

#plot positions
xyplot(RP979_GvsP[,6]~RP979_GvsP[,3]|LG, data=RP979_GvsP, xlab="RP979 genetic linkage map", ylab="CS physical map")
xyplot(RP979_GvsP[,6]~RP979_GvsP[,3]|chr, data=RP979_GvsP, xlab="RP979 genetic linkage map", ylab="CS physical map")

#plot order
plot(RP979_GvsP[,4],RP979_GvsP[,7], type="p",pch=15, cex=0.5, xlim=c(1, 6625), ylim=c(1,6625), xlab="RP979 genetic linkage map", ylab="CS physical map", font.lab=2,xaxt='n', yaxt='n', xaxs = "i", yaxs = "i")
abline(v=473.5)
abline(v=640.5)
abline(v=1252.5)
abline(v=1797.5)
abline(v=2276.5)
abline(v=3006.5)
abline(v=3446.5)
abline(v=3807.5)
abline(v=4273.5)
abline(v=4856.5)
abline(v=5252.5)
abline(v=5713.5)
abline(v=6274.5)

abline(h=464.5)
abline(h=666.5)
abline(h=1270.5)
abline(h=1823.5)
abline(h=2252.5)
abline(h=3056.5)
abline(h=3473.5)
abline(h=3852.5)
abline(h=4264.5)
abline(h=4881.5)
abline(h=5297.5)
abline(h=5730.5)
abline(h=6283.5)

axis(1, c(236.5,556.5,946,1524.5,2036.5,2641,3226.5,3627.5, 4041,4565.5, 5055, 5483.5, 5994.5, 6450), labels=c("1A", "1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B","7A","7B"),las=1)
axis(2, c(232.5,565.5,968,1546.5,2037.5,2654,3264.5,3663, 4059,4573.5, 5090, 5514.5, 6007.5, 6454.5), labels=c("1A", "1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B","7A","7B"),las=2)

#heat map of constructed linkage map
plot.rf(RP979na50.map.pE18.bm)
heatMap(RP979na50.map.pE18.bm, lmax=20)
heatMap(RP979na50.map.pE218bm, chr, mark, what = c("both", "lod", "rf"), lmax = 12, rmin = 0, markDiagonal = FALSE, color =rev(colorRampPalette(brewer.pal(11,"Spectral"))(256)))


#QTL mapping
RP979na50.RustQ <- calc.genoprob(RP979na50, step=1, error.prob=0.01)
RP979na50.RustQ <- calc.genoprob(RP979na50.map.pE18, step=1, error.prob=0.01)
RP979na50.RustQ <- calc.genoprob(RP979na50.map.pE18.bm, step=1, error.prob=0.01)

# single-QTL genome scan with "scanone"
out.hk.TTKSK <- scanone(RP979na50.RustQ, method="hk", pheno.col=2)
out.hk.TRTTF <- scanone(RP979na50.RustQ, method="hk", pheno.col=3)
out.hk.TMLKC <- scanone(RP979na50.RustQ, method="hk", pheno.col=4)

summary(out.hk.TTKSK)
summary(out.hk.TRTTF)
summary(out.hk.TMLKC)

plot(out.hk.TTKSK)
plot(out.hk.TRTTF)
plot(out.hk.TMLKC)

operm.hk <- scanone(RP979na50.RustQ, method="hk", n.perm=1000)
summary(operm.hk, alpha=0.05)

