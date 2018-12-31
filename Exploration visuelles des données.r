
library("IRdisplay")
display_html('<style> h1,h2,h3,h4,h5,h6 {text-align:center}</style>')

libraryPlus=function(p){ if(!require(p,character.only = T)){install.packages(p,character.only = T)}; library(p,character.only = T)}
librariesPlus=function(...){for(i in list(...)){libraryPlus(i)}}
librariesPlus(
    "FactoMineR",
    "factoextra",
    "dimRed",
    "lle",
    "Rtsne",
    "igraph",
    "TSclust",
    "vegan",
    "ade4",
    "MASS",
    "ggplot2",
    "gridExtra",
    "dplyr",
    "tidyverse",
    'caret',
    "virtualspecies",
    "devtools",
    "mctest",
    "memisc",
    "kohonen",
    "ggrepel",
    "mclust",
    "corrplot")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if(!require("RDRToolbox")) BiocManager::install("RDRToolbox", version = "3.8")
library("RDRToolbox")
#install.packages("devtools")
#options(unzip = "internal")
#install_github("Displayr/flipMultivariates")library(mctest)



rawdata = read.csv2("pomeroy-2002-v2_database.txt",sep="\t",header = T,row.names = 1)
data = read.csv2("pomeroy-2002-v2_database.txt",sep="\t",header = F,row.names = 1,skip = 2)
classes =read.csv2("pomeroy-2002-v2_database.txt",sep="\t",header = T,nrows = 1)[-1]
classes=t(classes)
data=t(data)
data=data.frame(classes,data,row.names = colnames(raw))
data_raw = scale(data[,-1],center = TRUE,scale = TRUE)                
data_labels=data[,1]
train.label = as.factor(data_labels)
colors = rainbow(length(unique(train.label)))

dim(data)

unique(data_labels)

table(data_labels)

res.pca.poomeroy = FactoMineR::PCA(data_raw,scale.unit = FALSE,ncp = 34)

fviz_eig(res.pca.poomeroy,addlabels = FALSE,choice = "eigenvalue",ncp = Inf,geom = c("line"), linecolor = "red",ggtheme = theme(panel.background = element_rect(fill = 'grey', size = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(size = 0.1, linetype = "solid",colour = "black"))) + geom_abline(slope=0,intercept = 1) + geom_label_repel(aes(label=round(res.pca.poomeroy$eig[,"eigenvalue"],2))) 

data.frame(Pourcentage=res.pca.poomeroy$eig[,3],Cumulation=1:41) %>% ggplot(aes(x=Pourcentage,y=Cumulation)) + geom_bar(stat = "identity") + geom_hline(yintercept = 34,color="blue")

# La somme est nulle
print(sum(res.pca.poomeroy$eig[,"eigenvalue"] < 1))

fviz_pca_ind(res.pca.poomeroy,
             geom.ind = "point", # Montre les points seulement (mais pas le "text")
             col.ind = classes, # colorer by groups
             #palette = c("#00AFBB", "#E7B800", "#FC4E07","#f698db","#322c61","#1d2e36"),
             addEllipses = TRUE, # Ellipses de concentration
             legend.title = "Groups"
             )

fviz_pca_var(res.pca.poomeroy, labels=F,geom = "point",alpha=0.3)

# Cercle de corrélation
fviz_pca_var(res.pca.poomeroy, select.var = list(cos2=0.8), repel=T,labelsize=3)

dl=apply(res.pca.poomeroy$var$cos2[,c(1,2)],1,sum)

names_0.8 = names(dl[dl>0.8])

fviz_pca_biplot(res.pca.poomeroy,
                select.var=list(cos2=0.8),
                geom.ind = c("point"),
                geom.var = c("point","text"),
                repel=T,
                col.ind = classes)

corrplot.mixed(cor(res.pca.poomeroy$ind$coord, method = "spearman"), number.cex = .7)

## Affichage des dimension de l'ACP réduite
print(dim(res.pca.poomeroy$ind$coord))

#### LDA
dis2 <- MASS::lda(class~., data=cbind(res.pca.poomeroy$ind$coord,data.frame(class=data_labels)))

dis2$counts

## Pourcentage de discrimination de chaque axes
prop.lda = dis2$svd^2/sum(dis2$svd^2)
prop.lda*100

## Dimension de la LDA 
dim(dis2$scaling)

plot(dis2,col=c("red","green","blue","yellow")[as.numeric(data_labels)],dimen=2)

#Pomeroy MDS 
fit <- cmdscale( dist(x = data_raw,method = "euclidean"), k = 3,eig=TRUE)

### on affiche le GOF de la MDS Poomery
print(fit$GOF[1])

scatterplot3d(fit$points,pch = 16, color = colors[train.label],
              grid=TRUE, box=FALSE,main="MDS",xlab="",ylab="",zlab="")
legend("bottom", legend = levels(train.label),
      col =  unique(colors[train.label]), pch = 16, xpd = TRUE, horiz = TRUE,inset = -0.20)

data_lle_3d = LLE(data=data_raw, dim=3, k=41)

scatterplot3d(data_lle_3d,pch = 16, color = colors[train.label],
              grid=TRUE, box=FALSE,main="LLE",xlab="",ylab="",zlab="")
legend("bottom", legend = levels(train.label),
      col =  unique(colors[train.label]), pch = 16, xpd = TRUE, horiz = TRUE,inset = -0.25)

data_dim1to10_ISOMAP = Isomap(data=res.pca.poomeroy$ind$coord, dims=1:10, k=41, plotResiduals=TRUE,verbose =FALSE,mod=TRUE)

data_isomap_3d = Isomap(data=data_raw,dims=3, k=40,mod=FALSE)

plot(data_isomap_3d$dim3,t='n', main="ISOMAP",xlab="",ylab="")
text(data_isomap_3d$dim3, labels = train.label,col=colors[train.label])

scatterplot3d(data_isomap_3d$dim3,pch = 16, color = colors[train.label],
              grid=TRUE, box=FALSE,main="ISOMAP",xlab="",ylab="",zlab="")
legend("bottom", legend = levels(train.label),
      col =  unique(colors[train.label]), pch = 16, xpd = TRUE, horiz = TRUE,inset = -0.15)

rawdata= read.csv2("gordon-2002_database.txt",sep="\t",header = TRUE,row.names = 1)
data= read.csv2("gordon-2002_database.txt",sep="\t",header = FALSE,row.names = 1,skip = 2,dec=".")
classes=read.csv2("gordon-2002_database.txt",sep="\t",header = TRUE,nrows = 1)[-1]
classes=t(classes)
data=t(data)
data=data.frame(classes,data,row.names = colnames(raw))
data_raw = scale(data[,-1],center = TRUE,scale=TRUE)
data_labels=data[,1]
train.label = as.factor(data_labels)
colors = rainbow(length(unique(train.label)))

print(dim(data))

res.pca.gordon = FactoMineR::PCA(data_raw,ncp = Inf)

data.frame(Pourcentage=res.pca.gordon$eig[,3],Cumulation=1:180) %>% ggplot(aes(x=Pourcentage,y=Cumulation)) + geom_bar(stat = "identity") + geom_hline(yintercept = 137,color="blue")

fviz_pca_ind(res.pca.gordon,
             geom.ind = "point", # Montre les points seulement (mais pas le "text")
             col.ind = classes, # colorer by groups
             #palette = c("#00AFBB", "#E7B800", "#FC4E07","#f698db","#322c61","#1d2e36"),
             addEllipses = TRUE, # Ellipses de concentration
             legend.title = "Groups"
             )

# Cercle de corrélation
fviz_pca_var(res.pca.gordon, select.var = list(cos2=0.6), repel=T,labelsize=3)

dl=apply(res.pca.gordon$var$cos2[,c(1,2)],1,sum)
names_0.8 = names(dl[dl>0.6])

fviz_pca_biplot(res.pca.gordon,
                select.var=list(cos2=0.6),
                geom.ind = c("point"),
                geom.var = c("point","text"),
                repel=T,
                col.ind = classes)

### PCA sur les 137 premières composantes
res.pca.gordon = FactoMineR::PCA(data_raw,ncp = 137)

fviz_eig(res.pca.gordon,addlabels = FALSE,choice = "eigenvalue",ncp = Inf,geom = c("line"), linecolor = "red",ggtheme = theme(panel.background = element_rect(fill = 'grey', size = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(size = 0.1, linetype = "solid",colour = "black"))) + geom_abline(slope=0,intercept = 1) 

dis2 <- lda(res.pca.gordon$ind$coord, data_labels)

## Proportion des classes
print(dis2$prior)

## Count Gordon
print(dis2$counts)

# Séparation linéaire 
plot(dis2)

## Pourcentage de discrimination de chaque axes
prop.lda = dis2$svd^2/sum(dis2$svd^2)
prop.lda*100

### MDS Gordon
fit <- cmdscale( dist(x = data_raw,method = "euclidean"), k = 2,eig = TRUE)

# GOF de la MDS sur Gordon
print(fit$GOF[1])

plot(fit$points,t='n', main="MDS",xlab="",ylab="")
text(fit$points, labels = train.label,col=colors[train.label])

### LLE 2D 180 voisins
data_lle_2d = LLE(data=data_raw, dim=2, k=180)
plot(data_lle_2d,t='n', main="LLE",xlab="",ylab="")
text(data_lle_2d, labels = classes,col=colors[train.label])

## ISOMAP Plot Residual 1:10
### Coude k=2
data_dim1to10_ISOMAP = Isomap(data=data_raw, dims=1:10, k=179, plotResiduals=TRUE,verbose =F,mod=FALSE)

data_isomap_2d = Isomap(data=data_raw, dims=2, k=179,mod=FALSE)

plot(data_isomap_2d$dim2,t='n', main="ISOMAP",xlab="",ylab="")
text(data_isomap_2d$dim2, labels = classes,col=colors[train.label])
