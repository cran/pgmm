\name{pgmmEM}
\alias{pgmmEM}
\title{Model-Based Clustering & Classification Using PGMMs}
\description{
Carries out model-based clustering or classification using parsimonious Gaussian mixture models. AECM algorithms are used for parameter estimation. The BIC or the ICL is used for model selection.}
\usage{
pgmmEM(x,rG=1:2,rq=1:2,class=NULL,icl=FALSE,zstart=2,cccStart=TRUE,loop=3,zlist=NULL,
		modelSubset=NULL,seed=123456,tol=0.1,relax=FALSE)
}
\arguments{
  \item{x}{
A matrix or data frame such that rows correspond to observations and columns correspond to variables.
}
  \item{rG}{
The range of values for the number of components.
}
  \item{rq}{
The range of values for the number of factors.
}
  \item{class}{
If \code{NULL} then model-based clustering is performed. If a vector with length equal to the number of observations, then model-based classification is performed. In this latter case, the ith entry of \code{class} is either zero, indicating that the component membership of observation i is unknown, or it corresponds to the component membership of observation i. See Examples below.
}
  \item{icl}{
If \code{TRUE} then the ICL is used for model selection. Otherwise, the BIC is used for model selection.
}
  \item{zstart}{
A number that controls what starting values are used: (\code{1}) Random; (\code{2}) k-means; or 
(\code{3}) user-specified via \code{zlist}.
}
  \item{cccStart}{
If \code{TRUE} then random starting values are put through the CCC model and the resulting group memberships are used as starting values for the models specified in \code{modelSubset}. Only relevant for \code{zstart=1}. See Examples.
}
  \item{loop}{
A number specifying how many different random starts should be used. Only relevant for \code{zstart=1}.
}
  \item{zlist}{
A list comprising vectors of initial classifications such that \code{zlist[[g]]} gives the g-component starting values. Only relevant for \code{zstart=3}. See Examples.
}
  \item{modelSubset}{
A vector of strings giving the models to be used.
}
  \item{seed}{
A number giving the pseudo-random number seed to be used.
}
  \item{tol}{
A number specifying the epsilon value for the convergence criteria used in the AECM algorithms. For each algorithm, the criterion is based on the difference between the log-likelihood at an iteration and an asymptotic estimate of the log-likelihood at that iteration. This asymptotic estimate is based on the Aitken acceleration and details are given in the References. Values of \code{tol} greater than the default are not accepted.
}
  \item{relax}{
By default, the number of factors cannot exceed half the number of variables. Setting \code{relax=TRUE} relaxes this constraint and allows the number of factors to take any value up to three-quarters the number of variables.
}
}
\details{
The data \code{x} are either clustered using the PGMM approach of McNicholas & Murphy (2005, 2008, 2010) or classified using the method described by McNicholas (2010). In either case, all 12 covariance structures given by McNicholas & Murphy (2010) are available. Parameter estimation is carried out using AECM algorithms, as described in McNicholas et al. (2010). Either the BIC or the ICL is used for model-selection. The number of AECM algorithms to be run depends on the range of values for the number of components \code{rG}, the range of values for the number of factors \code{rq}, and the number of models in \code{modelSubset}. Starting values are very important to the successful operation of these algorithms and so care must be taken in the interpretation of results. 
}
\value{
An object of class \code{pgmm} is a list with components:
\item{map}{A vector of integers, taking values in the range \code{rG}, indicating the maximum \emph{a posteriori} classifications for the best model.}
\item{model}{A string giving the name of the best model.}
\item{g}{The number of components for the best model.}
\item{q}{The number of factors for the best model.}
\item{zhat}{A matrix giving the raw values upon which \code{map} is based.}
\item{load}{The factor loadings matrix (Lambda) for the best model.}
\item{noisev}{The Psi matrix for the best model.}
\item{plot_info}{A list that stores information to enable \code{plot}.}
\item{summ_info}{A list that stores information to enable \code{summary}.}
In addition, the object will contain one of the following, depending on the value of \code{icl}.
\item{bic}{A number giving the BIC for each model.}
\item{icl}{A number giving the ICL for each model.}
}
\note{
Dedicated \code{print}, \code{plot}, and \code{summary} functions are available for objects of class \code{pgmm}.
}
\author{
Paul D. McNicholas [aut, cre], Aisha ElSherbiny [aut], K. Raju Jampani [ctb], Aaron McDaid [aut], Brendan Murphy [aut], Larry Banks [ctb]

Maintainer: Paul D. McNicholas <mcnicholas@math.mcmaster.ca>
}
\references{
Paul D. McNicholas and T. Brendan Murphy (2010). Model-based clustering of microarray expression data via latent Gaussian mixture models. \emph{Bioinformatics} \bold{26}(21), 2705-2712.

Paul D. McNicholas (2010). Model-based classification using latent Gaussian mixture models. \emph{Journal of Statistical Planning and Inference} \bold{140}(5), 1175-1181.

Paul D. McNicholas, T. Brendan Murphy, Aaron F. McDaid and Dermot Frost (2010). Serial and parallel implementations of model-based clustering via parsimonious Gaussian mixture models. \emph{Computational Statistics and Data Analysis} \bold{54}(3), 711-723.

Paul D. McNicholas and T. Brendan Murphy (2008). Parsimonious Gaussian mixture models. \emph{Statistics and Computing} \bold{18}(3), 285-296.

Paul D. McNicholas and T. Brendan Murphy (2005). Parsimonious Gaussian mixture models. Technical Report 05/11, Department of Statistics, Trinity College Dublin.
}
 \examples{
  \dontrun{
# Wine clustering example with three random starts and the CUU model.
 data("wine")
 x<-wine[,-1]
 x<-scale(x)
 wine_clust<-pgmmEM(x,rG=1:4,rq=1:4,zstart=1,loop=3,modelSubset=c("CUU"))
 table(wine[,1],wine_clust$map)

# Wine clustering example with custom starts and the CUU model.
 data("wine")
 x<-wine[,-1]
 x<-scale(x)
 hcl<-hclust(dist(x)) 
 z<-list()
 for(g in 1:4){ 
	 z[[g]]<-cutree(hcl,k=g)
 } 
 wine_clust2<-pgmmEM(x,1:4,1:4,zstart=3,modelSubset=c("CUU"),zlist=z)
 table(wine[,1],wine_clust2$map)
 print(wine_clust2)
 summary(wine_clust2)

# Olive oil classification by region (there are three regions), with two-thirds of
# the observations taken as having known group memberships, using the CUC, CUU and 
# UCU models. 
 data("olive") 
 x<-olive[,-c(1,2)] 
 x<-scale(x) 
 cls<-olive[,1]
 for(i in 1:dim(olive)[1]){
	 if(i\%\%3==0){cls[i]<-0}
 }
 olive_class<-pgmmEM(x,rG=3:3,rq=4:6,cls,modelSubset=c("CUC","CUU", 
  "CUCU"),relax=TRUE)
 cls_ind<-(cls==0) 
 table(olive[cls_ind,1],olive_class$map[cls_ind])

# Another olive oil classification by region, but this time suppose we only know
# two-thirds of the labels for the first two areas but we suspect that there might 
# be a third or even a fourth area. 
 data("olive") 
 x<-olive[,-c(1,2)] 
 x<-scale(x) 
 cls2<-olive[,1]
 for(i in 1:dim(olive)[1]){
   if(i\%\%3==0||i>420){cls2[i]<-0}
 }
 olive_class2<-pgmmEM(x,2:4,4:6,cls2,modelSubset=c("CUU"),relax=TRUE)
 cls_ind2<-(cls2==0) 
 table(olive[cls_ind2,1],olive_class2$map[cls_ind2])
 
# Coffee clustering example using k-means starting values for all 12
# models with the ICL being used for model selection instead of the BIC.
 data("coffee")
 x<-coffee[,-c(1,2)]
 x<-scale(x)
 coffee_clust<-pgmmEM(x,rG=2:3,rq=1:3,zstart=2,icl=TRUE)
 table(coffee[,1],coffee_clust$map)
 plot(coffee_clust)
 plot(coffee_clust,onlyAll=TRUE)
  }
  
# Coffee clustering example using k-means starting values for the UUU model, i.e., the
# mixture of factor analyzers model, for G=2 and q=1.
 data("coffee")
 x<-coffee[,-c(1,2)]
 x<-scale(x)
 coffee_clust_mfa<-pgmmEM(x,2:2,1:1,zstart=2,modelSubset=c("UUU"))
 table(coffee[,1],coffee_clust_mfa$map)
 }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{cluster}
\keyword{classif}
