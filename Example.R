##################################################################################################
# This is an example of running MLVPCA                                                           #
# The data:                                                                                      #
#     The sample data has 5 variates, 50 time points/variate,                                    #
#     100 subjects, and 5 repeated measures/subject.                                             #
# Model assumption:                                                                              #
#     Let's assume a two-way nested functional ANOVA model for the data,                         #                                      # 
#      with unstructured correlation structure among repeated measures.                          #
# Goal:                                                                                          #
#     Identify between subject level principal components and within subject                     #
#     level principal components, which highlight important time points and important variables. #
# Output of interest:                                                                            #
#     1. Eigenfunctions: directions that best explains the functional variation of the data      #
#     2. Principal component scores: reduced-dimensional predictors that can be further used to  #
#        predict clinical outcomes.                                                              #
##################################################################################################

load("data1.RData")

source("MLVPCA.r")
source("fun_GetCovarianceEst.r")
source("fun_SelectTuningPar.r")
source("fun_SeqADMM.r")
source("fun_other.r")

install.packages("Matrix")
install.packages("pracma")
install.packages("magic")

###############################################################################################
# Arguments:                                                                                  #
# data: a dataset where each row is individual time series that concatenates m length-p       #
#       vectors. Each column contains L2 repeated measures nested within L1 subjects.         #
# m: number of variates                                                                       #
# L1: number of subjects                                                                      #
# L2: number of repeated measures                                                             #
# k: the number of principal components to estimate for each level                            #
# model: the type of model to use. Default is "2WayNested.spatial", can also take "2WayNested"#
# parallel: whether to use parallel computing                                                 #
# no_cores: number of cores to take if using parallel computing                               #
# maxiter: maximum iteration to allow for                                                     #
# nfold: number of folds to divide the data for cross-validation selecting tuning parameter   #
# FOR MORE OPTIONS PLEASE SEE MANUAL FOR MORE DETAILS                                         #
###############################################################################################

result1 <- MLVPCA(data1,m=3, L1=100, L2=5, k=3, model="2WayNested.spatial",
                                  parallel = F, no_cores=1, maxiter=100, nfold=5)

##############################
#                            #
# Result 1 -- Eigenfunctions #
#   (shown in plots)         #
##############################

par(mfrow=c(2,3))

for (j in 1:3){
  phi <- result1$phi.z[,j]
  p <- length(phi)
  plot(phi, type="l", col = "red", cex=0.75, pch=20, cex.lab = 1.5,
       xlim=c(0,p), ylab=bquote(phi[.(j)]), xlab = "", xaxt = 'n')
  abline(v=seq(0,p,by=p/(3)), col=rep("black", 4),lwd=1)
  axis(side=1,at=c(25,75,125),labels=c("var1","var2","var3"))
  FVE.z=round(result1$FVE.z,3)
  title(bquote(paste( "Subject level ",phi[.(j)], ", explain " , .(FVE.z[j]), " variance")), outer=F,cex=2)
}

for (j in 1:3){
  phi <- result1$phi.w[,j]
  p <- length(phi)
  plot(phi, type="l", col = "red", cex=0.75, pch=20, cex.lab = 1.5,
       xlim=c(0,p), ylab=bquote(phi[.(j)]), xlab = "", xaxt = 'n' )
  abline(v=seq(0,p,by=p/(3)), col=rep("black", 4),lwd=1)
  axis(side=1,at=c(25,75,125),labels=c("var1","var2","var3"))
  FVE.w=round(result1$FVE.w,3)
  title(bquote(paste( "Electrode level ",phi[.(j)], ", explain " , .(FVE.w[j]), " variance")), outer=F,cex=2)
}

##############################
#                            #
# Result 2 -- PC scores      #
#                            #
##############################

# Subject level PC scores: each subject has 3 scores for PC1, PC2, and PC3

Zscore=unique(result1$PCscore.z[,1:3])
colnames(Zscore)=c("z1","z2","z3")
head(Zscore)

#           z1          z2         z3
# [1,]  0.38074056  0.27595976  0.4970751
# [2,]  0.93503404 -0.43811180  0.1603260
# [3,] -0.07241067  0.49980037  0.2631204
# [4,] -0.10044029 -0.62224119 -0.6081468
# [5,] -0.17498456  0.02835256  0.1670053
# [6,] -1.03942810 -0.05912763  0.2519630



# Within Subject level PC scores: each subject has 5 scores for per component, 15 scores in total

Wscore=cbind(t(matrix(result1$PCscore.w[,1],nrow=5)),
             t(matrix(result1$PCscore.w[,2],nrow=5)),
             t(matrix(result1$PCscore.w[,3],nrow=5)))
colnames(Wscore)=c(paste("w1.",as.character(1:5),sep=""),
                   paste("w2.",as.character(1:5),sep=""),
                   paste("w3.",as.character(1:5),sep=""))
head(Wscore)

#          w1.1       w1.2       w1.3       w1.4       w1.5        w2.1        w2.2       w2.3
# [1,] -1.55709716 -1.3437967 -1.4695118 -0.5228291 -0.1733608 -0.14772516  0.06280980 -0.4894211
# [2,]  0.77397074 -0.1331736  0.6428406  0.7990022  0.9739577 -0.02669401  0.27263267 -0.4047446
# [3,] -0.25366500 -0.3966171 -0.1977496  0.1114438  1.0743954  0.13374195 -0.09409238 -0.2100285
# [4,]  0.28369566  1.1094355 -0.7646597 -1.7083238 -0.2008515  0.82643981  0.32041283  1.2457603
# [5,]  0.14269262 -0.5528265  0.2486981  0.7448324  2.0617321 -0.69211861 -1.49163288 -0.4843820
# [6,] -0.01652084 -0.8477890  0.1165059 -0.9342613  0.5605362 -0.28021014  0.28346981  0.6631859
#          w2.4        w2.5       w3.1        w3.2        w3.3        w3.4       w3.5
# [1,] -0.03115104 -0.05123226  0.1120238  0.10561988 -0.09926533 -0.04296053  0.8287233
# [2,]  0.40561371  0.05332901  0.4313566  0.49167785  0.02314517 -0.01681057  0.5105388
# [3,]  0.06106034 -0.45813193 -0.9059911 -0.12462067  0.28604317  0.41218880 -0.1186995
# [4,]  0.61349583  0.10692350  0.2680816  0.63505674  0.96264927  0.54239414 -0.3739689
# [5,] -0.34815172  0.39085420 -0.2425786  0.01673745  0.27386781  0.06662542 -0.2876591
# [6,]  0.82478833  0.68425741 -0.1724768  0.19644754 -0.24904077 -0.23413316 -0.5819663
