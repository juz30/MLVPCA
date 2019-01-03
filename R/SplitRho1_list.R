

SplitRho1_list=function(X,op){


  if (op$model=="2WayNested"|op$model=="2WayNested.spatial"){

    L1=op$levels[1]
    L2=op$levels[2]
    nsplit=floor(L1/op$nfold)
    # starting value
    Vsplit.w=0
    Vsplit.z=0
    i1.w=0
    i1.z=0
    nrho1=length(op$rho1Seq.list$rho1Seq.w) # WE REQUEST SEQUENCE OF Z AND W  ARE THE SAME.

    # To get CV <S,H> for each of the 10 candidate rho, using parallel computing or not

    if (option$parallel==T){

      # parallel computing
      no_cores=op$no_cores
      cl=makeCluster(no_cores,type="FORK")

      EachElement=function(x){
        rho1.w=op$rho1Seq.list$rho1Seq.w[x]
        rho1.z=op$rho1Seq.list$rho1Seq.z[x]
        V.w=0
        V.z=0
        # For each rho, do cv nfold times, and aggregate test tr(t(H)*S) to get V. Then compare V for every rho
        # Split data by L1 (subject), if has more levels, split the data by the uppest level (to make the data structure complete within train and test)
        for (ifold in 1:op$nfold){
          idtest=((ifold-1)*nsplit*L2+1):(ifold*nsplit*L2)
          Xtrain=X[-idtest,]
          Xtest=X[idtest,]
          G.train=MultilevelS(Xtrain,L1-nsplit,L2,op)
          G.test=MultilevelS(Xtest,nsplit,L2,op)

          # Within subject
          K.train.w=G.train$G.w-rho1.w*op$SmoothD
          projH.w=seqADMM_0330(K.train.w, 1, op$PrevPi, op$PrevPi_d,0,0)
          K.test.w=G.test$G.w
          V.w=V.w+sum(projH.w*K.test.w)
          # Between subject
          K.train.z=G.train$G.z-rho1.z*op$SmoothD
          projH.z=seqADMM_0330(K.train.z, 1, op$PrevPi, op$PrevPi_d,0,0)
          K.test.z=G.test$G.z
          V.z=V.z+sum(projH.z*K.test.z)
        }
        cat("current total sum",V.w, V.z, "\t")
        return(c(V.w,V.z))
      }

      Vmatrix=matrix(parSapply(cl, 1:nrho1, EachElement),byrow = F, ncol=nrho1, nrow=2)
      stopCluster(cl)

      i1.w=which.max(Vmatrix[1,])
      i1.z=which.max(Vmatrix[2,])

    } else{

      for (i in 1:nrho1){

      rho1.w=op$rho1Seq.list$rho1Seq.w[i]
      rho1.z=op$rho1Seq.list$rho1Seq.z[i]
      V.w=0
      V.z=0
      # For each rho, do cv nfold times, and aggregate test tr(t(H)*S) to get V. Then compare V for every rho
      # Split data by L1 (subject), if has more levels, split the data by the uppest level (to make the data structure complete within train and test)
      for (ifold in 1:op$nfold){
        idtest=((ifold-1)*nsplit*L2+1):(ifold*nsplit*L2)
        Xtrain=X[-idtest,]
        Xtest=X[idtest,]
        G.train=MultilevelS(Xtrain,L1-nsplit,L2,op)
        G.test=MultilevelS(Xtest,nsplit,L2,op)

        # Within subject
        K.train.w=G.train$G.w-rho1.w*op$SmoothD
        projH.w=seqADMM_0330(K.train.w, 1, op$PrevPi, op$PrevPi_d,0,0)
        K.test.w=G.test$G.w
        V.w=V.w+sum(projH.w*K.test.w)
        # Between subject
        K.train.z=G.train$G.z-rho1.z*op$SmoothD
        projH.z=seqADMM_0330(K.train.z, 1, op$PrevPi, op$PrevPi_d,0,0)
        K.test.z=G.test$G.z
        V.z=V.z+sum(projH.z*K.test.z)
      }
        cat("current total sum",V.w, V.z, "\t")

        # If the V for the current rho is larger than the V for the previous rho, we replace Vsplit and il (to track down V and i).
      if (V.w>=Vsplit.w){
        Vsplit.w=V.w
        i1.w=i
      }
      if (V.z>=Vsplit.z){
        Vsplit.z=V.z
        i1.z=i
      }
    }


    }

      # Select best rho based on i1.w and i1.z

      if (i1.w<nrho1){
        rho_1.w=op$rho1Seq.list$rho1Seq.w[i1.w+1]
      }else {
        rho_1.w=op$rho1Seq.list$rho1Seq.w[i1.w]
      }

      if (i1.z<nrho1){
        rho_1.z=op$rho1Seq.list$rho1Seq.z[i1.z+1]
      }else{
        rho_1.z=op$rho1Seq.list$rho1Seq.z[i1.z]
      }

      # rho_1.w=op$rho1Seq.list$rho1Seq.w[i1.w]
      # rho_1.z=op$rho1Seq.list$rho1Seq.z[i1.z]

    return(rho_1=list("rho_1.w"=rho_1.w, "rho_1.z"=rho_1.z))

  }



}
