mi_ap_cond<-function(data1,data2,data3){
  # CMI(TF;TG|M)
  # all data inputs should be a vector class
  # data1 for TF
  # data2 for TG
  # data3 for M
  DIM = 3
  N1 <- length(data1)
  N2 <- length(data2)
  N3 <- length(data3)
  if (N1!=N2 & N1!=N2 & N2!=N3) {
    print("Size should be the same")
    break
  }
  N <- N1
  
  B1 <-sort(data1,index.return = T)
  B2 <-sort(data2,index.return = T)
  B3 <-sort(data3,index.return = T)
  b1<-B1$ix
  b2<-B2$ix
  b3<-B3$ix  
  idat <- cbind(b1,b2,b3)
  ydat=idat
  
  for(d in 1:DIM){
    ydat[idat[,d],d]=c(1:N);  
  }
  ddim=2^DIM; dim2=2*DIM;
  rm(idat)
  
  poc<-c()
  kon<-c()
  xcor=0; npar=1; poc[1]=1; kon[1]=N; poradi=c(1:N);
  NN=rep(0,ddim); marg=matrix(0,nrow = 8*ddim,ncol = dim2);
  marg[1,]=c(rep(1,DIM),   N*rep(1,DIM));
  Imm=rbind(0,1);  
  
  for(d in 2:DIM){
    Imm=rbind(cbind(matrix(0,nrow = nrow(Imm),ncol = ncol(Imm)) , Imm), cbind(matrix(1,nrow = nrow(Imm),ncol = ncol(Imm)),  Imm))
  }
  Imm<-Imm[,-1]
  chi2=c(0,   7.81,   13.9,   25.0,   42.0);
  run=0;
  while(npar>0){
    run=run+1;
    apoc=poc[npar]; akon=kon[npar];
    apor=poradi[apoc:akon]; Nex=length(apor);
    ave=floor((marg[npar,1:DIM]+marg[npar,(DIM+1):dim2])/2);
    J=(ydat[apor,]<= matrix(1,nrow = Nex, ncol =1)%*%ave)*1;
    I=matrix(0,nrow = Nex, ncol = ddim);
    amarg=matrix(1,nrow = ddim,ncol = 1) %*% marg[npar,];
    for(d in 1:ddim){
      I[,d]=matrix(1,nrow = Nex,1);
      for(k in 1:DIM){
        if(Imm[d,k]!=0){
          I[,d]=I[,d] & !J[,k];
          amarg[d,k]=ave[k]+1;
        }else{
          I[,d]=I[,d] & J[,k];
          amarg[d,k+DIM]=ave[k];
        }
        
      }
    }
    
    NN=apply(I,2,sum);
    tst = ddim*sum((NN - Nex/ddim * rep(1,ddim))^2)/Nex;
    
    if(tst>chi2[DIM] | run==1){
      npar=npar-1;
      for(ind in 1:ddim){
        if(NN[ind] > ddim){
          npar=npar+1;
          akon=apoc+NN[ind]-1;
          poc[npar]=apoc; kon[npar]=akon;
          marg[npar,]=amarg[ind,];
          poradi[apoc:akon]=apor[which(I[,ind]!=0)];
          apoc=akon+1;          
        }else{
          if(NN[ind] > 0){
            Nxx=prod(amarg[ind,(DIM+1):dim2]-amarg[ind,1:DIM]+rep(1,DIM));
            Nz=amarg[ind,6]-amarg[ind,3]+1;
            Jx = ((ydat[,1]>=amarg[ind,1]) & (ydat[,1]<=amarg[ind,4]))*1;
            Jy = ((ydat[,2]>=amarg[ind,2]) & (ydat[,2]<=amarg[ind,5]))*1;
            Jz = ((ydat[,3]>=amarg[ind,3]) & (ydat[,3]<=amarg[ind,6]))*1;
            Nxz=sum(Jx & Jz);
            Nyz=sum(Jy & Jz);
            
            cond = (NN[ind]*Nz)/(Nxz*Nyz);
            if(is.infinite(cond)){
              cond =1
            } 
            if(cond==0){
              cond =1
            }              
            xcor=xcor+NN[ind]*log(cond);            
          }          
        }
        
      }
    }else{
      Nxx=prod(marg[npar,(DIM+1):dim2]-marg[npar,1:DIM] + rep(1,DIM));
      Nz=marg[npar,6]-marg[npar,3]+1;
      Jx = ((ydat[,1]>=marg[npar,1]) & (ydat[,1]<=marg[npar,4]))*1;
      Jy = ((ydat[,2]>=marg[npar,2]) & (ydat[,2]<=marg[npar,5]))*1;
      Jz = ((ydat[,3]>=marg[npar,3]) & (ydat[,3]<=marg[npar,6]))*1;
      Nxz=sum(Jx & Jz); 
      Nyz=sum(Jy & Jz); 
      cond = (Nex*Nz)/(Nxz*Nyz);
      if(is.infinite(cond)){
        cond =1
      } 
      if(cond==0){
        cond =1
      }       
      xcor=xcor+Nex*log(cond);
      npar=npar-1;      
    }    
  } 
  return(xcor/N)
}


