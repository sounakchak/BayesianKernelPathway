
#function
NPKM<-function(dat1,genexp,Snki,Rnk,pathprob,geneprob,n,G,abeta,aalpha,a,b,mu,eta,drt){
  library(tmvtnorm)
  is.wholenumber <- function(z, tol = .Machine$double.eps^0.5)  abs(z - round(z)) < tol
  c=rep(1,n)
  yc=cbind(c,dat1)
  yc=t(yc)
  n=length(yc[1,])
  z=y=(as.matrix(yc[2,]))
  y=t(y)
  censor=which(yc[1,]==0)
  k=length(Snki[,1])
  p=length(Snki[1,])
  pwygmno=rowSums(Snki!=0)
  longestpwy=max(rowSums(Snki!=0))
  
  #initial pathway and genes
  repeat{
    PHI=0
    pwygene=matrix(0,nrow=longestpwy,ncol=k)
    gm=rep(0,length(Snki[1,]))
    
    phi=rbinom(k,1,pathprob)
    for(i in 1:k){
      PHI[((i-1)*n+1):(i*n)]=rep(phi[i],n)
    }
    inipwy=which(phi==1)
    for (i in 1:length(inipwy)){
      if (pwygmno[inipwy][i]>round(p*geneprob/sum(phi)-2)){
        pwygene[1:(round(p*geneprob/sum(phi)-2)),inipwy[i]]=sample(which(Snki[inipwy[i],]!=0),round(p*geneprob/sum(phi)-2))
        gm[pwygene[1:(round(p*geneprob/sum(phi)-2)),inipwy[i]]]=1
      }else{
        pwygene[1:pwygmno[inipwy][i],inipwy[i]]=which(Snki[inipwy[i],]!=0)
        gm[pwygene[1:pwygmno[inipwy][i],inipwy[i]]]=1
      }
    }
    
    pwygene=apply(pwygene,2,sort)
    inigmind=which(gm==1)
    for (i in 1:sum(gm)){
      if (sum(Snki[,inigmind[i]]!=0)>1){
        h=which(Snki[,inigmind[i]]!=0)
        h1=h[which(match(h,inipwy,nomatch=0)!=0)]
        if (length(h1)>1){
          maxg=max(gmnumber=colSums(pwygene[,h1]!=0))
          h2=h1[-which(pwygene[(longestpwy-maxg+1):longestpwy,h1]==inigmind[i],arr.ind = T)[,2]]
          pwygene[1,h2]=inigmind[i]
          pwygene[,h2]=apply(as.matrix(pwygene[,h2]),2,sort)
        }
      }
    }
    cr2=1
    subpwygene=pwygene[(longestpwy-5):longestpwy,which(phi!=0)]
    for (i in 1:(dim(subpwygene)[2]-1)){
      for (j in (i+1):dim(subpwygene)[2]){
        if (sum(subpwygene[,i]==subpwygene[,j])==6){
          cr2=0
          break
        }
      }
      if (cr2==0){break}
    }
    if (cr2==1){break}
  }#repeat ends
  
  #begin
  
  simphi=matrix(,nrow=k,ncol=G)
  simphi[,1]=phi
  simgm=matrix(,nrow=p,ncol=G)
  simgm[,1]=gm
  aK=matrix(0,nrow=n,ncol=n)
  theta=rep(2,k) #######theta's initial value#####
  simtheta=matrix(0,nrow=k,ncol=G)
  simtheta[,1]=theta
  simy=matrix(,nrow=n,ncol=G)
  simy[,1]=y
  lower=z[censor]
  tdf=a*2
  s2=b*2/tdf
  u=rep(0,length(lower))
  
  oldK=matrix(0,nrow=n, ncol=n*sum(phi))
  oldpwy=which(phi!=0)
  
  for (i in 1:length(oldpwy)){
    x=genexp[pwygene[,oldpwy[i]],]
    for(j in 1:n){
      for(l in j:n){
        oldK[j,(l+(i-1)*n)]=oldK[l,(j+(i-1)*n)]=exp(-sqrt(sum((x[,j]-x[,l])^2))/theta[oldpwy[i]])
      }
    }
  }
  oldKPHI=(oldK)
  J=matrix(1,nrow=n,ncol=n)
  oldSigma=diag(n)+aalpha*J+abeta*oldKPHI%*%t(oldKPHI)
  involdSigma=solve(oldSigma)
  
  R=as.matrix(Rnki)
  oldpwygene=pwygene
  
  logr2=-(n/2+a)*log(y%*%involdSigma%*%t(y)+2*b)+sum(simphi[,1])*log(pathprob)+(k-sum(simphi[,1]))*log(1-pathprob)+mu*sum(gm)+eta*t(gm)%*%R%*%gm
  
  ##MCMC
  for (g in 2:G){
    
    #proy=y
    #proy[censor]=rtmvt(1,u,(s2*oldSigma[censor,censor]),tdf,lower,algorithm="gibbs")
    #logr1=-(n/2+a)*log(2*b+proy%*%involdSigma%*%t(proy))+sum(simphi[,g-1])*log(pathprob)+(k-sum(simphi[,g-1]))*log(1-pathprob)+mu*sum(simgm[,g-1])+eta*t(simgm[,g-1])%*%R%*%simgm[,g-1]
    #logr=logr1-logr2
    #if((log(runif(1)))<logr) {
    #   y=proy
    #   logr2=logr1
    #   simy[,g]=proy
    #}else{simy[,g]=y}
    
    
    repeat{
      move=sample(1:3,1)
      ###########move 1########
      if (move==1){
        if (sum(simphi[,g-1])>2) {ar=sample(1:2,1)}else{ar=1}
        ########add pahtway and a gene##########
        if (ar==1){
          ######identify subset of pathway no in the model nor of their genes####
          #subgm=(simgm[,g-1]==0) #find out the gene that isn't in the model
          subSnki=Snki[-which(simphi[,g-1]==1),] #SubSnki is the Snki only with non-selected pwy
          subSnki=subSnki[-which(colSums(t(subSnki)*simgm[,g-1])!=0),] # only the non-selected pwy with no genes in the model
          if (dim(subSnki)[1]>0){
            apwyind=which(rownames(Snki)==sample(rownames(subSnki),1)) #randomly select 1 pwy
            ageneind=sample(which(Snki[apwyind,]!=0),1) #randomly choose one gene
            prophi=simphi[,g-1]
            prophi[apwyind]=1 #prophi with selected 1 pwy
            progm=simgm[,g-1]
            progm[ageneind]=1 #progm with selected 1 gene
            propwygene=oldpwygene
            propwygene[longestpwy,apwyind]=ageneind # add the information to pwygene
            subSnki1=Snki[which(simphi[,g-1]==1),] #subSnki1 is the Snki only with the selected pwy
            proK=oldK
            if (sum(subSnki1[,ageneind])>0){
              subpwy1=which(match(rownames(Snki),rownames(subSnki1)[which(subSnki1[,ageneind]!=0)],nomatch=0)!=0) #indicator of the selected pwy that has the same new selected gene
              for (i in 1:length(subpwy1)){
                propwygene[1,subpwy1[i]]=ageneind # add the gene in the selected pwy
                propwygene[,subpwy1[i]]=sort(propwygene[,subpwy1[i]]) # sort the gene in those pwy
                x=genexp[propwygene[,subpwy1[i]],]
                loc=which(subpwy1[i]==which(simphi[,g-1]==1))
                for(j in 1:n){
                  for(l in j:n){
                    proK[j,(l+(loc-1)*n)]=proK[l,(j+(loc-1)*n)]=exp(-sqrt(sum((x[,j]-x[,l])^2))/theta[subpwy1[i]])
                  }
                }
              }
            }
            
            x=genexp[ageneind,]
            for(j in 1:n){
              for(l in j:n){
                aK[j,l]=aK[l,j]=exp(-sqrt(sum((x[,j]-x[,l])^2))/theta[apwyind])
              }
            }
            
            proKPHI=matrix(0,nrow=n,ncol=sum(prophi)*n)
            aKloc=which(apwyind==which(prophi==1))
            if (aKloc==1) {proKPHI=cbind(aK,proK)}
            if (aKloc==sum(prophi)) {proKPHI=cbind(proK,aK)}
            if (aKloc>1 & aKloc<sum(prophi)){
              proKPHI[,1:((aKloc-1)*n)]=proK[,1:((aKloc-1)*n)]
              proKPHI[,((aKloc-1)*n+1):(aKloc*n)]=aK
              proKPHI[,(aKloc*n+1):(sum(prophi)*n)]=proK[,-(1:((aKloc-1)*n))]
            }
            
            
            
            #####Calculate the number of patyway belongs to cond1 and condId1####
            subpwy2=which(colSums(propwygene!=0)==1) # pathway that has 1 gene
            removeind=0
            for (l in 1:length(subpwy2)){
              pwywsg=which(propwygene==propwygene[longestpwy,subpwy2[l]],arr.ind=TRUE)[,2] #pathway with the same remove gene
              if (length(pwywsg)>1){
                pgsrm=propwygene[,pwywsg] #pathway has the same gene as the remove pathway
                pgsrm[which(pgsrm==propwygene[longestpwy,subpwy2[l]])] =0         #remove gene from pathyways
                pgsrmed=as.matrix(pgsrm[,-which(colSums(pgsrm)==0)])#delete the rm pathway
                for (i in 1:length(pgsrmed[1,])){
                  for (j in 1:k){
                    if (sum(sort(pgsrmed[,i])!=propwygene[,j])==0) {
                      removeind[l]=0
                      break
                    } else{removeind[l]=1}
                  }
                  if (removeind[l]==0) {break}
                }
              }else {removeind[l]=1}
            }
            
            
            
            
            proSigma=diag(n)+aalpha*J+abeta*proKPHI%*%t(proKPHI)
            invproSigma=solve(proSigma)
            logr1=-(n/2+a)*log(2*b+y%*%invproSigma%*%t(y))+sum(prophi)*log(pathprob)+(k-sum(prophi))*log(1-pathprob)+mu*sum(progm)+eta*t(progm)%*%R%*%progm
            logr11=logr1+log(pathprob)+log(dim(subSnki)[1])-log(sum(removeind))
            logr=logr11-logr2
            if((log(runif(1)))<logr) {
              simphi[,g]=prophi
              simgm[,g]=progm
              logr2=logr1
              oldK=proKPHI
              oldpwygene=propwygene
              oldSigma=proSigma
              involdSigma=invproSigma
            }else{
              simphi[,g]=simphi[,g-1]
              simgm[,g]=simgm[,g-1]
            }
            break
          }
        }else{
          ######remove pathway and its only gene###########
          #####Calculate the number of patyway belongs to cond1 and condId1####
          subpwy2=which(colSums(oldpwygene!=0)==1) # pathway that has 1 gene
          removeind=0
          if (length(subpwy2)>0){
            for (l in 1:length(subpwy2)){
              pwywsg=which(oldpwygene==oldpwygene[longestpwy,subpwy2[l]],arr.ind=TRUE)[,2] #pathway with the same remove gene
              if (length(pwywsg)>1){
                pgsrm=oldpwygene[,pwywsg] #pathway has the same gene as the remove pathway
                pgsrm[which(pgsrm==oldpwygene[longestpwy,subpwy2[l]])] =0         #remove gene from pathyways
                pgsrmed=as.matrix(pgsrm[,-which(colSums(pgsrm)==0)])#delete the rm pathway
                for (i in 1:length(pgsrmed[1,])){
                  for (j in 1:k){
                    if (sum(sort(pgsrmed[,i])!=oldpwygene[,j])==0) {
                      removeind[l]=0
                      break
                      
                    } else{removeind[l]=1}
                  }
                  if (removeind[l]==0) {break}
                }
              }else {removeind[l]=1}
            }
            
            
            if (sum(removeind)>=1){
              subpwy3=subpwy2[which(removeind==1)] #the pwy that can be removed fix the cond1 and condId1
              if(length(subpwy3)>1){rpwyind=sample(subpwy3,1)}else{rpwyind=subpwy3} #randomly choose one pwy
              rgeneind=oldpwygene[longestpwy,rpwyind]  #choose the removed gene
              prophi=simphi[,g-1]
              prophi[rpwyind]=0
              progm=simgm[,g-1]
              progm[rgeneind]=0
              propwygene=oldpwygene
              propwygene[longestpwy,rpwyind]=0
              pwywsg=which(propwygene==rgeneind,arr.ind=TRUE)[,2] # other pwy that has the removed gene
              rloc=which(rpwyind==which(simphi[,g-1]==1))
              proK=oldK[,-(((rloc-1)*n+1):(rloc*n))]
              if (length(pwywsg)>0){
                for (i in 1:length(pwywsg)){
                  propwygene[which(propwygene[,pwywsg[i]]==rgeneind),pwywsg[i]]=0 # rmove the gene
                  propwygene[,pwywsg[i]]=sort(propwygene[,pwywsg[i]])   # sort that removed gene pwy
                  x=genexp[propwygene[which(propwygene[,pwywsg[i]]!=0),pwywsg[i]],]
                  loc=which(pwywsg[i]==which(prophi==1))
                  for(j in 1:n){
                    for(l in j:n){
                      proK[j,(l+(loc-1)*n)]=proK[l,(j+(loc-1)*n)]=exp(-sqrt(sum((x[,j]-x[,l])^2))/theta[pwywsg[i]])
                    }
                  }
                }
              }
              
              
              
              proKPHI=(proK)
              proSigma=diag(n)+aalpha*J+abeta*proKPHI%*%t(proKPHI)
              invproSigma=solve(proSigma)
              logr1=-(n/2+a)*log(2*b+y%*%invproSigma%*%t(y))+sum(prophi)*log(pathprob)+(k-sum(prophi))*log(1-pathprob)+mu*sum(progm)+eta*t(progm)%*%R%*%progm
              logr11=logr1-log(pathprob)-log(length(which(colSums(propwygene!=0)==0)))+log(sum(removeind))
              logr=logr11-logr2
              if((log(runif(1)))<logr) {
                simphi[,g]=prophi
                simgm[,g]=progm
                logr2=logr1
                oldK=proKPHI
                oldpwygene=propwygene
                oldSigma=proSigma
                involdSigma=invproSigma
                break
              }else{
                simphi[,g]=simphi[,g-1]
                simgm[,g]=simgm[,g-1]
                break
              }
            }
          }
        }
      }
      
      if(move==2){
        if (sum(simgm[,g-1])>3) {ar=sample(1:2,1)}else{ar=1}
        if (ar==1){
          subpwy1=which(colSums(oldpwygene!=0)!=0)#subways selected in the model
          subpwy2=subpwy1[which(colSums(as.matrix(oldpwygene[,subpwy1])!=0)<pwygmno[subpwy1])] # set G
          if (length(subpwy2)>1){apwyind=sample(subpwy2,1)}else{apwyind=subpwy2}#randomly select a pwy
          gmin=oldpwygene[which(oldpwygene[,apwyind]!=0),apwyind] # gene of that pwy already in the model
          gmoff=which(match(which(Snki[apwyind,]!=0),gmin,nomatch=0)==0)#gene that are available to choose
          ageneind=sample(which(Snki[apwyind,]!=0)[gmoff],1)#randomly choose a gene
          prophi=simphi[,g-1]
          progm=simgm[,g-1]
          progm[ageneind]=1
          pwysgm=which(Snki[,ageneind]!=0)#pwy that have the same choosen gene
          pwysgm=pwysgm[which(match(pwysgm,subpwy2,nomatch=0)!=0)] # the previous pwy are in the model
          propwygene=oldpwygene
          propwygene[1,pwysgm]=ageneind
          propwygene[,pwysgm]=apply(as.matrix(propwygene[,pwysgm]),2,sort)
          proK=oldK
          loc=match(pwysgm,which(prophi==1))
          for (i in loc){
            x=genexp[propwygene[,which(prophi==1)[i]],]
            for(j in 1:n){
              for(l in j:n){
                proK[j,(l+(i-1)*n)]=proK[l,(j+(i-1)*n)]=exp(-sqrt(sum((x[,j]-x[,l])^2))/theta[which(prophi==1)[i]])
              }
            }
          }
          proKPHI=(proK)
          proSigma=diag(n)+aalpha*J+abeta*proKPHI%*%t(proKPHI)
          invproSigma=solve(proSigma)
          ###conditions cond2th and condI2th####
          subpwy3=which(colSums(propwygene!=0)>1)#pro subways that has more than 1 gene selected
          subpwy4=which(colSums(propwygene!=0)==1)#pro subways that has 1 gene selected
          
          gmno=colSums(as.matrix(propwygene[,subpwy3]!=0)) # subpwy3's gene #
          dffgene=matrix(0,nrow=max(gmno),ncol=length(subpwy3)) # genes can be selected >1 pwy 
          matchgene=integer(0)# the gene that doesn't fix condI2ga
          l=1
          for (i in 1:length(subpwy3)){
            if (sum(gmno[i]-gmno==1)>0){
              h=which(gmno[i]-gmno==1) # the number gene of ith pwy compare to the rest pwy, the location of the difference is 1.
              for(j in h){
                h1=match(propwygene[which(propwygene[,subpwy3[i]]!=0),subpwy3[i]],propwygene[which(propwygene[,subpwy3[j]]!=0),subpwy3[j]],nomatch=0) #the difference between genes in two pwy.
                if(sum(h1==0)==1){
                  matchgene[l]=propwygene[which(propwygene[,subpwy3[i]]!=0)[which(h1==0)],subpwy3[i]] #matchgene that can't be selected
                  l=l+1
                }
              }
            }
          }
          subpwy5lc=which(colSums(as.matrix(propwygene[,subpwy3]!=0))==2) #subpwy5lc is location of the subpwy3 has 2 genes
          if (length(subpwy4)>0){
            if (length(subpwy5lc)>0){
              for (i in subpwy5lc){
                h=propwygene[which(propwygene[,subpwy3[i]]!=0),subpwy3[i]]
                if (sum(match(h,propwygene[longestpwy,subpwy4],nomatch=0)!=0)==1){
                  h1=match(h,propwygene[longestpwy,subpwy4],nomatch=0)!=0
                  matchgene[l]=h[which(h1==F)]
                  l=l+1
                }
              }
            }
            
            matchgene=unique(c(matchgene,propwygene[longestpwy,subpwy4]))
          }
          
          if (length(matchgene)>0){
            for (i in 1:length(subpwy3)){
              h=propwygene[which(propwygene[,subpwy3[i]]!=0),subpwy3[i]]
              h1=h[which(match(h,matchgene,nomatch=0)==0)]
              if (length(h1)>0){dffgene[1:length(h1),i]=h1}
              
            }
          }else{dffgene=propwygene[(longestpwy-max(gmno)+1):longestpwy,subpwy3]}
          
          ###cond2ga and condI2ga
          condI2ga=colSums(as.matrix(dffgene)!=0) #genes that fit condI2ga in set cond2
          condI2=sum(condI2ga==0)
          condI2ga=condI2ga[which(condI2ga!=0)]
          
          logr1=-(n/2+a)*log(2*b+y%*%invproSigma%*%t(y))+sum(prophi)*log(pathprob)+(k-sum(prophi))*log(1-pathprob)+mu*sum(progm)+eta*t(progm)%*%R%*%progm
          logr11=logr1+log(length(subpwy2))+log(sum(1/condI2ga))
          logr11=logr11-log(length(subpwy3)-condI2)-log(sum(1/(rowSums(as.matrix(Snki[subpwy2,])!=0)-colSums(as.matrix(oldpwygene[,subpwy2])!=0))))
          logr=logr11-logr2
          if((log(runif(1)))<logr) {
            simphi[,g]=prophi
            simgm[,g]=progm
            logr2=logr1
            oldK=proKPHI
            oldpwygene=propwygene
            oldSigma=proSigma
            involdSigma=invproSigma
            break
          }else{
            simphi[,g]=simphi[,g-1]
            simgm[,g]=simgm[,g-1]
            break
          }
        }else{  #ar==2
          subpwy3=which(colSums(oldpwygene!=0)>1)#old subways that has more than 1 gene selected
          if (length(subpwy3)>0){
            subpwy4=which(colSums(oldpwygene!=0)==1)# old subways that has 1 gene selected
            
            gmno=colSums(as.matrix(oldpwygene[,subpwy3]!=0)) # subpwy3's gene #
            dffgene=matrix(0,nrow=max(gmno),ncol=length(subpwy3)) # genes can be selected >1 pwy 
            matchgene=integer(0)# the gene that doesn't fix condI2ga
            l=1
            for (i in 1:length(subpwy3)){
              if (sum(gmno[i]-gmno==1)>0){
                h=which(gmno[i]-gmno==1) # the number gene of ith pwy compare to the rest pwy, the location of the difference is 1.
                for(j in h){
                  h1=match(oldpwygene[which(oldpwygene[,subpwy3[i]]!=0),subpwy3[i]],oldpwygene[which(oldpwygene[,subpwy3[j]]!=0),subpwy3[j]],nomatch=0) #the difference between genes in two pwy.
                  if(sum(h1==0)==1){
                    matchgene[l]=oldpwygene[which(oldpwygene[,subpwy3[i]]!=0)[which(h1==0)],subpwy3[i]]
                    l=l+1
                  }
                }
              }
            }
            
            subpwy5lc=which(colSums(as.matrix(oldpwygene[,subpwy3]!=0))==2) #subpwy5lc is location of the subpwy3 has 2 genes
            if (length(subpwy4)>0){
              if (length(subpwy5lc)>0){
                for (i in subpwy5lc){
                  h=oldpwygene[which(oldpwygene[,subpwy3[i]]!=0),subpwy3[i]]
                  if (sum(match(h,oldpwygene[longestpwy,subpwy4],nomatch=0)!=0)==1){
                    h1=match(h,oldpwygene[longestpwy,subpwy4],nomatch=0)!=0
                    matchgene[l]=h[which(h1==F)]
                    l=l+1
                  }
                }
              }
              
              matchgene=unique(c(matchgene,oldpwygene[longestpwy,subpwy4]))
            }
            
            if (length(matchgene)>0){
              for (i in 1:length(subpwy3)){
                h=oldpwygene[which(oldpwygene[,subpwy3[i]]!=0),subpwy3[i]]
                h1=h[which(match(h,matchgene,nomatch=0)==0)]
                if (length(h1)>0){dffgene[1:length(h1),i]=h1}
                
              }
            }else{dffgene=oldpwygene[(longestpwy-max(gmno)+1):longestpwy,subpwy3]}
            
            
            
            ###cond2ga and condI2ga
            condI2ga=colSums(as.matrix(dffgene)!=0) #genes that fit condI2ga, which can be selected, given the pwy is seleced
            subpwy1=subpwy3[which(condI2ga!=0)] # pwy that fit the cond2 and condI2
            if(length(subpwy1)>0){
              if (length(subpwy1)>1){rpwyind=sample(subpwy1,1)}else{rpwyind=subpwy1} #randomly choose the pwy
              gminrpwy=dffgene[,which(subpwy3==rpwyind)] #genes in the pwy with 0
              subgenes=gminrpwy[which(gminrpwy!=0)] # genes that can be removed in the pwy
              if (length(subgenes)>1){rgeneind=sample(subgenes,1)}else{rgeneind=subgenes} #randomly select a gene to removed
              pwysgm=which(oldpwygene==rgeneind,arr.ind = T) # the pwy in the model that have the same gene
              prophi=simphi[,g-1]
              progm=simgm[,g-1]
              progm[rgeneind]=0
              propwygene=oldpwygene
              propwygene[pwysgm]=0
              propwygene[,pwysgm[,2]]=apply(as.matrix(propwygene[,pwysgm[,2]]),2,sort)
              proK=oldK
              loc=match(pwysgm[,2],which(prophi==1))
              for (i in loc){
                x=genexp[propwygene[,which(prophi==1)[i]],]
                for(j in 1:n){
                  for(l in j:n){
                    proK[j,(l+(i-1)*n)]=proK[l,(j+(i-1)*n)]=exp(-sqrt(sum((x[,j]-x[,l])^2))/theta[which(prophi==1)[i]])
                  }
                }
              }
              proKPHI=(proK)
              proSigma=diag(n)+aalpha*J+abeta*proKPHI%*%t(proKPHI)
              invproSigma=solve(proSigma)         
              logr1=-(n/2+a)*log(2*b+y%*%invproSigma%*%t(y))+sum(prophi)*log(pathprob)+(k-sum(prophi))*log(1-pathprob)+mu*sum(progm)+eta*t(progm)%*%R%*%progm
              ####conditions
              subpwy2=which(colSums(propwygene!=0)!=0)# pwy that are in the model
              subpwy2=subpwy2[which(colSums(propwygene[,subpwy2]!=0)<pwygmno[subpwy2])] #set G
              condI2ga=condI2ga[which(condI2ga!=0)] # delete the gene number equal to 0. 
              
              logr11=logr1+log(length(subpwy1))+ log(sum(1/(rowSums(Snki[subpwy2,]!=0)-colSums(propwygene[,subpwy2]!=0))))
              logr11=logr11-log(length(subpwy2))-log(sum(1/condI2ga))
              logr=logr11-logr2
              if((log(runif(1)))<logr) {
                simphi[,g]=prophi
                simgm[,g]=progm
                logr2=logr1
                oldK=proKPHI
                oldpwygene=propwygene
                oldSigma=proSigma
                involdSigma=invproSigma
                break
              }else{
                simphi[,g]=simphi[,g-1]
                simgm[,g]=simgm[,g-1]
                break
              }
            }
          }
        }
      } 
      
      if(move==3){
        if (sum(simphi[,g-1])>2) {ar=sample(1:2,1)}else{ar=1}
        if (ar==1){
          subpwy1=which(colSums(t(Snki)*simgm[,g-1])!=0) # all the pwy that have at least one gene in the model
          subpwy2=subpwy1[which(match(subpwy1,which(simphi[,g-1]!=0),nomatch=0)==0)] # the pwy are not in the model, but at least one of their genes are in the model
          if (length(subpwy2)>0){
            subpwygene=matrix(0,nrow=longestpwy,ncol=length(subpwy2))#pwy not in the model
            pwyoff=rep(1,length(subpwy2)) #indicator of the subpwy2 that can be selected
            for (i in 1:length(subpwy2)){ 
              h=which((t(Snki)*simgm[,g-1])[,subpwy2[i]]!=0,arr.ind = T)
              subpwygene[(longestpwy-length(h)+1):longestpwy,i]=h
              if (sum(colSums(as.matrix(subpwygene[,i]==oldpwygene))==longestpwy)==1){pwyoff[i]=0}
            }
            if (sum(pwyoff)>0){
              subpwy3=subpwy2[which(pwyoff==1)]
              if (length(subpwy3)>1){apwyind=sample(c(subpwy3),1)}else{apwyind=subpwy3}
              prophi=simphi[,g-1]
              progm=simgm[,g-1]
              prophi[apwyind]=1
              propwygene=oldpwygene
              propwygene[,apwyind]=subpwygene[,which(subpwy2==apwyind)]
              proK=oldK
              x=genexp[propwygene[,apwyind],]
              for(j in 1:n){
                for(l in j:n){
                  aK[j,l]=aK[l,j]=exp(-sqrt(sum((x[,j]-x[,l])^2))/theta[apwyind])
                }
              }
              proKPHI=matrix(0,nrow=n,ncol=sum(prophi)*n)
              aKloc=which(apwyind==which(prophi==1))
              if (aKloc==1) {proKPHI=cbind(aK,proK)}
              if (aKloc==sum(prophi)) {proKPHI=cbind(proK,aK)}
              if (aKloc>1 & aKloc<sum(prophi)){
                proKPHI[,1:((aKloc-1)*n)]=proK[,1:((aKloc-1)*n)]
                proKPHI[,((aKloc-1)*n+1):(aKloc*n)]=aK
                proKPHI[,(aKloc*n+1):(sum(prophi)*n)]=proK[,-(1:((aKloc-1)*n))]
              }
              
              proSigma=diag(n)+aalpha*J+abeta*proKPHI%*%t(proKPHI)
              invproSigma=solve(proSigma)
              logr1=-(n/2+a)*log(2*b+y%*%invproSigma%*%t(y))+sum(prophi)*log(pathprob)+(k-sum(prophi))*log(1-pathprob)+mu*sum(progm)+eta*t(progm)%*%R%*%progm
              ####cond3######
              cond3=0 # number of the pwy fits the cond3
              gene2=which(rowSums((t(Snki[which(prophi==1),])*progm)!=0)>1) # genes that have more than 1 time in the model.
              for (i in 1:sum(prophi)){
                h=propwygene[,which(prophi==1)[i]]
                if (sum(match(h[which(h!=0)],gene2,nomatch=0)!=0)==length(h[which(h!=0)])){cond3=cond3+1}
              }
              logr11=logr1+log(length(subpwy3))-log(cond3)
              logr=logr11-logr2
              if((log(runif(1)))<logr) {
                simphi[,g]=prophi
                simgm[,g]=progm
                logr2=logr1
                oldK=proKPHI
                oldpwygene=propwygene
                oldSigma=proSigma
                involdSigma=invproSigma
                break
              }else{
                simphi[,g]=simphi[,g-1]
                simgm[,g]=simgm[,g-1]
                break
              }
            }
          }
        }else{  #ar==2
          gene2=which(rowSums((t(Snki[which(simphi[,g-1]==1),])*simgm[,g-1])!=0)>1) # genes that have more than 1 time in the model.
          subpwy4=integer(0) # pwy that fit condI3
          j=1
          for (i in 1:sum(simphi[,g-1])){
            h=oldpwygene[,which(simphi[,g-1]==1)[i]]
            if (sum(match(h[which(h!=0)],gene2,nomatch=0)!=0)==length(h[which(h!=0)])){
              subpwy4[j]=which(simphi[,g-1]==1)[i]
              j=j+1
            }
          }
          if (length(subpwy4)>0){
            if (length(subpwy4)>1){rpwyind=sample(subpwy4,1)}else{rpwyind=subpwy4}
            prophi=simphi[,g-1]
            progm=simgm[,g-1]
            prophi[rpwyind]=0
            propwygene=oldpwygene
            propwygene[,rpwyind]=0
            rloc=which(rpwyind==which(simphi[,g-1]==1))
            proK=oldK[,-(((rloc-1)*n+1):(rloc*n))]
            proKPHI=proK
            proSigma=diag(n)+aalpha*J+abeta*proKPHI%*%t(proKPHI)
            invproSigma=solve(proSigma)
            logr1=-(n/2+a)*log(2*b+y%*%invproSigma%*%t(y))+sum(prophi)*log(pathprob)+(k-sum(prophi))*log(1-pathprob)+mu*sum(progm)+eta*t(progm)%*%R%*%progm
            #########condI3#####
            subpwy1=which(colSums(t(Snki)*progm)!=0) # all the pwy that have at least one gene in the model
            subpwy2=subpwy1[which(match(subpwy1,which(prophi!=0),nomatch=0)==0)] # the pwy are not in the model, but at least one of their genes are in the model
            if (length(subpwy2)>0){
              subpwygene=matrix(0,nrow=longestpwy,ncol=length(subpwy2))#pwy not in the model
              pwyoff=rep(1,length(subpwy2)) #indicator of the subpwy2 that can be selected
              for (i in 1:length(subpwy2)){ 
                h=which((t(Snki)*simgm[,g-1])[,subpwy2[i]]!=0,arr.ind = T)
                subpwygene[(longestpwy-length(h)+1):longestpwy,i]=h
                if (sum(colSums(as.matrix(subpwygene[,i]==propwygene))==longestpwy)==1){pwyoff[i]=0}
              }
            }
            logr11=logr1+log(length(subpwy4))-log(sum(pwyoff))
            logr=logr11-logr2
            if((log(runif(1)))<logr) {
              simphi[,g]=prophi
              simgm[,g]=progm
              logr2=logr1
              oldK=proKPHI
              oldpwygene=propwygene
              oldSigma=proSigma
              involdSigma=invproSigma
              break
            }else{
              simphi[,g]=simphi[,g-1]
              simgm[,g]=simgm[,g-1]
              break
            }
            
          }         
        }
      }
    } #repeat
    
    ###update theta########
    protheta=theta
    proK=oldK
    j=1
    for(i in which(simphi[,g]==1)){
      if (theta[i]-.5>0) {LB=theta[i]-.5}else{LB=0}
      if (theta[i]+.5<5) {UB=theta[i]+.5}else{UB=5}
      protheta[i]=runif(1,LB,UB)
      proK[,(((j-1)*n+1):(j*n))]=(oldK[,(((j-1)*n+1):(j*n))]^(theta[i]))^(1/protheta[i])
      proKPHI=(proK)
      proSigma=diag(n)+aalpha*J+abeta*proKPHI%*%t(proKPHI)
      invproSigma=solve(proSigma)
      logr1=-(n/2+a)*log(2*b+y%*%invproSigma%*%t(y))+sum(simphi[,g])*log(pathprob)+(k-sum(simphi[,g]))*log(1-pathprob)+mu*sum(simgm[,g])+eta*t(simgm[,g])%*%R%*%simgm[,g]
      logr=logr1-logr2
      if ((log(runif(1)))<logr) {
        theta[i]=protheta[i]
        logr2=logr1
        oldK=proKPHI
        oldSigma=proSigma
        involdSigma=invproSigma
      }
      j=j+1
    }
    simtheta[,g]=theta
    
    
    
    
    if (is.wholenumber(g/2000)){
      if (mean(colSums(simphi[,(g-1000):g]))>6 & mu>-5){
        mu=mu-.1
      }else if (mean(colSums(simphi[,(g-1000):g]))<2.5 & mu< -1.7){
        mu=mu+.1
      }
    }
    
    
    
    
    if (is.wholenumber(g/100)){
      simphio=simphi[,1:g]
      simgmo=simgm[,1:g]
      drt=paste(directory,"simoutput.RData",sep="")
      save(simphio,simgmo,simtheta, mu, eta, simy,  file=drt)
      
    }
    
  } # G
  
}








#data preparation
dat1=as.data.frame(read.table(file = "simY.txt"))
genexp=read.table(file = "simgene.txt")
Snki=read.table(file = "simSnki.txt")
Rnki=read.table(file = "simRnki.txt")
directory="simulation" #folder where the output will be
pathprob=.1 # percentage of the pathways expected to be important
geneprob=0.03 # percentage of the genes expected to be important
n=100 # number of subjects

G=1000 # iteration
abeta=.1 # prior for the Beta 
aalpha=10^2 # prior for the intercept, reasonably larger number
a=3 # augment data no use in simulation
b=.6 # augment data no use in simulation
mu=-2.8 # prior for the gene selection parameter
eta=.08 # prior for the gene selection parameter


NPKM(dat1,genexp,Snki,Rnk,pathprob,geneprob,n,G,abeta,aalpha,a,b,mu,eta,directory)
