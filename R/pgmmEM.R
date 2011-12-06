.packageName<-'pgmm'

init_load<-function(x4,z4,G4,p4,q4){
    sampcov<-list(); 
    for(g in 1:G4){sampcov[[g]]<-matrix(0,p4,p4);}
    mu<-matrix(0,G4,p4)
    n<-colSums(z4)
    pi<-colSums(z4)/dim(x4)[1]
    for(g in 1:G4){mu[g,]<-apply(x4*z4[,g],2,sum)/(n[g]);}
    for(g in 1:G4){sampcov[[g]]<-cov.wt(x4,wt=(z4[,g]),center=mu[g,])$cov*(sum(z4[,g]^2)-1)/sum(z4[,g]);}
    lambda<-list();
    lambda1_tmp<-c(rep(0,p4*q4*G4));
    s<-1;
    for(g in 1:G4){
        evec<-eigen(sampcov[[g]])$vector
        eval<-eigen(sampcov[[g]])$values
        for(i in 1:p4){
            for(j in 1:q4){
                lambda1_tmp[s]<-sqrt(eval[j])*evec[i];
                s<-s+1;
            }
        }
    }
    lambda[["sep"]]<-lambda1_tmp
    
    psi<-list()
    lam_mat<-list()
    k4<-0
    for(g in 1:G4){
        lam_mat[[g]]<-matrix(lambda1_tmp[(1+k4):(g*p4*q4)],nrow=p4,ncol=q4,byrow=TRUE)
        k4<-p4*q4*g
    }
    temp_p6<-c(rep(0,p4))
    for(g in 1:G4){
        temp_p6 <- temp_p6 + pi[g]*abs(diag(sampcov[[g]]-lam_mat[[g]]%*%t(lam_mat[[g]])))
    }
    psi[[6]]<-temp_p6
    psi[[5]]<-sum(psi[[6]])/p4
    psi_tmp<-matrix(0,G4,p4)
    for (g in 1:G4){
        psi_tmp[g,]<-abs(diag(sampcov[[g]]-lam_mat[[g]]%*%t(lam_mat[[g]])))
    } 
    psi[[7]]<-rowMeans(psi_tmp)
    psi[[8]]<-as.vector(t(psi_tmp))
    
    stilde<-matrix(0,p4,p4)
    for(g in 1:G4){
        stilde<-stilde+pi[g]*sampcov[[g]];
    }
    evec<-eigen(stilde)$vector
    eval<-eigen(stilde)$values
    lambda_tilde<-matrix(0,p4,q4)
    s<-1
    for(i in 1:p4){
        for(j in 1:q4){
            lambda1_tmp[s]<-sqrt(eval[j])*evec[i];
            s<-s+1
        }
    }
    lambda[["tilde"]]<-lambda1_tmp
    
    lam_mat[[1]]<-matrix(lambda1_tmp[1:(p4*q4)],nrow=p4,ncol=q4,byrow=TRUE)
    psi[[2]]<-abs(diag(stilde-lam_mat[[1]]%*%t(lam_mat[[1]])))
    psi[[1]]<-sum(psi[[2]])/p4
    psi_tmp<-matrix(0,G4,p4)
    for (g in 1:G4){
        psi_tmp[g,]<-abs(diag(sampcov[[g]]-lam_mat[[1]]%*%t(lam_mat[[1]])))
    } 
    psi[[3]]<-rowMeans(psi_tmp)
    psi[[4]]<-as.vector(t(psi_tmp))
    psi[[9]]<-psi[[3]]
    psi[[10]]<-psi[[7]]
    psi[[11]]<-psi[[1]]
    psi[[12]]<-psi[[5]]
    
    lambda[["psi"]]<-psi
    lambda
}
endPrint<-function(icl,zstart,loop,m_best,q_best,G_best,bic_best,class_ind){
    start_names<-c("NA","k-means","custom")
    if(!class_ind){
        if(!icl){
            if(zstart==1){
                if(loop==1){
                    cat("Based on 1 random start, the best model (BIC) for the range of factors and components used is a ", m_best, " model with q = ", q_best, 
                        " and G = ", G_best, ".\nThe BIC for this model is ", bic_best, ".", sep="")
                }else{
                    cat("Based on ", loop, " random starts, the best model (BIC) for the range of factors and components used is a ", m_best, " model with q = ", q_best, " and G = ", G_best, ".\nThe BIC for this model is ", bic_best, ".", sep="")
                }
            }else{
                cat("Based on ", start_names[zstart], " starting values, the best model (BIC) for the range of factors and components used is a ", m_best, " model with q = ", q_best, " and G = ", G_best, ".\nThe BIC for this model is ", bic_best, ".", sep="")
            }
        }else{
            if(zstart==1){
                if(loop==1){
                    cat("Based on 1 random start, the best model (ICL) for the range of factors and components used is a ", m_best, " model with q = ", q_best, 
                        " and G = ", G_best, ".\nThe ICL for this model is ", bic_best, ".", sep="")
                }else{
                    cat("Based on ", loop, " random starts, the best model (ICL) for the range of factors and components used is a ", m_best, " model with q = ", q_best, " and G = ", G_best, ".\nThe ICL for this model is ", bic_best, ".", sep="")
                }
            }else{
                cat("Based on ", start_names[zstart], " starting values, the best model (ICL) for the range of factors and components used is a ", m_best, " model with q = ", q_best, " and G = ", G_best, ".\nThe ICL for this model is ", bic_best, ".", sep="")
            }
        }
    }else{
        if(!icl){
            cat("Based on the labelled and unlabelled data provided, the best model (BIC) for the range of factors and components used is a ", m_best, " model with q = ", q_best, " and G = ", G_best, ".\nThe BIC for this model is ", bic_best, ".", sep="")
        }else{
            cat("Based on the labelled and unlabelled data provided, the best model (ICL) for the range of factors and components used is a ", m_best, " model with q = ", q_best, " and G = ", G_best, ".\nThe ICL for this model is ", bic_best, ".", sep="")
        }
    }        
}
pgmmEM<-function(x,class=NULL,icl=FALSE,zstart=2,cccStart=TRUE,loop=3,zlist=NULL,qmax=2,qmin=1,Gmax=2,Gmin=1,modelSubset=NULL,seed=123456,tol=0.1,relax=FALSE){
	set.seed(seed)
	x<-as.matrix(x)
	run_pgmm<-function(x,z,bic,cls,q,p,G,N,model,cluster,lambda,psi,TOL=tol){
		p4<-.C("pgmm_c",as.double(x),as.double(z),as.double(bic),as.integer(cls),as.integer(q),as.integer(p),as.integer(G),as.integer(N),
			   as.integer(model),as.integer(cluster),as.double(lambda),as.double(psi),as.double(TOL),PACKAGE="pgmm")
		list(p4[[2]],p4[[3]])
	}
	is.int<-function(no){
		abs(no-round(no))<1e-15
	}
	if((tol>0.1)||(!is.double(tol))){stop("Invalid entry for tol.")}
	models_all<-c("CCC","CCU","CUC","CUU","UCC","UCU","UUC","UUU","CCUU","UCUU","CUCU","UUCU")
	model_num<-1:12
	names(model_num)<-models_all
	if(is.null(modelSubset)){modelSubset<-models_all}
	bic_out<-list()
	if(!is.null(zlist)||zstart==3){if(!is.list(zlist)){stop("Expected a list for zlist.")}}
	if(!is.null(class)){if(any(!is.double(class))&any(!is.int(class))){stop("The vector class may contain integers only.")}}
	G_offset<-Gmin-1
	q_offset<-qmin-1
	N<-dim(x)[1]
	p<- dim(x)[2]
	x1<-as.vector(t(x))
	bic_max<--Inf
	bic_best<--Inf
	if(!is.null(zlist)){
		for(g1 in Gmin:Gmax){
			if(g1>1){
				if((any(!is.integer(zlist[[g1]])))){stop("Each element of zlist (G>1) must contain integers only.");}
				if(length(zlist[[g1]])!=N){stop("Each element of zlist (G>1) must have length equal to the number of samples.");}
			}
		}
	}
	test_paras<-c(qmax,qmin,Gmax,Gmin)
	if(!is.null(class))if(length(class)!=N){stop("The vector class must have length equal to the number of samples.");}
	if(any(!is.double(test_paras))&any(!is.integer(test_paras))){stop("The parameters qmax, qmin, Gmax, Gmin take integer values only.")}
	if(Gmax<Gmin){stop("Gmax<Gmin");}
	if(qmax<qmin){stop("qmax<qmin");}
	if(!relax){
        if(qmax>p/2){stop("qmax>p/2 is not allowed; set relax=TRUE to relax this constraint.");}
    }else{
        if(qmax>(3*p/4)){stop("qmax>3p/4 is not allowed.");}
    }
	if(!is.double(x)){stop("All elements of x must have type double.");}
	if(!is.logical(cccStart)){stop("cccStart takes logical values only.");}
	if(is.null(class)){class<-rep(0,N);class_ind<-0;}else{class_ind<-1;}
	for(mod in modelSubset){
		bic_temp<-matrix(NA,Gmax-Gmin+1,qmax-qmin+1)
		rownames(bic_temp)<-c(Gmin:Gmax)
		colnames(bic_temp)<-c(qmin:qmax)	
		bic_out[[mod]]<-bic_temp
	}
	if(class_ind){
        if(max(class)>Gmin){stop("The parameter Gmin cannot be less than max(class).")}
		for(g1 in Gmin:Gmax){
			zt<-matrix(0,N,g1)
			cls_ind<-(class==0)
			for (i in 1:N){
				if(cls_ind[i]){zt[i,]<-1/g1}
				else{zt[i,class[i]]<-1}	
			}
			z1<-as.vector(t(zt))
			LAMBDA<-list()
			for (q1 in qmin:qmax){
				LAMBDA[[q1-q_offset]]<-init_load(x,zt,g1,p,q1)
			}		
			for (m in modelSubset){
				for (q1 in qmin:qmax){
                    if(substr(m,1,1)=="C"){
                        lam_temp<-LAMBDA[[q1-q_offset]]$tilde
                    }else{
                        lam_temp<-LAMBDA[[q1-q_offset]]$sep
                    }
                    psi_temp<-LAMBDA[[q1-q_offset]]$psi[[model_num[m]]]
                    temp<-run_pgmm(x1,z1,0,class,q1,p,g1,N,model_num[m],class_ind,lam_temp,psi_temp)
					bic_out[[m]][g1-G_offset,q1-q_offset]<-temp[[2]]
					if(!is.nan(temp[[2]])){
						if(icl&(g1>1)){
                            z_mat_tmp<-matrix(temp[[1]],nrow=N,ncol=g1,byrow=TRUE)
                            mapZ<-rep(0,N)
                            for(i9 in 1:N){
                                mapZ[i9]<-which(z_mat_tmp[i9,1:g1]==max(z_mat_tmp[i9,1:g1]))
                            }
                            icl2<-0
                            for(i9 in 1:N){
                                icl2<-icl2+log(z_mat_tmp[i9,mapZ[i9]])
                            }
                            bic_out[[m]][g1-G_offset,q1-q_offset]<-bic_out[[m]][g1-G_offset,q1-q_offset]+icl2
                        }
                        if(temp[[2]]>bic_max){
							z_best<-temp[[1]];bic_best<-temp[[2]];
							bic_max<-bic_best;G_best<-g1;q_best<-q1;
							m_best<-m;
						}
					}
				}
			}
		}
	}else{
		bic_start<-matrix(NA,Gmax-Gmin+1,qmax-qmin+1)
		if(zstart==1){
			for(l in 1:loop){
				for(g1 in Gmin:Gmax){
					z<-matrix(0,N,g1)
                    for(i in 1:N){
                        sum<-0
                        for (j in 1:g1){
                            z[i,j]<-runif(1,0,1);
                            sum<-sum+z[i,j];
                        }
                        for (j in 1:g1){
                            z[i,j]<-z[i,j]/sum
                        }
                    }
					z1<-as.vector(t(z))
					LAMBDA<-list()
					for (q1 in qmin:qmax){
						LAMBDA[[q1-q_offset]]<-init_load(x,z,g1,p,q1)
					}
					if(cccStart){
                        bic_ccc_max<--Inf
						for (q1 in qmin:qmax){
							temp<-run_pgmm(x1,z1,0,class,q1,p,g1,N,1,class_ind,LAMBDA[[q1-q_offset]]$tilde,LAMBDA[[q1-q_offset]]$psi[[1]])
							bic_start[g1-G_offset,q1-q_offset]<-temp[[2]]
							if(!is.nan(temp[[2]])){
								if(icl&(g1>1)){
									z_mat_tmp<-matrix(temp[[1]],nrow=N,ncol=g1,byrow=TRUE)
									mapZ<-rep(0,N)
									for(i9 in 1:N){
										mapZ[i9]<-which(z_mat_tmp[i9,1:g1]==max(z_mat_tmp[i9,1:g1]))
									}
									icl1<-0
									if(icl){
                                        for(i9 in 1:N){
                                            icl1<-icl1+log(z_mat_tmp[i9,mapZ[i9]])
                                        }
                                    }
 									bic_start[g1-G_offset,q1-q_offset]<-bic_start[g1-G_offset,q1-q_offset]+icl1
								}
								if(bic_start[g1-G_offset,q1-q_offset]>bic_ccc_max){
									z_init_best<-temp[[1]];
                                    bic_ccc_max<-bic_start[g1-G_offset,q1-q_offset]
								}
							}
						}
						z_init_mat<-matrix(z_init_best,nrow=N,ncol=g1,byrow=TRUE);						
						for (q1 in qmin:qmax){
							LAMBDA[[q1-q_offset]]<-init_load(x,z_init_mat,g1,p,q1)
						}
						z1<-as.vector(t(z_init_mat))
					}
                    for (m in modelSubset){
                        for (q1 in qmin:qmax){
                            if(substr(m,1,1)=="C"){
                                lam_temp<-LAMBDA[[q1-q_offset]]$tilde
                            }else{
                                lam_temp<-LAMBDA[[q1-q_offset]]$sep
                            }
                            psi_temp<-LAMBDA[[q1-q_offset]]$psi[[model_num[m]]]
                            temp<-run_pgmm(x1,z1,0,class,q1,p,g1,N,model_num[m],class_ind,lam_temp,psi_temp)
                            bic_out[[m]][g1-G_offset,q1-q_offset]<-temp[[2]]
                            if(!is.nan(temp[[2]])){
                                if(icl&(g1>1)){
                                    z_mat_tmp<-matrix(temp[[1]],nrow=N,ncol=g1,byrow=TRUE)
                                    mapZ<-rep(0,N)
                                    for(i9 in 1:N){
                                        mapZ[i9]<-which(z_mat_tmp[i9,1:g1]==max(z_mat_tmp[i9,1:g1]))
                                    }
                                    icl2<-0
                                    for(i9 in 1:N){
                                        icl2<-icl2+log(z_mat_tmp[i9,mapZ[i9]])
                                    }
                                    bic_out[[m]][g1-G_offset,q1-q_offset]<-bic_out[[m]][g1-G_offset,q1-q_offset]+icl2
                                }
                                if(temp[[2]]>bic_max){
                                    z_best<-temp[[1]];bic_best<-temp[[2]];
                                    bic_max<-bic_best;G_best<-g1;q_best<-q1;
                                    m_best<-m;
                                }
                            }
                        }
                    }
				}
			}
		}else if((zstart==2)||(zstart==3)){
			for(g1 in Gmin:Gmax){
				z<-matrix(0,N,g1)
				if(zstart==3){
					if(g1==1){z_ind<-c(rep(1,N))}
					else {z_ind<-zlist[[g1]]}
				}
                if(zstart==2){
					if(g1==1){z_ind<-c(rep(1,N))}
					else{set.seed(123456);z_ind<-kmeans(x,g1,nstart=5)$cluster;}
				}
				for (i in 1:N){
					z[i,z_ind[i]]<-1	
				}
				LAMBDA<-list()
				for (q1 in qmin:qmax){
					LAMBDA[[q1-q_offset]]<-init_load(x,z,g1,p,q1)
				}
				z1<-as.vector(t(z))
				for (m in modelSubset){
					for (q1 in qmin:qmax){
						if(substr(m,1,1)=="C"){
							lam_temp<-LAMBDA[[q1-q_offset]]$tilde
						}else{
							lam_temp<-LAMBDA[[q1-q_offset]]$sep
						}
						psi_temp<-LAMBDA[[q1-q_offset]]$psi[[model_num[m]]]
                        temp<-run_pgmm(x1,z1,0,class,q1,p,g1,N,model_num[m],class_ind,lam_temp,psi_temp)
						bic_out[[m]][g1-G_offset,q1-q_offset]<-temp[[2]]
						if(!is.nan(temp[[2]])){
							if(icl&(g1>1)){
								z_mat_tmp<-matrix(temp[[1]],nrow=N,ncol=g1,byrow=TRUE)
								mapZ<-rep(0,N)
								for(i9 in 1:N){
									mapZ[i9]<-which(z_mat_tmp[i9,1:g1]==max(z_mat_tmp[i9,1:g1]))
								}
								icl2<-0
								for(i9 in 1:N){
									icl2<-icl2+log(z_mat_tmp[i9,mapZ[i9]])
								}
								bic_out[[m]][g1-G_offset,q1-q_offset]<-bic_out[[m]][g1-G_offset,q1-q_offset]+icl2
							}
							if(temp[[2]]>bic_max){
								z_best<-temp[[1]];bic_best<-bic_out[[m]][g1-G_offset,q1-q_offset];
								bic_max<-bic_best;G_best<-g1;q_best<-q1;
								m_best<-m;
							}
						}
					}
				}
			}
		}else{
            stop("Invalid entry for zstart: 1 random; 2 k-means; 3 user-specified list.")
        }
	}
	z_mat<-matrix(z_best,nrow=N,ncol=G_best,byrow=TRUE)
	class_best<-rep(0,N)
	for(i in 1:N){
		class_best[i]<-which(z_mat[i,1:G_best]==max(z_mat[i,1:G_best]))
	}
    endPrint(icl,zstart,loop,m_best,q_best,G_best,bic_best,class_ind)
	if(!icl){
		foo<-list(map=class_best,model=m_best,g=G_best,q=q_best,bic=bic_out,z_hat=z_mat,plot_info=list(Gmin,Gmax,modelSubset,icl),summ_info=list(icl,zstart,loop,bic_best,class_ind))
	}else{
		foo<-list(map=class_best,model=m_best,g=G_best,q=q_best,icl=bic_out,z_hat=z_mat,plot_info=list(Gmin,Gmax,modelSubset,icl),summ_info=list(icl,zstart,loop,bic_best,class_ind))
	}
    class(foo)<-"pgmm"
    foo
}
summary.pgmm<-function(object,...){
    a<-object$summ_info[[1]]
    b<-object$summ_info[[2]]
    c<-object$summ_info[[3]]
    d<-object$model
    e<-object$q
    f<-object$g
    k<-object$summ_info[[4]]
    m<-object$summ_info[[5]]
    endPrint(a,b,c,d,e,f,k,m)
}
print.pgmm<-function(x,...){
    bicl<-x$plot_info[[4]]
    if(!bicl){
        cat("BIC for each model, number of components (rows), and number of factors (columns).\n")
        print.default(x$bic)
    }else{
        cat("ICL for each model, number of components (rows), and number of factors (columns).\n")
        print.default(x$icl)
    }
}
plot.pgmm<-function(x,onlyAll=FALSE,...){
    x$plot_info[[3]]->models
    x$plot_info[[4]]->icl1
	if(length(models)<3){
		par(mfrow=c(1,2),ask=FALSE)
	}else if(length(models)<5){
		par(mfrow=c(2,2),ask=FALSE)
	}else if(length(models)<7){
		par(mfrow=c(2,3),ask=FALSE)
	}else if(length(models)<10){
		par(mfrow=c(3,3),ask=FALSE)
	}else{
		par(mfrow=c(3,4),ask=FALSE)
	}
	if(icl1){
        bicl<-x$icl
        ylabel<-"ICL"
    }else{
        bicl<-x$bic
        ylabel<-"BIC"
    }
    for(k in models){
        n<-which(models==k)
        matplot(c(x$plot_info[[1]]:x$plot_info[[2]]),bicl[[n]],type="b",ylab=ylabel,xlab="G",main=k,xaxt="n")
        axis(1,at=c(x$plot_info[[1]]:x$plot_info[[2]]))
	}
    
    if(!onlyAll){
		for(k in models){
	        par(mfrow=c(1,1),ask=TRUE)
	        n<-which(models==k)
    	    matplot(c(x$plot_info[[1]]:x$plot_info[[2]]),bicl[[n]],type="b",ylab=ylabel,xlab="G",main=k,xaxt="n")
        	axis(1,at=c(x$plot_info[[1]]:x$plot_info[[2]]))
        }
	}
}