library(RobustGaSP)
library(plot3D)
library(latex2exp)
library(Rcpp)
library(RcppEigen)
library(GpGp)
library(spNNGP)
library(laGP)
source('https://raw.githubusercontent.com/katzfuss-group/scaledVecchia/master/vecchia_scaled.R')
sourceCpp(file='src/functions.cpp')  ###
source('functions_lattice.R')




Y_full_mat=as.matrix(read.csv(file='satellite_data/t3_ascending_obs_matrix1_full.csv',header=F))
x1=(read.csv(file='satellite_data/t3_ascending_x_coordinate1_full.csv',header=F))[[1]]
x2=(read.csv(file='satellite_data/t3_ascending_y_coordinate1_full.csv',header=F))[[1]]


data_full = expand.grid(x1=x1,x2=x2)
data_full$value=as.vector(t(matrix(Y_full_mat/100,length(x2),length(x1))))

n1 = length(x1)
n2 = length(x2)
N=n1*n2

par(mgp=c(2,1,0), mar=c(3,3,4,2)+.1)#
zlim=range(data_full$value,na.rm=T)
image2D(matrix(data_full$value,n1,n2),x=x1,y=x2,xlab=expression(x[1]),ylab=expression(x[2]),
        clab='m/yr',main=paste('true'),zlim=zlim)


check_small_region = T
if(check_small_region){
  # select the size of the region 
  
  x1_sub_ind = c(101,300) ##200x200
  x2_sub_ind = c(101,300)
  

  
  sub_ind = which(data_full$x1 >= x1[x1_sub_ind[1]] & data_full$x1 <= x1[x1_sub_ind[2]] &
                    data_full$x2 >= x2[x2_sub_ind[1]] & data_full$x2 <= x2[x2_sub_ind[2]])
  
  
  x1 = x1[x1 >= x1[x1_sub_ind[1]] & x1 <= x1[x1_sub_ind[2]]]
  x2 = x2[x2 >= x2[x2_sub_ind[1]] & x2 <= x2[x2_sub_ind[2]]]
  
  data_full = data_full[sub_ind,]
  
  n1 = length(x1)
  n2 = length(x2)
  N=n1*n2
  
  #plot
  zlim=range(data_full$value,na.rm=T)
  image2D(matrix(data_full$value,n1,n2),x=x1,y=x2,xlab=expression(x[1]),ylab=expression(x[2]),
          clab='m/yr',main=paste('true'),zlim=zlim)
}



################ I add 0.25
settings = expand.grid(missing_prop=seq(0.05,0.25,0.05),missing_type = c("mask","random"))[,2:1]

NRMSE_rec = array(NA, c(7,2,dim(settings)[1]), 
                             dimnames = list(c("IKF","Vecchia1","Vecchia2","Svecchia1","Svecchia2","NNGP","laGP"),
                                             c("NRMSE_full","NRMSE_miss"),
                                             apply(settings,1,function(x) paste0(x[1],x[2])))) #vector(mode="list",dim(settings)[1])

time_rec = array(NA, c(7,2,dim(settings)[1]), 
                 dimnames = list(c("IKF","Vecchia1","Vecchia2","Svecchia1","Svecchia2","NNGP","laGP"),
                                 c("time_model","time_pred"),
                                 apply(settings,1,function(x) paste0(x[1],x[2])))) #vector(mode="list",dim(settings)[1])


pred_rec = array(NA, c(7,N,dim(settings)[1]), 
                 dimnames = list(c("IKF","Vecchia1","Vecchia2","Svecchia1","Svecchia2","NNGP","laGP"),
                                 NULL,
                                 apply(settings,1,function(x) paste0(x[1],x[2])))) #vector(mode="list",dim(settings)[1])

IKF_param_rec = matrix(NA,nr=dim(settings)[1],nc=3, dimnames = list(apply(settings,1,function(x) paste0(x[1],x[2])),c('beta1','beta2','tau')))

Y_obs_NA_rec = array(NA, c(dim(settings)[1],N), 
                     dimnames = list(apply(settings,1,function(x) paste0(x[1],x[2])),
                                     NULL)) #vector(mode="list",dim(settings)[1])

for(i_s in 1:dim(settings)[1]){
  print(settings[i_s,])
  
  missing_type = settings[i_s,1]
  missing_prop = settings[i_s,2]
  
  set.seed(i_s)
  ####################### select missing region ##############################
  full_ind = which(!is.na(data_full$value))
  
  if(missing_type == 'random'){
    N0=round(N*(1-missing_prop)) # assume only observe N0 samples
    obs_ind=sort(sample(full_ind,N0)) # need to be sorted
    
    miss_ind = setdiff(full_ind,obs_ind) # already sorted
  }else if(missing_type == 'mask'){
    #center = c(-3000,-2000)
    #center = c(.7*min(data_full$x1)+.3*max(data_full$x1),
    #           .7*min(data_full$x2)+.3*max(data_full$x2))
    center = c(.5*min(data_full$x1)+.5*max(data_full$x1),
               .5*min(data_full$x2)+.5*max(data_full$x2))
    
    dist_x1 = mean(diff(unique(data_full$x1)))
    dist_x2 = mean(diff(unique(data_full$x2)))
    
    radius = sqrt(missing_prop * N / pi *dist_x1*dist_x2)
    distances = apply(data_full, 1, function(x) sqrt(sum((x[1:2]-center)^2)))
    miss_ind = which(distances < radius) # need to be sorted
    obs_ind = setdiff(full_ind,miss_ind) # already sorted
    
    N0 = length(obs_ind)
  }
  non_obs_ind = setdiff(1:N,obs_ind)
  
  
  obs_input = as.matrix(data_full[obs_ind,1:2])
  obs_output = as.matrix(data_full[obs_ind,3])
  
  full_input = as.matrix(data_full[,1:2])
  full_output = as.matrix(data_full[,3])
  
  
  Y_obs_NA = full_output
  Y_obs_NA[-obs_ind] = NA
  Y_obs_NA_rec[i_s,] = Y_obs_NA

  #  par(mgp=c(2,1,0), mar=c(3,3,4,2)+.1)#
  #  image2D(matrix(full_output,n1,n2),x=x1,y=x2,xlab=expression(x[1]),ylab=expression(x[2]),
  #          clab='m/yr',main=paste('true'),zlim=zlim)
  #  image2D(matrix(Y_obs_NA,n1,n2),x=x1,y=x2,xlab=expression(x[1]),ylab=expression(x[2]),
  #          clab='m/yr',main=paste('obs'),zlim=zlim)
  #  title(TeX(paste0(missing_prop*100,'% ', missing_type,' missing')),line = 0.5,cex.main=0.8)#dev.off()

  
  ########################## IKF ##########################
  if(missing_type == 'random'){
    split_method = "random"
  }else if(missing_type == 'mask'){
    split_method = "disk"
  }
  
  kernel_type='matern_5_2'
  tilde_nu=0.1

  
  # first setting with param_ini=log(c(1e-4,1e-4,1))
  print("IKF")
  
  time_IKF=system.time({
  param_ini=log(c(1e-4,1e-4,1)) ##result look similar
  #param_ini=log(c(1e-4,1e-4,1e-2)) #initial value
  #param_ini=log(c(1e-4,1e-4,100)) ##
  zero_mean = TRUE
  
  res = satellite_est_param(param_ini,valid_prop=0.2,split_method=split_method,
                              N=N,N0=N0,full_input=full_input,full_output=full_output,x1=x1,x2=x2,obs_ind=obs_ind,kernel_type=kernel_type,zero_mean = zero_mean,
                              tilde_nu=tilde_nu, CGtol=0.001*sd(obs_output),CGmaxIte=100,maxit=100,print_par=F)
  if(zero_mean){
    param = res
    theta = 0
  }else{
    param = res$param
    theta = res$theta
  }
  
  })
  
  IKF_param_rec[i_s,] = param
  
  time_IKF_pred=system.time({
    pred_Z_CG=missing_prediction(param=param, Y_obs=obs_output, input1=x1, input2=x2, non_obs_ind=non_obs_ind, theta = theta, 
                                           kernel_type=kernel_type, tilde_nu=tilde_nu, tol = 1e-10, maxIte = 1000)
    
    
  })
  pred_rec[1,,i_s] = pred_Z_CG
  
  time_rec[1,1,i_s] = time_IKF[3]
  time_rec[1,2,i_s] = time_IKF_pred[3]

  #NRMSE of full obs
  NRMSE_rec[1,1,i_s]=sqrt(mean((pred_Z_CG[full_ind]-full_output[full_ind])^2))/sd(full_output[full_ind])
  #NRMSE of testing
  NRMSE_rec[1,2,i_s]=sqrt(mean((pred_Z_CG[miss_ind]-full_output[miss_ind])^2,na.rm=T))/sd(full_output[full_ind],na.rm = T) ###change this to  full_ind, otherwise not good for comparison 
  
  
  # par(mgp=c(2,1,0), mar=c(3,3,4,2)+.1)#
  # image2D(matrix(full_output,n1,n2),x=x1,y=x2,xlab=expression(x[1]),ylab=expression(x[2]),
  #         clab='m/yr',main=paste('true'),zlim=zlim)
  # image2D(matrix(pred_Z_CG,n1,n2),x=x1,y=x2,xlab=expression(x[1]),ylab=expression(x[2]),
  #         clab='m/yr',main=paste('obs'),zlim=zlim)
  # title(TeX(paste0(missing_prop*100,'% ', missing_type,' missing')),line = 0.5,cex.main=0.8)#dev.off()
  
  
  
  ########################## Vecchia ##########################
  # first setting with m_seq = c(10, 30)
  print("Vecchia1")
  
  time_Vecchia1=system.time({
    fit_Vecchia <- fit_model(y=obs_output, locs=obs_input,X=rep(1,N0), 
                             "matern25_isotropic",m_seq = c(10,30),silent=F)
  })
  
  time_Vecchia_pred1=system.time({
    pred_Vecchia=predictions(fit = fit_Vecchia,locs_pred=full_input,
                             X_pred=rep(1,N))
  })
  
  pred_rec[2,,i_s] = pred_Vecchia
  
  time_rec[2,1,i_s] = time_Vecchia1[3]
  time_rec[2,2,i_s] = time_Vecchia_pred1[3]
  
  #overall NRMSE
  NRMSE_rec[2,1,i_s]=sqrt(mean((pred_Vecchia[full_ind]-full_output[full_ind])^2))/sd(full_output[full_ind])
  #missing NRMSE
  NRMSE_rec[2,2,i_s]=sqrt(mean((pred_Vecchia[miss_ind]-full_output[miss_ind])^2,na.rm=T))/sd(full_output[full_ind],na.rm = T)
  
  # par(mgp=c(2,1,0), mar=c(3,3,4,2)+.1)#
  # image2D(matrix(full_output,n1,n2),x=x1,y=x2,xlab=expression(x[1]),ylab=expression(x[2]),
  #         clab='m/yr',main=paste('true'),zlim=zlim)
  # image2D(matrix(pred_Vecchia,n1,n2),x=x1,y=x2,xlab=expression(x[1]),ylab=expression(x[2]),
  #         clab='m/yr',main=paste('Vecchia'),zlim=zlim)
  # title(TeX(paste0(missing_prop*100,'% ', missing_type,' missing')),line = 0.5,cex.main=0.8)#dev.off()
  
  
  # second setting with m_seq = c(30,90)
  print("Vecchia2")
  
  time_Vecchia2=system.time({
    fit_Vecchia <- fit_model(y=obs_output, locs=obs_input,X=rep(1,N0), 
                             "matern25_isotropic",m_seq = c(30,90),silent=F) 
  })
  
  time_Vecchia_pred2=system.time({
    pred_Vecchia=predictions(fit = fit_Vecchia,locs_pred=full_input,
                             X_pred=rep(1,N),m = 150)
  })
  
  pred_rec[3,,i_s] = pred_Vecchia
  
  time_rec[3,1,i_s] = time_Vecchia2[3]
  time_rec[3,2,i_s] = time_Vecchia_pred2[3]
  
  #overall NRMSE
  NRMSE_rec[3,1,i_s]=sqrt(mean((pred_Vecchia[full_ind]-full_output[full_ind])^2))/sd(full_output[full_ind])
  #missing NRMSE
  NRMSE_rec[3,2,i_s]=sqrt(mean((pred_Vecchia[miss_ind]-full_output[miss_ind])^2,na.rm=T))/sd(full_output[full_ind],na.rm = T)
  
  
  ########################## scaled Vecchia ##########################
  # first setting with m_s = c(30)
  print("Svecchia1")
  
  time_SVecchia1=system.time({
    fit_scaled_vecchia=fit_scaled(obs_output,obs_input,ms=c(30),trend='zero', nu = 2.5,nug=NULL)
  })
  
  time_SVecchia_pred1=system.time({
    pred_scaled_vecchia=predictions_scaled(fit_scaled_vecchia,full_input,m=60)
  })
  pred_rec[4,,i_s] = pred_scaled_vecchia
  
  time_rec[4,1,i_s] = time_SVecchia1[3]
  time_rec[4,2,i_s] = time_SVecchia_pred1[3]
  
  
  #overall NRMSE
  NRMSE_rec[4,1,i_s]=sqrt(mean((pred_scaled_vecchia[full_ind]-full_output[full_ind])^2))/sd(full_output[full_ind])
  #missing NRMSE
  NRMSE_rec[4,2,i_s]=sqrt(mean((pred_scaled_vecchia[miss_ind]-full_output[miss_ind])^2,na.rm=T))/sd(full_output[full_ind])
  
  
  # second setting with m_s = c(90)
  print("Svecchia2")
  
  time_SVecchia2=system.time({
    fit_scaled_vecchia=fit_scaled(obs_output,obs_input,ms=c(90),trend='zero', nu = 2.5,nug=NULL)
  })
  
  time_SVecchia_pred2=system.time({
    pred_scaled_vecchia=predictions_scaled(fit_scaled_vecchia,full_input,m=150)
  })
  pred_rec[5,,i_s] = pred_scaled_vecchia
  
  time_rec[5,1,i_s] = time_SVecchia2[3]
  time_rec[5,2,i_s] = time_SVecchia_pred2[3]
  
  
  #overall NRMSE
  NRMSE_rec[5,1,i_s]=sqrt(mean((pred_scaled_vecchia[full_ind]-full_output[full_ind])^2))/sd(full_output[full_ind])
  #missing NRMSE
  NRMSE_rec[5,2,i_s]=sqrt(mean((pred_scaled_vecchia[miss_ind]-full_output[miss_ind])^2,na.rm=T))/sd(full_output[full_ind],na.rm = T)
  
  
  ########################## NNGP ##########################
  print("NNGP")
  cov.model <- "matern" #"matern", "exponential"
  
  sigma.sq <- 5e-5
  
  sigma.sq.IG <- c(2, sigma.sq)
  
  g <-  5
  theta.alpha <- as.matrix(expand.grid("phi"=seq(.1e-4, 0.1, length.out=g), 
                                       "alpha"=seq(1e-5/sigma.sq, 1e-4/sigma.sq, length.out=g),
                                       "nu"=seq(2.5, 2.5, by=1)))
  
  n.neighbors = 30
  
  time_NNGP=system.time({
    m.c <- spConjNNGP(obs_output~1, coords=obs_input, n.neighbors = n.neighbors,
                      k.fold = 5, score.rule = "crps",
                      n.omp.threads = 1,
                      theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG,
                      cov.model = cov.model,verbose=F)
  })
  
  
  ##prediction
  theta.alpha <- m.c$theta.alpha
  #names(theta.alpha) <- c("phi", "alpha", "nu")
  
  time_NNGP_pred=system.time({
    p_NNGP = spConjNNGP(obs_output~1, coords=obs_input, n.neighbors = n.neighbors,
                        X.0=matrix(1,N), coords.0=full_input,
                        n.omp.threads = 1,
                        theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG,
                        cov.model = cov.model,verbose=F)
  })
  pred_rec[6,,i_s] = p_NNGP$y.0.hat
  
  time_rec[6,1,i_s] = time_NNGP[3]
  time_rec[6,2,i_s] = time_NNGP_pred[3]
  
  #overall NRMSE
  NRMSE_rec[6,1,i_s]=sqrt(mean((p_NNGP$y.0.hat[full_ind]-full_output[full_ind])^2))/sd(full_output[full_ind])
  #missing NRMSE
  NRMSE_rec[6,2,i_s]=sqrt(mean((p_NNGP$y.0.hat[miss_ind]-full_output[miss_ind])^2,na.rm=T))/sd(full_output[full_ind],na.rm = T)
  
  
  
  ########################## laGP ##########################
  print("laGP")
  
  time_laGP=system.time({
    m_laGP_nn = aGP(X=obs_input, Z=obs_output, XX=full_input, method="nn", verb=0, end = 70)
  })
  pred_rec[7,,i_s] = m_laGP_nn$mean
  
  time_rec[7,1,i_s] = time_laGP[3]
  
  
  #overall NRMSE
  NRMSE_rec[7,1,i_s]=sqrt(mean((m_laGP_nn$mean[full_ind]-full_output[full_ind])^2))/sd(full_output[full_ind])
  #missing NRMSE
  NRMSE_rec[7,2,i_s]=sqrt(mean((m_laGP_nn$mean[miss_ind]-full_output[miss_ind])^2,na.rm=T))/sd(full_output[full_ind],na.rm = T)
  

}


#save.image("satellite_results.RData")


##### plot
load("satellite_results.RData")


#NRMSE
#NRMSE_array = abind(NRMSE_rec, along = 3)

NRMSE_df = melt(NRMSE_rec, varnames = c("Method", "Metric", "Scenario"), value.name = "NRMSE")

NRMSE_df = separate(NRMSE_df, Scenario, into = c("Type", "Prop"), sep = "(?<=[a-zA-Z])(?=[0-9])", convert = TRUE)



type = "mask" #c("mask","random")

pdf(paste0("plots/satellite_NRMSE_",type,".pdf"),width=4,height=1.8)
NRMSE_df %>%
  filter(Type == type & Metric == "NRMSE_miss") %>%
  filter(Method != "Vecchia1" & Method != "Svecchia1") %>%
  ggplot(aes(Prop,NRMSE,color = Method,linetype = Method)) + 
  geom_line()+
  geom_point(size = 1)+
  #scale_color_observable()+
  scale_color_manual(values = c("#4269D0FF", "#EFB118FF", "#FF725CFF", "#6CC5B0FF", "#3CA951FF"),
                     labels = c("IKF-CG","Vecchia","SVecchia","NNGP","laGP"))+
  scale_linetype_discrete(labels = c("IKF-CG","Vecchia","SVecchia","NNGP","laGP"))+
  #ylim(c(0.05,0.075))+
  xlab("Proportion")+
  theme_linedraw()+
  theme(legend.title = element_text(size = 10),legend.text=element_text(size=9),
        legend.key.size = unit(1, 'lines'),
        #axis.title.y=element_blank(),axis.text.y = element_text(size=10), #, angle=20
        axis.text=element_text(size=10),axis.title=element_text(size=10),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()


#Time
#time_array = abind(time_rec, along = 3)

time_mat = apply(time_rec,c(1,3),function(x) sum(x,na.rm=T))

time_df = melt(time_mat, varnames = c("Method", "Scenario"), value.name = "Time")

time_df = separate(time_df, Scenario, into = c("Type", "Prop"), sep = "(?<=[a-zA-Z])(?=[0-9])", convert = TRUE)

pdf(paste0("plots/satellite_time_",type,".pdf"),width=4,height=1.8)
time_df %>%
  filter(Type == type) %>%
  filter(Method != "Vecchia1" & Method != "Svecchia1") %>%
  ggplot(aes(Prop,Time,color = Method,linetype = Method)) + 
  geom_line()+
  geom_point(size = 1)+
  #scale_color_observable()+
  scale_color_manual(values = c("#4269D0FF", "#EFB118FF", "#FF725CFF", "#6CC5B0FF", "#3CA951FF"),
                     labels = c("IKF-CG","Vecchia","SVecchia","NNGP","laGP"))+
  scale_linetype_discrete(labels = c("IKF-CG","Vecchia","SVecchia","NNGP","laGP"))+
  #ylim(c(0.05,0.075))+
  xlab("Proportion")+
  ylab("Time (s)")+
  theme_linedraw()+
  theme(legend.title = element_text(size = 10),legend.text=element_text(size=9),
        legend.key.size = unit(1, 'lines'),
        #axis.title.y=element_blank(),axis.text.y = element_text(size=10), #, angle=20
        axis.text=element_text(size=10),axis.title=element_text(size=10),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()


i_s = 10


dat = data.frame(x1 = data_full$x1, x2 = data_full$x2, truth = data_full$value,
                 obs = Y_obs_NA_rec[i_s,], IKF_pred = pred_rec[1,,i_s])
#full_ind=which(!is.na(dat$truth))
dat$IKF_pred[-full_ind] = NA

zlim = range(dat$truth,na.rm = T)

pdf(paste0("plots/satellite_obs",i_s,".pdf"),width=4,height=2.5)
ggplot(dat) + geom_tile(aes(x1,x2,fill=obs)) +
  theme_bw() +
  scale_fill_distiller(palette="RdYlBu",name="",limits=zlim) + 
  ggtitle("Observation")+
  xlab(expression(x[1])) +
  ylab(expression(x[2]))+
  theme(axis.text=element_text(size=11),axis.title=element_text(size=12),
        plot.title = element_text(hjust = 0.5,size=13),
        legend.text=element_text(size=10),
        #panel.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()

pdf(paste0("plots/satellite_pred",i_s,".pdf"),width=4,height=2.5)
ggplot(dat) + geom_tile(aes(x1,x2,fill=IKF_pred)) +
  theme_bw() +
  scale_fill_distiller(palette="RdYlBu",name="",limits=zlim) + 
  ggtitle("IKF-CG Prediction")+
  xlab(expression(x[1])) +
  ylab(expression(x[2]))+
  theme(axis.text=element_text(size=11),axis.title=element_text(size=12),
        plot.title = element_text(hjust = 0.5,size=13),
        legend.text=element_text(size=10),
        #panel.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()

pdf("plots/satellite_truth.pdf",width=4,height=2.5)
ggplot(dat) + geom_tile(aes(x1,x2,fill=truth)) +
  theme_bw() +
  scale_fill_distiller(palette="RdYlBu",name="",limits=zlim) + 
  ggtitle("Truth")+
  xlab(expression(x[1])) +
  ylab(expression(x[2]))+
  theme(axis.text=element_text(size=11),axis.title=element_text(size=12),
        plot.title = element_text(hjust = 0.5,size=13),
        legend.text=element_text(size=10),
        #panel.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()





