library(RobustGaSP)
library(plot3D)
library(Rcpp)
library(RcppEigen)
library(GpGp)
library(spNNGP)
library(laGP)
library(GPvecchia)
source('https://raw.githubusercontent.com/katzfuss-group/scaledVecchia/master/vecchia_scaled.R')
sourceCpp(file='src/functions.cpp')  ###
source('functions_lattice.R')

n_repeat = 20

settings = data.frame(missing_type = c("mask","mask","random"), 
                      missing_prop = c(0.2, 0.2, 0.2))

time_rec = NRMSE_rec = IKF_param_rec = obs_IKF_pred_rec = vector(mode="list",dim(settings)[1])

names(time_rec) = names(NRMSE_rec) = names(IKF_param_rec) = names(obs_IKF_pred_rec) = paste0('setting',1:dim(settings)[1])


for(i_s in 1:dim(settings)[1]){
  print(settings[i_s,])
  missing_type = settings[i_s,1]
  missing_prop = settings[i_s,2]
  
  NRMSE_rec[[i_s]] = array(NA,c(7,2,n_repeat),dimnames = list(c("IKF","Vecchia1","Vecchia2","Svecchia1","Svecchia2","NNGP","laGP"),
                                                              c("NRMSE_full","NRMSE_miss"),paste0("rep",1:n_repeat)))
  
  time_rec[[i_s]] = array(NA,c(7,2,n_repeat),dimnames = list(c("IKF","Vecchia1","Vecchia2","Svecchia1","Svecchia2","NNGP","laGP"),
                                                             c("time_model","time_pred"),paste0("rep",1:n_repeat)))
  
  IKF_param_rec[[i_s]] = array(NA,c(n_repeat,3),dimnames = list(paste0("rep",1:n_repeat),c("beta1","beta2","tau")))
  
  
  
  for(i_repeat in 1:n_repeat){
    print(i_repeat)
    
    set.seed(i_repeat+(i_s-1)*n_repeat)
    
    
    ########################## simulate data ##########################
    ## data size
    n1 = 100
    n2 = 100
    N = n1*n2
    
    sigma0=10
    
    ## equally space 
    input1=seq(-5,10,length.out=n1) ##row
    input2=seq(0,15,length.out=n2) ##column
    
    input=as.matrix(expand.grid(input1=input1,input2=input2))
    
    
    Z_full =  apply(input, 1, branin)
    Y_full = Z_full + sigma0*rnorm(N)
    
    
    if(missing_type == 'random'){
      N0=round(N*(1-missing_prop)) # assume only observe N0 samples
      obs_ind=sort(sample(1:N,N0)) # need to be sorted
      
      inside_ind = setdiff(1:N,obs_ind) # already sorted
    }else if(missing_type == 'mask'){
      #disk_center=c(5,6)
      
      if(i_s == 1){
        disk_center = c(.3*min(input1)+.7*max(input1),
                        .5*min(input2)+.5*max(input2))
      }else if (i_s == 2){
        disk_center = c(.7*min(input1)+.3*max(input1),
                        .7*min(input2)+.3*max(input2))
      }
      
      
      dist_x1 = mean(diff(input1))
      dist_x2 = mean(diff(input2))
      
      radius = sqrt(missing_prop * N / pi *dist_x1*dist_x2)
      distances = apply(input, 1, function(x) sqrt(sum((x-disk_center)^2)))
      inside_ind = which(distances < radius) # need to be sorted
      obs_ind = setdiff(1:N,inside_ind) # already sorted
      
      N0 = length(obs_ind)
    }
    
    
    Y_obs_NA=Y_full
    Y_obs_NA[inside_ind] = NA
    
    
    Y_obs = Y_full[-inside_ind]
    N0 = N-length(inside_ind)
    
    if(i_repeat == 1){
      obs_IKF_pred_rec[[i_s]] = array(NA,c(2,N),dimnames = list(c("obs","pred"),NULL))
      obs_IKF_pred_rec[[i_s]][1,] = Y_obs_NA
    }
    
    
    ########################## IKF ##########################
    print("IKF")
    
    kernel_type='matern_5_2'
    tilde_nu=0.1
    
    time_IKF=system.time({

      param_ini=log(c(1,1,1))
      
      param = branin_est_param(Y_obs_NA=Y_obs_NA, obs_ind=obs_ind, inside_ind=inside_ind, N=N, N0=N0, valid_prop = 0.2, 
                               input1=input1, input2=input2, kernel_type=kernel_type, param_ini=param_ini,
                               tilde_nu=tilde_nu, CGtol = 0.001*sd(Y_obs), CGmaxIte=100,maxit=100,print_par=F)
    })
    
    IKF_param_rec[[i_s]][i_repeat,] = param
    
    time_IKF_pred=system.time({
      pred_Z_CG=missing_prediction(param=param, Y_obs=Y_obs, input1=input1, input2=input2, non_obs_ind=inside_ind, theta = 0,
                                   kernel_type=kernel_type, tilde_nu=tilde_nu, tol = 1e-1, maxIte = 1000)
      
      
    })
    
    if(i_repeat == 1){
      obs_IKF_pred_rec[[i_s]][2,] = pred_Z_CG
    }
    
    time_rec[[i_s]][1,1,i_repeat] = time_IKF[3]
    time_rec[[i_s]][1,2,i_repeat] = time_IKF_pred[3]
    
    #overall NRMSE
    NRMSE_rec[[i_s]][1,1,i_repeat]=sqrt(mean((pred_Z_CG-Z_full)^2))/sd(Z_full)
    
    #missing NRMSE
    NRMSE_rec[[i_s]][1,2,i_repeat]=sqrt(mean((pred_Z_CG[inside_ind]-Z_full[inside_ind])^2))/sd(Z_full)
    
    # image2D(matrix(Z_full,n1,n2),x=input1,y=input2,main='true mean',xlab='x1',ylab='x2',zlim=range(Y_full))
    # image2D(matrix(pred_Z_CG,n1,n2),x=input1,y=input2,main='predicted mean',xlab='x1',ylab='x2',zlim=range(Y_full))
    
    
    ########################## Vecchia ##########################
    # first setting with m_seq = c(10, 30)
    print("Vecchia1")
    
    time_Vecchia1=system.time({
      fit_Vecchia <- fit_model(y=Y_obs, locs=input[-inside_ind,],X=rep(1,N0), 
                               "matern25_isotropic",m_seq = c(10,30),silent=F)
    })
    
    time_Vecchia_pred1=system.time({
      pred_Vecchia=predictions(fit = fit_Vecchia,locs_pred=input,
                               X_pred=rep(1,N))
    })
    
    time_rec[[i_s]][2,1,i_repeat] = time_Vecchia1[3]
    time_rec[[i_s]][2,2,i_repeat] = time_Vecchia_pred1[3]
    
    #overall NRMSE
    NRMSE_rec[[i_s]][2,1,i_repeat]=sqrt(mean((pred_Vecchia-Z_full)^2))/sd(Z_full)
    #missing NRMSE
    NRMSE_rec[[i_s]][2,2,i_repeat]=sqrt(mean((pred_Vecchia[inside_ind]-Z_full[inside_ind])^2))/sd(Z_full)
    
    
    # second setting with m_seq = c(30,90)
    print("Vecchia2")
    
    time_Vecchia2=system.time({
      fit_Vecchia <- fit_model(y=Y_obs, locs=input[-inside_ind,],X=rep(1,N0), 
                               "matern25_isotropic",m_seq = c(30,90),silent=F) 
    })
    
    time_Vecchia_pred2=system.time({
      pred_Vecchia=predictions(fit = fit_Vecchia,locs_pred=input,
                               X_pred=rep(1,N),m = 150)
    })
    
    time_rec[[i_s]][3,1,i_repeat] = time_Vecchia2[3]
    time_rec[[i_s]][3,2,i_repeat] = time_Vecchia_pred2[3]
    
    #overall NRMSE
    NRMSE_rec[[i_s]][3,1,i_repeat]=sqrt(mean((pred_Vecchia-Z_full)^2))/sd(Z_full)
    #missing NRMSE
    NRMSE_rec[[i_s]][3,2,i_repeat]=sqrt(mean((pred_Vecchia[inside_ind]-Z_full[inside_ind])^2))/sd(Z_full)
    
    
    ########################## scaled Vecchia ##########################
    # first setting with m_s = c(30)
    print("Svecchia1")
    
    time_SVecchia1=system.time({
      fit_scaled_vecchia=fit_scaled(Y_obs,input[-inside_ind,],ms=c(30),trend='zero', nu = 2.5,nug=NULL)
    })
    
    time_SVecchia_pred1=system.time({
      pred_scaled_vecchia=predictions_scaled(fit_scaled_vecchia,input,m=60)
    })
    
    time_rec[[i_s]][4,1,i_repeat] = time_SVecchia1[3]
    time_rec[[i_s]][4,2,i_repeat] = time_SVecchia_pred1[3]
    
    
    #overall NRMSE
    NRMSE_rec[[i_s]][4,1,i_repeat]=sqrt(mean((pred_scaled_vecchia-Z_full)^2))/sd(Z_full)
    #missing NRMSE
    NRMSE_rec[[i_s]][4,2,i_repeat]=sqrt(mean((pred_scaled_vecchia[inside_ind]-Z_full[inside_ind])^2))/sd(Z_full)
    
    
    # second setting with m_s = c(90)
    print("Svecchia2")
    
    time_SVecchia2=system.time({
      fit_scaled_vecchia=fit_scaled(Y_obs,input[-inside_ind,],ms=c(90),trend='zero', nu = 2.5,nug=NULL)
    })
    
    time_SVecchia_pred2=system.time({
      pred_scaled_vecchia=predictions_scaled(fit_scaled_vecchia,input,m=150)
    })
    
    time_rec[[i_s]][5,1,i_repeat] = time_SVecchia2[3]
    time_rec[[i_s]][5,2,i_repeat] = time_SVecchia_pred2[3]
    
    
    #overall NRMSE
    NRMSE_rec[[i_s]][5,1,i_repeat]=sqrt(mean((pred_scaled_vecchia-Z_full)^2))/sd(Z_full)
    #missing NRMSE
    NRMSE_rec[[i_s]][5,2,i_repeat]=sqrt(mean((pred_scaled_vecchia[inside_ind]-Z_full[inside_ind])^2))/sd(Z_full)
    
    
    ########################## NNGP ##########################
    print("NNGP")
    
    cov.model <- "matern" #"matern", "exponential"
    
    sigma.sq <- 100
    
    sigma.sq.IG <- c(2, sigma.sq)
    
    g <-  5
    theta.alpha <- as.matrix(expand.grid("phi"=seq(.01, 0.5, length.out=g), 
                                         "alpha"=seq(0.01/sigma.sq, 0.2/sigma.sq, length.out=g),
                                         "nu"=seq(2.5, 2.5, by=1)))
    
    n.neighbors = 30
    
    time_NNGP=system.time({
      m.c <- spConjNNGP(Y_obs~1, coords=input[-inside_ind,], n.neighbors = n.neighbors,
                        k.fold = 5, score.rule = "crps",
                        n.omp.threads = 1,
                        theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG,
                        cov.model = cov.model,verbose=F)
    })
    
    
    ##prediction
    theta.alpha <- m.c$theta.alpha
    #names(theta.alpha) <- c("phi", "alpha", "nu")
    
    time_NNGP_pred=system.time({
      p_NNGP = spConjNNGP(Y_obs~1, coords=input[-inside_ind,], n.neighbors = n.neighbors,
                          X.0=matrix(1,N), coords.0=input,
                          n.omp.threads = 1,
                          theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG,
                          cov.model = cov.model,verbose=F)
    })
    
    
    time_rec[[i_s]][6,1,i_repeat] = time_NNGP[3]
    time_rec[[i_s]][6,2,i_repeat] = time_NNGP_pred[3]
    
    #overall NRMSE
    NRMSE_rec[[i_s]][6,1,i_repeat]=sqrt(mean((p_NNGP$y.0.hat-Z_full)^2))/sd(Z_full)
    #missing NRMSE
    NRMSE_rec[[i_s]][6,2,i_repeat]=sqrt(mean((p_NNGP$y.0.hat[inside_ind]-Z_full[inside_ind])^2))/sd(Z_full)
    
    
    
    ########################## laGP ##########################
    print("laGP")
    
    time_laGP=system.time({
      m_laGP_nn = aGP(X=input[-inside_ind,], Z=Y_obs, XX=input, method="nn", verb=0, end = 70)
    })
    
    time_rec[[i_s]][7,1,i_repeat] = time_laGP[3]
    
    
    #overall NRMSE
    NRMSE_rec[[i_s]][7,1,i_repeat]=sqrt(mean((m_laGP_nn$mean-Z_full)^2))/sd(Z_full)
    #missing NRMSE
    NRMSE_rec[[i_s]][7,2,i_repeat]=sqrt(mean((m_laGP_nn$mean[inside_ind]-Z_full[inside_ind])^2))/sd(Z_full)
    

  }
}



#save.image("brainin_simulation.RData")



#######plot
library(reshape2)
library(tidyverse)
library(plot3D)

load("brainin_simulation.RData")


# plot
i_s = 1
zlim = range(obs_IKF_pred_rec[[i_s]][1,],na.rm = T)

dat = data.frame(x1 = input[,1], x2 = input[,2], 
                 obs = obs_IKF_pred_rec[[i_s]][1,], IKF_pred = obs_IKF_pred_rec[[i_s]][2,], truth = Z_full)

pdf(paste0("plots/branin_obs",i_s,".pdf"),width=3.8,height=3)
ggplot(dat) + geom_raster(aes(x1,x2,fill=obs)) +
  theme_bw() +
  scale_fill_distiller(palette="RdYlGn",name="",limits=zlim) + 
  ggtitle("Observation")+
  xlab(expression(x[1])) +
  ylab(expression(x[2]))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.text=element_text(size=12),
        #panel.background = element_blank(),
        panel.grid.major = element_line(linewidth = .3),panel.grid.minor = element_blank())
dev.off()

pdf(paste0("plots/branin_pred",i_s,".pdf"),width=3.8,height=3)
ggplot(dat) + geom_raster(aes(x1,x2,fill=IKF_pred)) +
  theme_bw() +
  scale_fill_distiller(palette="RdYlGn",name="",limits=zlim) + 
  ggtitle("IKF-CG Prediction")+
  xlab(expression(x[1])) +
  ylab(expression(x[2]))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.text=element_text(size=12),
        #panel.background = element_blank(),
        panel.grid.major = element_line(linewidth = .3),panel.grid.minor = element_blank())
dev.off()

pdf("plots/branin_truth.pdf",width=3.8,height=3)
ggplot(dat) + geom_raster(aes(x1,x2,fill=truth)) +
  theme_bw() +
  scale_fill_distiller(palette="RdYlGn",name="",limits=zlim) + 
  ggtitle("Truth")+
  xlab(expression(x[1])) +
  ylab(expression(x[2]))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.text=element_text(size=12),
        #panel.background = element_blank(),
        panel.grid.major = element_line(linewidth = .3),panel.grid.minor = element_blank())
dev.off()




# mask1 NRMSE


NRMSE_mask1_df <- melt(NRMSE_rec[[1]], varnames = c("Method", "Metric", "Repeat"), value.name = "NRMSE")

pdf("plots/NRMSE_mask1.pdf",width=2.5,height=2)
NRMSE_mask1_df %>% 
  filter(Method != "Vecchia1" & Method != "Svecchia1") %>%
  filter(Metric == "NRMSE_miss") %>%
  ggplot(aes(Method, NRMSE))+
  geom_violin(color = "#0072B2",scale = "area")+ #aes(color=Method)
  scale_x_discrete(labels = c("Vecchia2" = "Vecchia", "Svecchia2" = "SVecchia", "IKF" = "IKF-CG"))+
  #scale_y_continuous(breaks=seq(0,0.9,0.3))+
  stat_summary(fun=mean, geom="point", shape=20, size=1,color = "#6A3D9A")+ #aes(color=Method)
  #scale_color_brewer(palette="Accent")+
  #geom_boxplot(width=0.1)+
  #geom_point(size = 0.2)+
  coord_flip()+
  theme_linedraw()+
  theme(legend.position = "none",legend.text=element_text(size=10),
        axis.title.y=element_blank(),#axis.text.y = element_text(size=10), #, angle=20
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()



time_mask1_df <- melt(apply(time_rec[[1]], c(1,3), function(x) sum(x,na.rm=T)), varnames = c("Method", "Repeat"), value.name = "time")

time_mask1_df %>% 
  filter(Method != "Vecchia1" & Method != "Svecchia1") %>%
  ggplot(aes(Method, time))+
  geom_violin(aes(color=Method),scale = "area")+
  scale_x_discrete(labels = c("Vecchia2" = "Vecchia", "Svecchia2" = "Svecchia"))+
  #facet_wrap(~Metric)
  theme_linedraw()+
  theme(legend.position = "none",axis.title.x=element_blank(),
        legend.text=element_text(size=10),axis.text=element_text(size=12),axis.title=element_text(size=13),
        panel.grid.major = element_line(colour = "gray80",linetype='dotted'),panel.grid.minor = element_blank())




# mask2 NRMSE
NRMSE_mask2_df <- melt(NRMSE_rec[[2]], varnames = c("Method", "Metric", "Repeat"), value.name = "NRMSE")

pdf("plots/NRMSE_mask2.pdf",width=2.5,height=2)
NRMSE_mask2_df %>% 
  filter(Method != "Vecchia1" & Method != "Svecchia1") %>%
  filter(Metric == "NRMSE_miss") %>%
  ggplot(aes(Method, NRMSE))+
  geom_violin(color = "#0072B2",scale = "area")+ #aes(color=Method)
  scale_x_discrete(labels = c("Vecchia2" = "Vecchia", "Svecchia2" = "SVecchia", "IKF" = "IKF-CG"))+
  #scale_y_continuous(breaks=seq(0,1.2,0.4))+
  stat_summary(fun=mean, geom="point", shape=20, size=1,color = "#6A3D9A")+ #aes(color=Method)
  #scale_color_brewer(palette="Accent")+
  #geom_boxplot(width=0.1)+
  #geom_point(size = 0.2)+
  coord_flip()+
  theme_linedraw()+
  theme(legend.position = "none",legend.text=element_text(size=10),
        axis.title.y=element_blank(),#axis.text.y = element_text(size=10), #, angle=20
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()

time_mask2_df <- melt(apply(time_rec[[2]], c(1,3), function(x) sum(x,na.rm=T)), varnames = c("Method", "Repeat"), value.name = "time")

time_mask2_df %>% 
  filter(Method != "Vecchia1" & Method != "Svecchia1") %>%
  ggplot(aes(Method, time))+
  geom_violin(aes(color=Method),scale = "area")+
  scale_x_discrete(labels = c("Vecchia2" = "Vecchia", "Svecchia2" = "Svecchia"))+
  #facet_wrap(~Metric)
  theme_linedraw()+
  theme(legend.position = "none",axis.title.x=element_blank(),
        legend.text=element_text(size=10),axis.text=element_text(size=12),axis.title=element_text(size=13),
        panel.grid.major = element_line(colour = "gray80",linetype='dotted'),panel.grid.minor = element_blank())



# random NRMSE
NRMSE_random_df <- melt(NRMSE_rec[[3]], varnames = c("Method", "Metric", "Repeat"), value.name = "NRMSE")

pdf("plots/NRMSE_random.pdf",width=2.5,height=2)
NRMSE_random_df %>% 
  filter(Method != "Vecchia1" & Method != "Svecchia1") %>%
  filter(Metric == "NRMSE_miss") %>%
  ggplot(aes(Method, NRMSE))+
  geom_violin(color = "#0072B2",scale = "area")+ #aes(color=Method)
  scale_x_discrete(labels = c("Vecchia2" = "Vecchia", "Svecchia2" = "SVecchia", "IKF" = "IKF-CG"))+
  #scale_y_continuous(breaks=seq(0,1.2,0.4))+
  stat_summary(fun=mean, geom="point", shape=20, size=1,color = "#6A3D9A")+ #aes(color=Method)
  #scale_color_brewer(palette="Accent")+
  #geom_boxplot(width=0.1)+
  #geom_point(size = 0.2)+
  coord_flip()+
  theme_linedraw()+
  theme(legend.position = "none",legend.text=element_text(size=10),
        axis.title.y=element_blank(),#axis.text.y = element_text(size=10), #, angle=20
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()



#mean time_total
time_mask1_mean = apply(apply(time_rec[[1]],c(1,3),function(x) sum(x,na.rm=T)),1,mean)
time_mask2_mean = apply(apply(time_rec[[2]],c(1,3),function(x) sum(x,na.rm=T)),1,mean)
time_random_mean = apply(apply(time_rec[[3]],c(1,3),function(x) sum(x,na.rm=T)),1,mean)

cbind(Disk1=time_mask1_mean,Disk2=time_mask2_mean,Random=time_random_mean)




