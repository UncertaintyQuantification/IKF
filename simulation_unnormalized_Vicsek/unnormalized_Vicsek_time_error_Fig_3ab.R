library(FastGaSP)
library(RobustGaSP)
source('../functions/functions_particle.R')

library(Rcpp)
library(RcppEigen)
sourceCpp(file='src/functions.cpp')  ###

kernel_type='matern_5_2'

# the setting of the simulation
n_t=100 #number of particles
T_sim_seq = seq(4,200,4)
sigma_0=.2


h = 0.1 ##time interval between each frame 
cut_r=.5 ##true cutoff



######## time comparison ######## 

n_repeat=20

time_rec_all = time_rec_pred = matrix(0,length(T_sim_seq),3,dimnames = list(NULL,c("IKF-CG","CG","Direct")))

####some common parameters
testing_n = 100
testing_input=as.numeric(seq(-1,1,length.out = testing_n))


testing_output=testing_input 


for(t in 1:length(T_sim_seq)){
  print(t)
  
  T_sim=T_sim_seq[t]

  
  for(i_repeat in 1: n_repeat){
    # set.seed(t)
    
    set.seed(t*n_repeat+i_repeat)
    
    # simulation
    vx_abs=0.5
    vy_abs=0.5
    v_abs=sqrt(vx_abs^2+vy_abs^2)
    
    
    sim = simulate_particle(v_abs=v_abs, n_t=n_t, T_sim=T_sim, 
                            h=h, cut_r=cut_r, sigma_0=sigma_0, model = 'unnormalized_Vicsek')
    
    D_y=sim@D_y
    N=n_t*T_sim*D_y
    p_bar=2 ##expected number of neighbor
    tol=10^{-8}*N*p_bar ##tolerance
    
    cut_r_max=1.5
    param = log(c(0.1,10,cut_r/(cut_r_max-cut_r)))
    
    t1=system.time({
      est_IKF = model_pred_fixed_param_IKF(sim,param=param,cut_r_max=cut_r_max,
                                           kernel_type=kernel_type, tol=tol,testing_input=testing_input)
    })
    
    time_rec_all[t,1]=time_rec_all[t,1]+t1[3]
    time_rec_pred[t,1]= time_rec_pred[t,1]+est_IKF$time_rec[3]
    

    # plot(est_IKF@predictions$mean,testing_output,type='l')
    # abline(0,1,col=2)
    # sqrt(mean((est_IKF@predictions$mean-testing_output)^2))/sd(testing_output)
    
    if(T_sim_seq[t]<=80){
      t2=system.time({
        est_CG = model_pred_fixed_param_CG(sim,param=param,cut_r_max=cut_r_max,
                                           kernel_type=kernel_type, tol=tol,testing_input=testing_input)
      })
      time_rec_all[t,2]=time_rec_all[t,2]+t2[3]
      time_rec_pred[t,2]=time_rec_pred[t,2]+est_CG$time_rec[3]
    }

    if(T_sim_seq[t]<=20){
      t3=system.time({
        est_direct = model_pred_fixed_param_direct(sim,param=param,cut_r_max=cut_r_max,
                                                   kernel_type=kernel_type,testing_input=testing_input)
      })
      time_rec_all[t,3]=time_rec_all[t,3]+t3[3]
      time_rec_pred[t,3]=time_rec_pred[t,3]+est_direct$time_rec[3]
    }

  }
  
  time_rec_all[t,]=time_rec_all[t,]/n_repeat
  time_rec_pred[t,]=time_rec_pred[t,]/n_repeat
}





######## error comparison ######## 

T_sim_seq_error=seq(4,20,1)

RMSE_rec = matrix(0,length(T_sim_seq_error),3,dimnames = list(NULL,c("IKF-CG","Direct","diff")))

n_repeat_error=20

for(t in 1:length(T_sim_seq_error)){
  print(t)
  T_sim = T_sim_seq_error[t]

  
  for(i_repeat in 1:n_repeat_error){
    set.seed((t-1)*n_repeat_error+i_repeat)
    
    # simulation
    vx_abs=0.5
    vy_abs=0.5
    v_abs=sqrt(vx_abs^2+vy_abs^2)
    
    sim = simulate_unnormalized_Vicsek(v_abs=v_abs, n_t=n_t, T_sim=T_sim, 
                                       h=h, cut_r=cut_r, sigma_0=sigma_0)
    
    D_y=sim@D_y
    N=n_t*T_sim*D_y
    p_bar=2 ##expected number of neighbor
    tol=10^{-8}*N*p_bar ##tolerance
    
    cut_r_max=1.5
    param = log(c(0.1,10,cut_r/(cut_r_max-cut_r)))
    
    
    est_IKF = model_pred_fixed_param_IKF(sim,param=param,cut_r_max=cut_r_max,
                                         kernel_type=kernel_type, tol=tol,testing_input=testing_input)
    
    RMSE_rec[t,1]=RMSE_rec[t,1]+sqrt(mean((est_IKF$pred_mean_fast-testing_output)^2))
    
    est_direct = model_pred_fixed_param_direct(sim,param=param,cut_r_max=cut_r_max,
                                               kernel_type=kernel_type,testing_input=testing_input)
    
    RMSE_rec[t,2]=RMSE_rec[t,2]+sqrt(mean((est_direct$pred_mean_direct-testing_output)^2))
    
    
    RMSE_rec[t,3]=RMSE_rec[t,3]+sqrt(mean((est_IKF$pred_mean_fast-est_direct$pred_mean_direct)^2))
  }
  RMSE_rec[t,]=RMSE_rec[t,]/n_repeat_error
}


#save.image(file='unnormalized_Vicsek_time_error.RData')

# load('unnormalized_Vicsek_time_error.RData')


##for plots
library(tidyverse)
N_seq = n_t*T_sim_seq*D_y
N_seq_error = n_t*T_sim_seq_error*D_y

time_rec=time_rec_all
#time_rec=time_rec_pred
time_rec[time_rec==0]=NA


time_rec[5,3]=NA #remove the last point of direct computation

dat1=data.frame(N=rep(N_seq,3),
                method=factor(rep(c("Direct computation","CG","IKF-CG"),each=length(N_seq)),levels =c( "Direct computation","CG","IKF-CG")),
                time=c(time_rec[,3],time_rec[,2],time_rec[,1]))

pdf('plots/time_comparison_unnormalized_Vicsek.pdf',height=3.5,width=5)
ggplot(dat1,aes(N,time,shape=method,fill=method))+geom_point(size = 3)+
  xlab(expression(tilde(N)))+ylab('time (s)')+
  scale_shape_manual(values=c(24,22,21))+
  scale_fill_manual(values=c("#fbb4ae","#decbe4", "#80b1d3"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = "top",#legend.position = c(.72, .8),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.text=element_text(size=15),
        #legend.box.background = element_rect(),legend.box.margin = margin(0,0,0,0),
        axis.text=element_text(size=15),axis.title=element_text(size=20),
        legend.title=element_blank())
dev.off()


dat2=data.frame(N=rep(N_seq_error,3),
                method=factor(rep(c('Direct computation','IKF-CG','Difference'),each=length(N_seq_error)),levels = c('Direct computation','IKF-CG','Difference')),
                nrmse=c(RMSE_rec[,2],RMSE_rec[,1],RMSE_rec[,3])/sd(testing_output))
pdf('plots/NRMSE_unnormalized_Vicsek.pdf',height=3.5,width=5)
ggplot(dat2,aes(N,nrmse,shape=method,fill=method,size=method))+geom_point()+
  xlab(expression(tilde(N)))+ylab('NRMSE')+
  scale_shape_manual(values=c(24,21,23))+
  scale_fill_manual(values=c("#fbb4ae","#80b1d3", "#ffffb3"))+
  scale_size_manual(values=c(3,2.5,2.8))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = "top",#legend.position = c(.3, .4),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.text=element_text(size=15),
        #legend.box.background = element_rect(),legend.box.margin = margin(0,0,0,0),
        axis.text=element_text(size=15),axis.title=element_text(size=20),
        legend.title=element_blank())
dev.off()

legend_dat=data.frame(x=rnorm(4),y=rnorm(4),method=factor(c('Direct computation','CG','IKF-CG','Difference'),levels = c('Direct computation','CG','IKF-CG','Difference')))
pdf('plots/legend.pdf',height=3.5,width=15)
ggplot(legend_dat,aes(x,y))+geom_point(aes(shape=method,fill=method),size = 3)+
  scale_shape_manual(values=c(24,22,21,23))+
  scale_fill_manual(values=c("#fbb4ae","#decbe4","#80b1d3", "#ffffb3"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = "top",#legend.position = c(.3, .4),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.text=element_text(size=15),
        #legend.box.background = element_rect(),legend.box.margin = margin(0,0,0,0),
        axis.text=element_text(size=15),axis.title=element_text(size=20),
        legend.title=element_blank())
dev.off()






