library(FastGaSP)
library(RobustGaSP)
source('../functions/functions_particle.R')

kernel_type='matern_5_2'

# the setting of the simulation
n_t=100 #number of particles
T_sim=5  # time steps
sigma_0_seq=c(.1,.2)


h = 0.1 ##time interval between each frame 
cut_r=.5 ##true cutoff


B=100
n_repeat=20

beta_est_rec=array(NA,c(3,n_repeat,length(sigma_0_seq)),
                   dimnames = list(c("est","LB","UB"),paste("Exp",1:n_repeat,sep = ""),c("sigma0.1","sigma0.2")))

tau_est_rec=array(NA,c(3,n_repeat,length(sigma_0_seq)),
                  dimnames = list(c("est","LB","UB"),paste("Exp",1:n_repeat,sep = ""),c("sigma0.1","sigma0.2")))

radius_est_rec=array(NA,c(3,n_repeat,length(sigma_0_seq)),
                     dimnames = list(c("est","LB","UB"),paste("Exp",1:n_repeat,sep = ""),c("sigma0.1","sigma0.2")))



for(s in 1:2){
  sigma_0=sigma_0_seq[s]
  for(i_repeat in 1:n_repeat){
    print(c(s,i_repeat))
    set.seed((s-1)*n_repeat+i_repeat+1)
    
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
    
    # fit
    cut_r_max=1.5
    param_ini = log(c(0.3,1000,0.3/(cut_r_max-0.3)))  # Initial parameter values
    

    est = fit(sim,param=param_ini,cut_r_max=cut_r_max, 
              kernel_type=kernel_type,tol=tol)
    
    
    rb = residual_bootstrap_unnormalized_Vicsek(est, sim, B = B,
                                                param_ini=param_ini, cut_r_max=cut_r_max,
                                                kernel_type = kernel_type, 
                                                tol = tol) 
    
    beta_est_rec[,i_repeat,s] = rb$beta_est_rec
    tau_est_rec[,i_repeat,s] = rb$tau_est_rec
    radius_est_rec[,i_repeat,s] = rb$radius_est_rec
    
    
  }
  
  
}


radius_est_rec

#save.image("unnormalized_Vicsek_residual_boostrap.RData")

#load("unnormalized_Vicsek_residual_boostrap.RData")








####### plot
library(tidyverse)
dat_radius0.1=data.frame(simulation=1:n_repeat,
                         est=radius_est_rec[1,,1],
                         LB=radius_est_rec[2,,1],
                         UB=radius_est_rec[3,,1])
dat_radius0.1$covered=factor(dat_radius0.1$LB<0.5 & dat_radius0.1$UB>0.5,levels = c("TRUE","FALSE"))
jpeg("plots/radius_CI_s0.1.jpeg", quality = 100,width = 1000, height = 600, res = 250)
ggplot(dat_radius0.1,aes(x=simulation,y=est))+
  geom_errorbar(aes(ymin=LB, ymax=UB,color=covered), width=.5)+
  geom_point(aes(fill=covered), shape=21, stroke=0)+
  scale_fill_manual(values=c("#037b83","#b69100")) +
  scale_colour_manual(values=c("#00AFBB","#E7B800")) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "purple")+
  xlab("Simulation")+ylab(expression(hat(r)))+
  ggtitle(expression(sigma[0]^2==0.1^2))+
  theme_linedraw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")
dev.off()

#(bootstrap residual) CI for radius with sigma_0=0.2
dat_radius0.2=data.frame(simulation=1:n_repeat,
                         est=radius_est_rec[1,,2],
                         LB=radius_est_rec[2,,2],
                         UB=radius_est_rec[3,,2])

dat_radius0.2$covered=factor(dat_radius0.2$LB<0.5 & dat_radius0.2$UB>0.5,levels = c("TRUE","FALSE"))

jpeg("plots/radius_CI_s0.2.jpeg", quality = 100,width = 1000, height = 600, res = 250)
ggplot(dat_radius0.2,aes(x=simulation,y=est))+
  geom_errorbar(aes(ymin=LB, ymax=UB,color=covered), width=.5)+
  geom_point(aes(fill=covered), shape=21, stroke=0)+
  scale_fill_manual(values=c("#037b83","#b69100")) +
  scale_colour_manual(values=c("#00AFBB","#E7B800")) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "purple")+
  xlab("Simulation")+ylab(expression(hat(r)))+
  ggtitle(expression(sigma[0]^2==0.2^2))+
  theme_linedraw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.position = c(.15, .22),#legend.direction="horizontal",
        legend.text=element_text(size=7),legend.title = element_text(size=9),
        legend.key.spacing.y = unit(-0.15, 'cm'),
        #legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.background=element_rect(fill = alpha("white", 0)),
        plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")
dev.off()




