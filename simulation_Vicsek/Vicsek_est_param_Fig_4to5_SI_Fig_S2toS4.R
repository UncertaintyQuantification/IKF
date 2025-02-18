library(reshape2)
library(ggplot2)
library(scales)
library(tidyverse)
library(FastGaSP)

####### choose the kernel type
kernel_type='matern_5_2'
#kernel_type='exp'


T_sim_seq = c(5,10)
n_t_seq = c(100, 300, 900)
#sigma_0=.5
sigma_0_seq=c(0.1,0.2)
#settings for all experiment
settings = expand.grid(T_sim=T_sim_seq ,n_t=n_t_seq,sigma_0=sigma_0_seq)
n_repeat=20


# save results
nrmse_rec = covered95_rec = length95_rec = 
  matrix(NA,dim(settings)[1],n_repeat,
         dimnames = list(paste("T",settings[,1],",n",settings[,2],",s",settings[,3],sep=""),
                         paste("Exp",1:n_repeat,sep = "")))

beta_rec = tau_rec = radius_rec = 
  matrix(NA,dim(settings)[1],n_repeat,
         dimnames = list(paste("T",settings[,1],",n",settings[,2],",s",settings[,3],sep=""),
                         paste("Exp",1:n_repeat,sep = "")))

# run simulation
for(i_s in 1:dim(settings)[1]){
  print(i_s)
  T_sim = settings[i_s,1]
  n_t = settings[i_s,2]
  sigma_0 = settings[i_s,3]
  
  D_y=1 ##even it is 2D, here y obs is 1 D for Vicsek
  N=n_t*T_sim*D_y ##this one is the number of observations, \tilde N in the paper
  cut_r=.5 ##true cutoff
  
  for(i_repeat in 1: n_repeat){
    set.seed((i_s-1)*n_repeat+i_repeat+1)
    
    # simulation
    vx_abs=0.5
    vy_abs=0.5
    v_abs=sqrt(vx_abs^2+vy_abs^2)
    
    sim = simulate_particle(v_abs=v_abs, n_t = n_t, T_sim = T_sim, 
                            cut_r = cut_r, sigma_0 = sigma_0,)
    
    # fit
    cut_r_max=1.5
    
    p_bar=2 ##expected number of neighbor
    tol=10^{-8}*N*p_bar ##tolerance

    testing_n = 200
    testing_input = as.numeric(seq(-pi,pi,(2*pi)/(testing_n-1)))
    testing_inputs = matrix(testing_input)
    
    param_ini = log(c(0.3,1000,0.3/(cut_r_max-0.3)))  # Initial parameter values
    
    
    est = fit(sim,param=param_ini,cut_r_max=cut_r_max, 
              kernel_type = kernel_type, tol = tol,
              testing_inputs = testing_inputs)
    
    beta_rec[i_s,i_repeat]=est@parameters['beta']
    tau_rec[i_s,i_repeat]=est@parameters['tau']
    radius_rec[i_s,i_repeat]=est@parameters['radius']
    
    pred_mean = est@predictions$mean
    lower95 = est@predictions$lower95
    upper95 = est@predictions$upper95
    
    testing_output=testing_input 
    
    NRMSE = mean( (pred_mean-testing_output)^2)/sd(testing_output)
    coverage95 = mean(testing_output<upper95 & testing_output>lower95)
    length95=mean(upper95-lower95)
    
    nrmse_rec[i_s,i_repeat]=NRMSE
    covered95_rec[i_s,i_repeat]=coverage95
    length95_rec[i_s,i_repeat]=length95
    
    print(c(nrmse_rec[i_s,i_repeat],covered95_rec[i_s,i_repeat], length95_rec[i_s,i_repeat]))
    
  }
  
  
}

#save.image('record_est_param_Vicsek_matern.RData')
#save.image('record_est_param_Vicsek_exp.RData')



####plot

#load("record_est_param_Viscek_matern_5_2_kernel_Feb_18.RData")
load("record_est_param_Viscek_exp_kernel_Feb_18.RData")

#exp_set$T_sim=factor(exp_set$T_sim, levels = c(5, 10),labels = c("T = 5", "T = 10"))
exp_set$T_sim=as.factor(exp_set$T_sim)
#exp_set$n_t=as.factor(exp_set$n_t)
exp_set$n_t=factor(exp_set$n_t, levels = c(100,300,900),
                   labels= c("n[p]==100","n[p]==300","n[p]==900"))
exp_set$sigma_0=factor(exp_set$sigma_0, levels = c(0.1,0.2),
                       labels= c("sigma[0]==0.1","sigma[0]==0.2"))


radius=cbind(exp_set,radius_rec)
radius_plot=melt(radius,id=colnames(exp_set))

nrmse=cbind(exp_set,nrmse_rec)
nrmse_plot=melt(nrmse,id=colnames(exp_set))

cov95=cbind(exp_set,data.frame(L=apply(length95_rec,1,mean),
                               deltaP=0.95-apply(covered95_rec,1,mean),
                               P=apply(covered95_rec,1,mean)))


######separate the plots for sigma = 0.1 and sigma = 0.2
chosen_sigma = 0.1
jpeg(paste("plots/Vicsek_radius_sigma",chosen_sigma,"_",kernel_type,".jpeg",sep=""), quality = 100,width = 1000, height = 580, res = 250)
radius_plot %>% filter(sigma_0 == paste0("sigma[0]==", chosen_sigma)) %>% 
  ggplot(aes(T_sim,value))+
  geom_boxplot(aes(fill=T_sim))+facet_grid(~ n_t, labeller = label_parsed)+
  geom_hline(yintercept=0.5, linetype="dashed", color = "purple")+
  scale_fill_manual(values=c("#8dd3c7","#fed9a6"))+
  xlab(expression(n[tau]))+
  ylab(expression(hat(r)))+
  #ggtitle(expression(sigma[0]==0.1))+
  ggtitle(bquote(sigma[0]^2 == .(chosen_sigma)^2))+
  #theme_bw()+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),axis.title=element_text(size=15),
        strip.background = element_rect(colour = "black", fill = "white",linewidth=0.5),
        strip.text.x = element_text(size = 11),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")#+
#scale_fill_brewer(palette="Pastel1")+
#labs(fill = "T")
dev.off()

jpeg(paste("plots/Vicsek_nrmse_sigma",chosen_sigma,"_",kernel_type,".jpeg",sep=""), quality = 100,width = 1000, height = 580, res = 250)
nrmse_plot %>% filter(sigma_0 == paste0("sigma[0]==", chosen_sigma)) %>% 
  ggplot(aes(T_sim,value))+
  geom_boxplot(aes(fill=T_sim))+facet_grid(~ n_t, labeller = label_parsed)+
  scale_fill_manual(values=c("#8dd3c7","#fed9a6"))+
  xlab(expression(n[tau]))+
  ylab("NRMSE")+
  scale_y_log10("NRMSE",breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits=range(nrmse_rec))+
  #ggtitle(expression(sigma[0]==0.1))+
  ggtitle(bquote(sigma[0]^2 == .(chosen_sigma)^2))+
  #theme_linedraw()+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),axis.title=element_text(size=15),
        strip.background = element_rect(colour = "black", fill = "white",linewidth=0.5),
        strip.text.x = element_text(size = 11),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")#+
#scale_fill_brewer(palette="Pastel1")+
#labs(fill = "T")
dev.off()

cov95_fix_sigma = cov95 %>% filter(sigma_0 == paste0("sigma[0]==", chosen_sigma))

ylim_L=c(0,max(cov95_fix_sigma$L))
#ylim_P=c(-0.05,0.03)
ylim_P=c(0.85,1.00)
b = diff(ylim_L)/diff(ylim_P)
a = ylim_L[1] - b*ylim_P[1]

jpeg(paste("plots/Vicsek_cov95_sigma",chosen_sigma,"_",kernel_type,".jpeg",sep=""), quality = 100,width = 1000, height = 580, res = 250)
cov95_fix_sigma %>% 
  ggplot(aes(x=T_sim))+
  #geom_bar(aes(y=L), stat="identity",fill="#FFCC99",alpha=0.9)+
  #geom_line(aes(y=a+deltaP*b,group = 1),color="#663300") +
  geom_bar(aes(y=L,fill=T_sim), stat="identity",alpha=0.9)+
  geom_point(aes(y=a+P*b,shape=T_sim,color=T_sim))+
  geom_hline(yintercept=a+0.95*b, linetype="dashed", color = "purple")+
  scale_fill_manual(values=c("#8dd3c7","#fed9a6"))+
  scale_color_manual(values=c("#1a6759","#e78600"))+
  xlab(expression(n[tau]))+
  facet_grid(~ n_t, labeller = label_parsed)+
  scale_y_continuous(name = expression(L(95*"%")),limits=ylim_L,
                     sec.axis = sec_axis(~(.-a)/b, name=expression(P(95*"%"))))+
  #ggtitle(expression(sigma[0]==0.1))+
  ggtitle(bquote(sigma[0]^2 == .(chosen_sigma)^2))+
  #theme_linedraw()+theme(plot.title = element_text(hjust = 0.5))+
  #scale_fill_brewer(palette="Pastel1")+scale_color_brewer(palette="Set1")+
  #labs(fill = "T",color="T",shape="T") + theme(legend.position = "none")
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),axis.title=element_text(size=15),
        strip.background = element_rect(colour = "black", fill = "white",linewidth=0.5),
        strip.text.x = element_text(size = 11),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
#theme_minimal()
dev.off()



