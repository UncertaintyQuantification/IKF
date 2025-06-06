library(FastGaSP)
library(RobustGaSP)
library(tidyverse)
source('../functions/functions_particle.R')

exp=6
n_repeat=20
i_repeat=3
set.seed((exp-1)*n_repeat+i_repeat+1)


kernel_type='matern_5_2'

# the setting of the simulation
n_t=900 #number of particles
T_sim=10  # time steps
sigma_0=.1


h = 0.1 ##time interval between each frame 
cut_r=.5 ##true cutoff


vx_abs=0.5
vy_abs=0.5
v_abs=sqrt(vx_abs^2+vy_abs^2)

sim = simulate_particle(v_abs=v_abs, n_t=n_t, T_sim=T_sim, 
                        h=h, cut_r=cut_r, sigma_0=sigma_0, model = 'unnormalized_Vicsek')



D_y=sim@D_y
N=n_t*T_sim*D_y

p_bar=2 ##expected number of neighbor
tol=10^{-8}*N*p_bar ##tolerance

testing_n = 200
testing_input=seq(-1,1,length.out=testing_n)
testing_inputs=matrix(testing_input)

cut_r_max=1.5
param_ini = log(c(0.3,1000,0.3/(cut_r_max-0.3)))  # Initial parameter values

est = fit(sim,param=param_ini,cut_r_max=cut_r_max, 
          kernel_type=kernel_type,tol=tol,testing_inputs = testing_inputs)


beta=est@parameters['beta']
tau=est@parameters['tau']
radius=est@parameters['radius']

pred_mean = est@predictions$mean
lower95 = est@predictions$lower95
upper95 = est@predictions$upper95

testing_output=testing_input 




dat_kernel=data.frame(input=rep(testing_input,2),
                      output=c(testing_output,pred_mean),
                      type=factor(rep(c("Truth","Prediction"),each=length(testing_input)),levels=c("Truth","Prediction")))
dat_CI=data.frame(x=c(testing_input,rev(testing_input)),
                  y=c(lower95,rev(upper95)))
pdf('plots/unnormalized_Vicsek_pred_CI.pdf',height=3.5,width=4.5)
ggplot(dat_CI)+geom_polygon(aes(x=x,y=y),fill="grey80")+
  geom_line(data=dat_kernel,aes(x=input,y=output,linetype=type,color=type))+
  xlab("d")+ylab("z(d)")+
  ggtitle("Interaction kernel")+
  scale_linetype_manual(values = c("solid","dashed"))+
  scale_color_manual(values=c("#fda6b1","#00618e"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = c(.8, .2),legend.key = element_rect(fill = "white", colour = "gray50"),
        legend.text=element_text(size=10),legend.title=element_blank(),
        #legend.box.background = element_rect(),legend.box.margin = margin(0,0,0,0),
        axis.text=element_text(size=12),axis.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5,size = 18))
dev.off()

plot_ind = 1:21
dat_kernel_insert=data.frame(input=rep(testing_input[plot_ind],2),
                             output=c(testing_output[plot_ind],pred_mean[plot_ind]),
                             type=factor(rep(c("Truth","Prediction"),each=length(testing_input[plot_ind])),levels=c("Truth","Prediction")))
dat_CI_insert=data.frame(x=c(testing_input[plot_ind],rev(testing_input[plot_ind])),
                         y=c(lower95[plot_ind],rev(upper95[plot_ind])))
pdf('plots/unnormalized_Vicsek_pred_insert_CI.pdf',height=2.8,width=3.6)
ggplot(dat_CI_insert)+geom_polygon(aes(x=x,y=y),fill="grey90")+
  geom_line(data=dat_kernel_insert,aes(x=input,y=output,linetype=type,color=type),size=0.6)+
  #xlab("d")+ylab("z(d)")+
  scale_linetype_manual(values = c("solid","dashed"))+
  scale_color_manual(values=c("#fda6b1","#00618e"))+
  scale_x_continuous(breaks = seq(-1, -0.8, by = 0.1))+
  scale_y_continuous(breaks = seq(-1, 0.8, by = 0.1))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 15),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
dev.off()




########### predict trajectory ########
T_test=50
weights = est@gp_weights
v_train_vec=est@training_data$training_velocity


p1=as.vector(rbind(sim@px_list[[T_sim+1]],sim@py_list[[T_sim+1]]))
v1=as.vector(rbind(sim@vx_list[[T_sim+1]],sim@vy_list[[T_sim+1]]))

set.seed((exp-1)*n_repeat+i_repeat+1)
m_unnormalized_Vicsek_test=unnormalized_Vicsek(p0=p1,v0=v1,n_t=n_t,T_sim=T_test,h=h,cut_r=cut_r,sigma_0=sigma_0)

input_pos_test = m_unnormalized_Vicsek_test$pos
v_test = m_unnormalized_Vicsek_test$v


est_param=c(beta,tau,radius)
set.seed((exp-1)*n_repeat+i_repeat+1)
m_unnormalized_Vicsek_pred=unnormalized_Vicsek_pred_trajec(est_param,v1,p1,v_abs,n_t,T_test,h,sigma_0,weights,v_train_vec,kernel_type = kernel_type)


pos_pred=m_unnormalized_Vicsek_pred$pos_pred
v_pred=m_unnormalized_Vicsek_pred$v_pred


## plot
plot_gap=20
plot_id=c(1:(n_t/plot_gap))*plot_gap
dat_path_truth=data.frame(id=rep(plot_id,each=T_test+1),
                          path=factor(rep(c("Truth"),each=(T_test+1)*length(plot_id)),levels=c("Truth","Prediction")),
                          s1=as.vector(t(input_pos_test[((plot_id-1)*2+1),])),
                          s2=as.vector(t(input_pos_test[(plot_id*2),])))
dat_arrow_truth=data.frame(id=plot_id,
                           x=input_pos_test[c((plot_id-1)*2+1),T_test],
                           y=input_pos_test[c(plot_id*2),T_test],
                           xend=input_pos_test[c((plot_id-1)*2+1),T_test]+(v_test[c((plot_id-1)*2+1),T_test+1])*0.3,
                           yend=input_pos_test[c(plot_id*2),T_test]+(v_test[c(plot_id*2),T_test+1])*0.3)

dat_path_pred=data.frame(id=rep(plot_id,each=T_test+1),
                         path=factor(rep(c("Prediction"),each=(T_test+1)*length(plot_id)),levels=c("Truth","Prediction")),
                         s1=as.vector(t(pos_pred[((plot_id-1)*2+1),])),
                         s2=as.vector(t(pos_pred[(plot_id*2),])))
dat_arrow_pred=data.frame(id=plot_id,
                          x=pos_pred[c((plot_id-1)*2+1),T_test],
                          y=pos_pred[c(plot_id*2),T_test],
                          xend=pos_pred[c((plot_id-1)*2+1),T_test]+(v_pred[c((plot_id-1)*2+1),T_test+1])*0.3,
                          yend=pos_pred[c(plot_id*2),T_test]+(v_pred[c(plot_id*2),T_test+1])*0.3)

pdf(paste("plots/unnormalized_Vicsek_pred_traj_",T_test,"_steps.pdf",sep=""),height=3.5,width=4.5)
ggplot(dat_path_truth,aes(s1,s2,group=id))+geom_path(aes(linetype=path,col=path),size=0.65)+
  geom_path(data=dat_path_pred,aes(linetype=path,col=path),size=0.6)+
  geom_segment(data=dat_arrow_truth,aes(x=x,y=y,xend=xend,yend=yend),col="#fda6b1",
               arrow = arrow(length = unit(0.03, "inches"),type="closed"))+
  geom_segment(data=dat_arrow_pred,aes(x=x,y=y,xend=xend,yend=yend),col="#00618e",
               arrow = arrow(length = unit(0.03, "inches"),type="closed"))+
  xlim(c(-2,33))+ylim(c(0,34))+
  xlab(expression(s[1]))+ylab(expression(s[2]))+
  ggtitle("Particle trajectories")+
  scale_linetype_manual(values = c("solid","dashed"))+
  scale_color_manual(values=c("#fda6b1","#00618e"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = c(.37, .89),legend.key = element_rect(fill = "white", colour = "gray50"),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.text=element_text(size=10),legend.title=element_blank(),
        #legend.box.background = element_rect(),legend.box.margin = margin(0,0,0,0),
        axis.text=element_text(size=12),axis.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5,size = 18))
dev.off()



