library(FastGaSP)



exp=6 # exp=6 corresponds to n_t=900, T_sim=10, sigma_0=0.1
n_repeat=20
i_repeat=1
set.seed((exp-1)*n_repeat+i_repeat+1)


kernel_type='matern_5_2'

n_t = 900
T_sim = 10
sigma_0=0.1

D_y=2 ##here y obs is 2 D
N=n_t*T_sim*D_y 
h = 0.1 ##simulation interval 
cut_r=.5 ##truth cutoff



vx_abs=0.5
vy_abs=0.5
v_abs=sqrt(vx_abs^2+vy_abs^2)

sim = simulate_particle(v_abs=v_abs, n_t = n_t, T_sim = T_sim, 
                        cut_r = cut_r, sigma_0 = sigma_0, model = 'two_interactions_Vicsek')


cut_r_max=1.5

p_bar=2 ##expected number of neighbor
tol=10^{-8}*N*p_bar ##tolerance


testing_n=200
testing_v_input=seq(min(c(unlist(sim@vx_list), unlist(sim@vy_list))),max(c(unlist(sim@vx_list), unlist(sim@vy_list))),length.out=testing_n)
testing_v_output=testing_v_input 

testing_d_input=seq(0.025,sim@radius,length.out=testing_n)
testing_f_output=f_Vicsek_variation(testing_d_input,r_max=sim@radius)

testing_inputs = cbind(testing_v_input, testing_d_input)

param_ini=log(c(0.01,10,10^4,10^4,0.8/(cut_r_max-0.8)))

est = fit(sim,param=param_ini,cut_r_max=cut_r_max, 
          kernel_type=kernel_type,tol=tol,testing_inputs = testing_inputs)

beta_v=est@parameters['beta_v']
beta_f=est@parameters['beta_f']
tau_v=est@parameters['tau_v']
tau_f=est@parameters['tau_f']
radius=est@parameters['radius']

pred_mean_v = est@predictions$mean_v
lower95_v = est@predictions$lower95_v
upper95_v = est@predictions$upper95_v

pred_mean_f = est@predictions$mean_f
lower95_f = est@predictions$lower95_f
upper95_f = est@predictions$upper95_f



# #save.image('Vicsek_variation_n900_T10_sig01_results.RData')
load('Vicsek_variation_n900_T10_sig01_results.RData')


dat_kernel1=data.frame(input=rep(testing_v_input,2),
                       output=c(testing_v_output,pred_mean_v),
                       type=factor(rep(c("Truth","Prediction"),each=length(testing_v_input)),levels=c("Truth","Prediction")))
dat_CI1=data.frame(x=c(testing_v_input,rev(testing_v_input)),
                   y=c(lower95_v,rev(upper95_v)))

#pdf('Vicsek_variation_pred_CI1.pdf',height=3.5,width=4.5)
ggplot(dat_CI1)+geom_polygon(aes(x=x,y=y),fill="grey80")+
  geom_line(data=dat_kernel1,aes(x=input,y=output,linetype=type,color=type))+
  xlab(expression(d[1]))+ylab(expression(z[1](d[1])))+
  ggtitle("First interaction kernel")+
  scale_linetype_manual(values = c("solid","dashed"))+
  scale_color_manual(values=c("#fda6b1","#00618e"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = c(.8, .2),legend.key = element_rect(fill = "white", colour = "gray50"),
        legend.text=element_text(size=10),legend.title=element_blank(),
        #legend.box.background = element_rect(),legend.box.margin = margin(0,0,0,0),
        axis.text=element_text(size=12),axis.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5,size = 18))
#dev.off()
plot_ind = 1:20
dat_kernel1_insert=data.frame(input=rep(testing_v_input[plot_ind],2),
                              output=c(testing_v_output[plot_ind],pred_mean_v[plot_ind]),
                              type=factor(rep(c("Truth","Prediction"),each=length(testing_v_input[plot_ind])),levels=c("Truth","Prediction")))
dat_CI1_insert=data.frame(x=c(testing_v_input[plot_ind],rev(testing_v_input[plot_ind])),
                          y=c(lower95_v[plot_ind],rev(upper95_v[plot_ind])))
#pdf('Vicsek_variation_pred_insert_CI1.pdf',height=2.8,width=3.6)
ggplot(dat_CI1_insert)+geom_polygon(aes(x=x,y=y),fill="grey90")+
  geom_line(data=dat_kernel1_insert,aes(x=input,y=output,linetype=type,color=type),size=0.6)+
  #xlab("d")+ylab("z(d)")+
  scale_linetype_manual(values = c("solid","dashed"))+
  scale_color_manual(values=c("#fda6b1","#00618e"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = "none",
        #legend.box.background = element_rect(),legend.box.margin = margin(0,0,0,0),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text=element_text(size=15),#axis.title=element_text(size=20),
        plot.title = element_text(hjust = 0.5,size = 15))
#dev.off()


dat_kernel2=data.frame(input=rep(testing_d_input,2),
                       output=c(testing_f_output,pred_mean_f),
                       type=factor(rep(c("Truth","Prediction"),each=length(testing_d_input)),levels=c("Truth","Prediction")))
dat_CI2=data.frame(x=c(testing_d_input,rev(testing_d_input)),
                   y=c(lower95_f,rev(upper95_f)))

#pdf('Vicsek_variation_pred_CI2.pdf',height=3.5,width=4.5)
ggplot(dat_CI2)+geom_polygon(aes(x=x,y=y),fill="grey80")+
  geom_line(data=dat_kernel2,aes(x=input,y=output,linetype=type,color=type))+
  #xlab("d")+ylab("z(d)")+
  xlab(expression(d[2]))+ylab(expression(z[2](d[2])))+
  ggtitle("Second interaction kernel")+
  scale_linetype_manual(values = c("solid","dashed"))+
  scale_color_manual(values=c("#fda6b1","#00618e"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = c(.8, .8),legend.key = element_rect(fill = "white", colour = "gray50"),
        legend.text=element_text(size=10),legend.title=element_blank(),
        #legend.box.background = element_rect(),legend.box.margin = margin(0,0,0,0),
        axis.text=element_text(size=12),axis.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5,size = 18))
#dev.off()
plot_ind = 1:20
dat_kernel2_insert=data.frame(input=rep(testing_d_input[plot_ind],2),
                              output=c(testing_f_output[plot_ind],pred_mean_f[plot_ind]),
                              type=factor(rep(c("Truth","Prediction"),each=length(testing_d_input[plot_ind])),levels=c("Truth","Prediction")))
dat_CI2_insert=data.frame(x=c(testing_d_input[plot_ind],rev(testing_d_input[plot_ind])),
                          y=c(lower95_f[plot_ind],rev(upper95_f[plot_ind])))
#pdf('Vicsek_variation_pred_insert_CI2.pdf',height=2.8,width=3.6)
ggplot(dat_CI2_insert)+geom_polygon(aes(x=x,y=y),fill="grey90")+
  geom_line(data=dat_kernel2_insert,aes(x=input,y=output,linetype=type,color=type),size=0.6)+
  #xlab("d")+ylab("z(d)")+
  scale_linetype_manual(values = c("solid","dashed"))+
  scale_color_manual(values=c("#fda6b1","#00618e"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = "none",
        #legend.box.background = element_rect(),legend.box.margin = margin(0,0,0,0),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text=element_text(size=15),#axis.title=element_text(size=20),
        plot.title = element_text(hjust = 0.5,size = 15))
#dev.off()






