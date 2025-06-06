library(FastGaSP)
library(RobustGaSP)
source('../functions/functions_particle.R')

n_t=100
T_sim_seq=seq(4,20,1)    # time steps
sigma_0=.2

h = 0.1 ##time interval between each frame 
cut_r=.5 ##true cutoff

tilde_nu = 0.1
maxIte=1000

log_det_rec=matrix(0,2,length(T_sim_seq),dimnames = list(c("Truth","Approx"),T_sim_seq))
log_lik_rec=matrix(0,2,length(T_sim_seq),dimnames = list(c("Truth","Approx"),T_sim_seq))


n_repeat=20

kernel_type='matern_5_2'



for(i_sim in 1:length(T_sim_seq)){
  print(i_sim)
  
  T_sim=T_sim_seq[i_sim]

  
  for(i_repeat in 1:n_repeat){
    set.seed(i_sim*n_repeat+i_repeat)
    
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
    
    est = fit(sim,param=param,cut_r_max=cut_r_max, est_param=FALSE,
              kernel_type=kernel_type,tol=tol, testing_input=NULL, compute_CI=FALSE)
    
    
    beta=exp(param[1])
    tau=exp(param[2])  ###sigma_2/sigma_2_0
    threshold_r=exp(param[3])/(1+exp(param[3])) * cut_r_max
    
    A_all_vec=est@training_data$A_v
    v_all_vec=est@training_data$training_velocity
    num_neighbors_all_vec=est@training_data$num_neighbors
    sort_d_all=sort(v_all_vec,index.return=T)
    N_tilde_all=length(v_all_vec) ###this is N_j in the paper
    
    
    delta_x_all=sort_d_all$x[-1]-sort_d_all$x[-N_tilde_all]
    
    
    output_all=as.vector(rbind(
      unlist(sim@vx_list[1+1:sim@T_time]),
      unlist(sim@vy_list[1+1:sim@T_time])
    ))
    
    d_sort_rev_ix=sort(sort_d_all$ix,index.return=T)$ix
    
    num_neighbors_v_vec=2*num_neighbors_all_vec
    A_full=matrix(0,D_y*n_t*T_sim,N_tilde_all) ##put zero for now, we only need to record nonzero
    count_y=1;
    neighbor_start=0
    for(t in 1:T_sim){
      for(i in 1:n_t){
        neighbor_index = neighbor_start+1:num_neighbors_v_vec[count_y]
        A_full[D_y*(count_y-1)+1:D_y,neighbor_index]=A_all_vec[(2*neighbor_index[1]-1):(2*neighbor_index[num_neighbors_v_vec[count_y]])]
        neighbor_start=neighbor_start+num_neighbors_v_vec[count_y]
        count_y=count_y+1
      }
    }
    
    
    R0=abs(outer(v_all_vec,v_all_vec,'-'))
    R=matern_5_2_funct(R0,beta_i=beta)
    
    #
    SS=A_full%*%R%*%t(A_full)*tau
    
    L_SS=t(chol(SS+diag(N)))
    true_log_det=2*sum(log(diag(L_SS)))
    
    
    
    
    
    l=20
    q=1
    Omega=matrix(rnorm(N*l),N,l)
    Y=matrix(rnorm(N*l),N,l)
    #Omega=matrix(rbinom(N*l,1,.5)*2-1,N,l)
    for(i in 1:q){
      if(i==1){
        #Y=SS%*%Omega
        #use IKF
        for(i_l in 1:l){
          A_t_O=A_t_times_x_particle(output=Omega[,i_l], A_all_v=A_all_vec, num_neighbors_vec=2*num_neighbors_all_vec, D_y=D_y, N_tilde=N_tilde_all)
          #A_t_O-t(A_full)%*%Omega[,i_l]
          #R_A_t_O_temp=R_times_z(param=log(c(beta,tilde_nu)), have_noise=T, delta_x=delta_x_all, z=A_t_O[sort_d_all$ix], kernel_type=kernel_type)-tilde_nu*A_t_O[sort_d_all$ix]
          R_A_t_O_temp=IKF(beta, tilde_nu,delta_x=delta_x_all,output=A_t_O[sort_d_all$ix], kernel_type=kernel_type)-tilde_nu*A_t_O[sort_d_all$ix]
          R_A_t_O=R_A_t_O_temp[d_sort_rev_ix]
          #R_A_t_O-R%*%t(A_full)%*%Omega[,i_l]
          A_R_A_t_O=A_times_x_particle(output=R_A_t_O, A_all_v=A_all_vec, num_neighbors_vec=2*num_neighbors_all_vec, D=D_y, N=N)
          #A_R_A_t_O-A_full%*%R%*%t(A_full)%*%Omega[,i_l]
          Y[,i_l]=A_R_A_t_O*tau
        }
      }else{
        #Y=SS%*%Y
        A_t_Y=A_t_times_x_particle(output=Y[,i_l], A_all_v=A_all_vec, num_neighbors_vec=2*num_neighbors_all_vec, D_y=D_y, N_tilde=N_tilde_all)
        #A_t_Y-t(A_full)%*%Y[,i_l]
        R_A_t_Y_temp=IKF(beta, tilde_nu,delta_x=delta_x_all,output=A_t_Y[sort_d_all$ix], kernel_type=kernel_type)-tilde_nu*A_t_Y[sort_d_all$ix]
        R_A_t_Y=R_A_t_Y_temp[d_sort_rev_ix]
        #R_A_t_Y-R%*%t(A_full)%*%Y[,i_l]
        A_R_A_t_Y=A_times_x_particle(output=R_A_t_Y, A_all_v=A_all_vec, num_neighbors_vec=2*num_neighbors_all_vec, D=D_y, N=N)
        #A_R_A_t_Y-A_full%*%R%*%t(A_full)%*%Y[,i_l]
        Y[,i_l]=A_R_A_t_Y*tau
      }
    }
    QR_Y=qr(Y)
    Q=qr.Q(QR_Y)
    T_mat=t(Q)%*%SS%*%Q
    #QR_Y$rank
    L_T=t(chol(T_mat+diag(l)))
    approx_log_det=2*sum(log(diag(L_T)))
    
    
    
    log_det_rec[1,i_sim]=log_det_rec[1,i_sim]+true_log_det
    log_det_rec[2,i_sim]=log_det_rec[2,i_sim]+approx_log_det
    
    
    
    m_CG=IKF_CG_particle(param=log(c(beta,tau)), kernel_type=kernel_type, delta_x_all=delta_x_all, output=output_all, 
                         A_all_v = A_all_vec, sort_d_all_ix=sort_d_all$ix,  num_neighbors_vec=2*num_neighbors_all_vec, tilde_nu=tilde_nu,
                         D=D_y,  N=N,   tol=tol,  maxIte = maxIte)
    
    ans_CG=m_CG[[1]]
    
    
    sigma_2_0_est = output_all%*%m_CG[[1]]/length(output_all)
    sigma_2_0_est_direct=output_all%*%solve(SS+diag(N),output_all)/length(output_all)
    
    true_log_lik=-true_log_det/2-N/2*log(2*pi)-N/2*log(sigma_2_0_est_direct)-1/2*N
    approx_log_lik=-approx_log_det/2-N/2*log(2*pi)-N/2*log(sigma_2_0_est)-1/2*N
    
    log_lik_rec[1,i_sim]=log_lik_rec[1,i_sim]+true_log_lik
    log_lik_rec[2,i_sim]=log_lik_rec[2,i_sim]+approx_log_lik
    
  }
  
  log_det_rec[,i_sim]=log_det_rec[,i_sim]/n_repeat
  log_lik_rec[,i_sim]=log_lik_rec[,i_sim]/n_repeat
}


#save.image(file='unnormalized_Vicsek_loglik.RData')

#load('unnormalized_Vicsek_loglik.RData')

dat1=data.frame(N=rep(T_sim_seq*n_t,3),
                method=factor(rep(c('Direct computation','IKF-CG','Difference'),each=length(T_sim_seq)),levels = c('Direct computation','IKF-CG','Difference')),
                log_lik=c(log_lik_rec[1,],log_lik_rec[2,],log_lik_rec[1,]-log_lik_rec[2,]))
pdf('plots/log_lik_unnormalized_Vicsek.pdf',height=3.5,width=5)
ggplot(dat1,aes(N,log_lik,shape=method,fill=method,size=method))+geom_point()+
  xlab(expression(tilde(N)))+ylab('log-likelihood')+
  scale_shape_manual(values=c(24,21,23))+
  scale_fill_manual(values=c("#fbb4ae","#80b1d3", "#ffffb3"))+
  scale_size_manual(values=c(3,2.5,2.8))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = "top", #legend.position = c(.3, .8),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.text=element_text(size=15),
        #legend.box.background = element_rect(),legend.box.margin = margin(0,0,0,0),
        axis.text=element_text(size=15),axis.title=element_text(size=20),
        legend.title=element_blank())
dev.off()


