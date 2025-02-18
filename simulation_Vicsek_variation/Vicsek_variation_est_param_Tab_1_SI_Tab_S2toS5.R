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


##noise 
noise_type='Gaussian'


nrmse_v_rec = covered95_v_rec = length95_v_rec = 
  matrix(NA,dim(settings)[1],n_repeat,
         dimnames = list(paste("T",settings[,1],",n",settings[,2],",s",settings[,3],sep=""),
                         paste("Exp",1:n_repeat,sep = "")))

nrmse_f_rec = covered95_f_rec = length95_f_rec = 
  matrix(NA,dim(settings)[1],n_repeat,
         dimnames = list(paste("T",settings[,1],",n",settings[,2],",s",settings[,3],sep=""),
                         paste("Exp",1:n_repeat,sep = "")))

beta_v_rec = beta_f_rec = tau_v_rec = tau_f_rec = radius_rec = 
  matrix(NA,dim(settings)[1],n_repeat,
         dimnames = list(paste("T",settings[,1],",n",settings[,2],",s",settings[,3],sep=""),
                         paste("Exp",1:n_repeat,sep = "")))


for(i_s in 1:dim(settings)[1]){
  print(i_s)
  T_sim = settings[i_s,1]
  n_t = settings[i_s,2]
  sigma_0 = settings[i_s,3]
  
  D_y=2 ##here y obs is 2 D
  N=n_t*T_sim*D_y ##this one is the number of observations 
  cut_r=.5 ##truth cutoff
  
  for(i_repeat in 1:n_repeat){
    set.seed((i_s-1)*n_repeat+i_repeat+1)
    
    # simulation
    vx_abs=0.5
    vy_abs=0.5
    v_abs=sqrt(vx_abs^2+vy_abs^2)
    
    sim = simulate_particle(v_abs=v_abs, n_t = n_t, T_sim = T_sim, 
                              cut_r = cut_r, sigma_0 = sigma_0, model = 'two_interactions_Vicsek')
    
    # fit
    cut_r_max=1.5
    
    p_bar=2 ##expected number of neighbor
    tol=10^{-8}*N*p_bar ##tolerance

    
    testing_n=200
    testing_v_input=seq(min(c(unlist(sim@vx_list), unlist(sim@vy_list))),max(c(unlist(sim@vx_list), unlist(sim@vy_list))),length.out=testing_n)
    testing_v_output=testing_v_input 
    
    testing_d_input=seq(0.025,sim@radius,length.out=testing_n)
    testing_f_output=f_Vicsek_variation(testing_d_input,r_max=sim@radius)
    
    testing_inputs = cbind(testing_v_input, testing_d_input)
    
    if(kernel_type=='matern_5_2'){
      param_ini=log(c(0.01,10,10^4,10^4,0.8/(cut_r_max-0.8)))
      
      est = fit(sim,param=param_ini,cut_r_max=cut_r_max, 
                kernel_type=kernel_type,tol=tol,testing_inputs = testing_inputs)
      
      
      if(est@parameters[5]<0.1 | est@parameters[5]>0.9){
        param_ini=log(c(0.01,10,10^4,10^4,0.3/(cut_r_max-0.3)))
        est = fit(sim,param=param_ini,cut_r_max=cut_r_max, 
                  kernel_type=kernel_type,tol=tol,testing_inputs = testing_inputs)
        
      }
    }else if(kernel_type=='exp'){
      param_ini=log(c(0.01,10,10^2,10^2,0.8/(cut_r_max-0.8))) 
      
      est = fit(sim,param=param_ini,cut_r_max=cut_r_max, 
                kernel_type=kernel_type,tol=tol,testing_inputs = testing_inputs)
      
      if(est@parameters[5]<0.1 | est@parameters[5]>0.9){
        param_ini=log(c(0.01,10,10^2,10^2,0.3/(cut_r_max-0.3))) 
        est = fit(sim,param=param_ini,cut_r_max=cut_r_max, 
                  kernel_type=kernel_type,tol=tol,testing_inputs = testing_inputs)
        
      }
      
    }
    
    
    pred_mean_v = est@predictions$mean_v
    lower95_v = est@predictions$lower95_v
    upper95_v = est@predictions$upper95_v
    
    pred_mean_f = est@predictions$mean_f
    lower95_f = est@predictions$lower95_f
    upper95_f = est@predictions$upper95_f
    
    NRMSE_v = mean( (pred_mean_v-testing_v_output)^2)/sd(testing_v_output)
    coverage95_v = mean(testing_v_output<upper95_v & testing_v_output>lower95_v)
    length95_v = mean(upper95_v-lower95_v)
    
    NRMSE_f = mean( (pred_mean_f-testing_f_output)^2)/sd(testing_f_output)
    coverage95_f = mean(testing_f_output<upper95_f & testing_f_output>lower95_f)
    length95_f = mean(upper95_f-lower95_f)
    
    
    nrmse_v_rec[i_s,i_repeat]=NRMSE_v
    nrmse_f_rec[i_s,i_repeat]=NRMSE_f
    covered95_v_rec[i_s,i_repeat]=coverage95_v
    covered95_f_rec[i_s,i_repeat]=coverage95_f
    length95_v_rec[i_s,i_repeat]=length95_v
    length95_f_rec[i_s,i_repeat]=length95_f
    
    beta_v_rec[i_s,i_repeat]=est@parameters['beta_v']
    beta_f_rec[i_s,i_repeat]=est@parameters['beta_f']
    tau_v_rec[i_s,i_repeat]=est@parameters['tau_v']
    tau_f_rec[i_s,i_repeat]=est@parameters['tau_f']
    radius_rec[i_s,i_repeat]=est@parameters['radius']
    
    print(c(nrmse_v_rec[i_s,i_repeat],nrmse_f_rec[i_s,i_repeat],
            covered95_v_rec[i_s,i_repeat], covered95_f_rec[i_s,i_repeat],
            length95_v_rec[i_s,i_repeat], length95_f_rec[i_s,i_repeat],
            radius_rec[i_s,i_repeat]))
    
  }
  
  print(c(mean(nrmse_v_rec[i_s,]),mean(nrmse_f_rec[i_s,]),
          mean(covered95_v_rec[i_s,]),mean(covered95_f_rec[i_s,]),
          mean(length95_v_rec[i_s,]),mean(length95_f_rec[i_s,])))
  
  
}


# save.image('record_est_param_Vicsek_variation_matern.RData')
# save.image('record_est_param_Vicsek_variation_exp.RData')


load("record_est_param_Vicsek_variation_matern_5_2_Feb_18.RData")
#load("record_est_param_Vicsek_variation_exp_Feb_18.RData")


apply(nrmse_v_rec, 1, mean)
apply(nrmse_f_rec, 1, mean)

apply(covered95_v_rec, 1, mean)
apply(covered95_f_rec, 1, mean)

apply(length95_v_rec, 1, mean)
apply(length95_f_rec, 1, mean)






