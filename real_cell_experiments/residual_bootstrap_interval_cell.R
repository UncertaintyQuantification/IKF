library(FastGaSP)
source('../functions/functions_particle.R')

particle_data = read.csv('cell_data.csv')

exp = trajectory_data(particle_data)

exp_sub=extract_time_window(exp, 1, 53) # use first 53 time frames to build the model

p_bar=2 ##expected number of neighbor
tol=10^{-8}*sum(sapply(exp_sub@particle_tracking,nrow))*p_bar

apolar_vicsek=T 
kernel_type='matern_5_2'


testing_n=200
param_rec=matrix(NA,2,3,dimnames = list(c("x","y"),c("beta","tau","r")))
testing_input_rec=matrix(NA,2,testing_n,dimnames = list(c("x","y")))
testing_pred_rec=array(NA,c(3,testing_n,2),dimnames = list(c("est","LB","UB"),NULL,c("x","y")))
fit_time = rep(NA, 2)

#bootstrap
do_bootstrap = FALSE
if(do_bootstrap){
  B=100
  param_boot = array(NA, c(B,3,2), dimnames = list(paste0('b',1:B), c("beta","tau","r"), c("x","y")))
  testing_pred_boot=array(NA, c(3,testing_n,B,2), dimnames = list(c("est","LB","UB"), NULL,paste0('b',1:B), c("x","y")))
}

set.seed(0)
for(i_direction in 1:2){
  print(i_direction)
  
  if(i_direction == 1){
    direction = "x"
  }else if(i_direction == 2){
    direction = "y"
  }
  
  if(direction == "x"){
    vx_pairs = get_consecutive_data(exp_sub, "vx")
    vx_list = vx_pairs$start
    vx_end_list = vx_pairs$end
    output_all=unlist(vx_end_list)

    sigma_2_record=sapply(vx_list, var)
    n_record = sapply(vx_list,length)
    
    testing_input = seq(min(unlist(vx_list)),max(unlist(vx_list)),length.out=testing_n)
    testing_inputs = matrix(as.numeric(testing_input))
  }else if(direction == "y"){
    vy_pairs = get_consecutive_data(exp_sub, "vy")
    vy_list = vy_pairs$start
    vy_end_list = vy_pairs$end
    output_all=unlist(vy_end_list)
    
    sigma_2_record=sapply(vy_list, var)
    n_record = sapply(vy_list,length)
    
    testing_input = seq(min(unlist(vy_list)),max(unlist(vy_list)),length.out=testing_n)
    testing_inputs = matrix(as.numeric(testing_input))
  }
  
  testing_input_rec[i_direction,]=testing_input
  
  cut_r_max = 20
  
  param_ini=log(c(1,10,10/(cut_r_max-10)))
  
  fit_time[i_direction]=system.time({
  est = fit(exp_sub,param_ini,cut_r_max=cut_r_max, #nx = 30, ny = 20,
            tol = tol, testing_inputs = testing_inputs, apolar_vicsek = apolar_vicsek, direction = direction)
  })[3]
  
  param_rec[i_direction,] = est@parameters
  testing_pred_rec[1,,i_direction] = est@predictions$mean
  testing_pred_rec[2,,i_direction] = est@predictions$lower95
  testing_pred_rec[3,,i_direction] = est@predictions$upper95
  
  if(do_bootstrap){
    rb = residual_bootstrap_with_prediction_cell(exp_data=exp_sub, est=est, B=B, 
                                                 n_record = n_record, output_all=output_all, sigma_2_record=sigma_2_record, 
                                                 param_ini=param_ini, cut_r_max=cut_r_max, 
                                                 kernel_type=kernel_type, tol=tol, 
                                                 testing_inputs=testing_inputs, apolar_vicsek = apolar_vicsek, direction=direction) 
    
    param_boot[,,i_direction] = rb$boot_est_param
    testing_pred_boot[,,,i_direction] = rb$boot_output
  }

}

if(do_bootstrap){
  pred_sd_boot = array(NA, c(testing_n,B,2), dimnames = list(NULL,paste0('b',1:B), c("x","y")))
  standard_q_seq=matrix(qnorm(seq(0.01,0.99,0.01),0,1),1)
  pred_CI_diff_q_boot = array(NA, c(testing_n,length(standard_q_seq),B,2), dimnames = list(NULL,paste0('q',seq(0.01,0.99,0.01)),paste0('b',1:B), c("x","y")))
  
  for(i_direction in 1:2){
    for(b in 1:B){
      est_sd = (testing_pred_boot['UB',,b,i_direction]-testing_pred_boot['LB',,b,i_direction])/2/qnorm(0.975)
      pred_sd_boot[,b,i_direction] = est_sd
      pred_CI_diff_q_boot[,,b,i_direction] = est_sd %*% standard_q_seq+matrix(testing_pred_boot['est',,b,i_direction],testing_n,99)
    }
  }
}

# > fit_time_max_100
# [1] 1171.223 1369.215
# > fit_time # cut_r_max = 50
# [1]  989.621 1147.115
# > fit_time # cut_r_max = 20
# [1] 905.479 839.456


#####save.image(file=paste0('cell_res_bootstrap_B_',B,'.RData0'))
load(paste0('cell_res_bootstrap_B_',100,'.RData'))



#save.image(file=paste0('cell_res_no_bootstrap.RData'))

################################ plot ################################
for(i_direction in 1:2){
  if(i_direction == 1){
    direction = "x"
  }else if(i_direction == 2){
    direction = "y"
  }
  if(direction=="x"){
    #hist_info=hist(vx_list[[T_time]],breaks = 30)
    hist_info=hist(unlist(vx_list),breaks = 30)
    main_text = "Horizontal direction"
  }else if(direction=="y"){
    #hist_info=hist(vy_list[[T_time]],breaks = 30)
    hist_info=hist(unlist(vy_list),breaks = 40)
    main_text = "Vertical direction"
  }
  testing_v_input=testing_input_rec[i_direction,]
  #LB95_mean=apply(testing_output_boot[[i_direction]], 2, function(x) quantile(x,0.025))
  #UB95_mean=apply(testing_output_boot[[i_direction]], 2, function(x) quantile(x,0.975))
  LB95_mean=rep(NA,testing_n)
  UB95_mean=rep(NA,testing_n)
  for(i in 1:testing_n){
    LB95_mean[i]=quantile(pred_CI_diff_q_boot[i,,,i_direction],0.025)
    UB95_mean[i]=quantile(pred_CI_diff_q_boot[i,,,i_direction],0.975)
  }
  
  
  ylim_kernel=range(min(-0.01,LB95_mean),max(0.01,UB95_mean))
  #ylim_kernel=range(min(-0.01,testing_pred_rec[,,i_direction]),max(0.01,testing_pred_rec[,,i_direction]))
  ylim_hist=c(0,max(hist_info$counts))
  b = diff(ylim_hist)/diff(ylim_kernel)
  a = ylim_hist[1] - b*ylim_kernel[1]
  
  jpeg(paste("../plots/kernel_bootstrap_interval_",direction,".jpeg",sep=""), quality = 100,width = 1000, height = 550, res = 250)
  par(mgp=c(2,1,0), mar=c(3,3.6,1,4.6)+.1)
  plot(hist_info,border=F,yaxt="n",ylab="",xlab="d",col="#ffeefc",xlim=c(-.02,.02),main=main_text)
  polygon(c(testing_v_input,rev(testing_v_input)), a+c(LB95_mean,rev(UB95_mean))*b, col = "grey80", border = "grey80")
  #polygon(c(testing_v_input,rev(testing_v_input)), a+c(testing_pred_rec[2,,i_direction],rev(testing_pred_rec[3,,i_direction]))*b, col = "grey80", border = F)
  lines(testing_v_input,a+testing_pred_rec[1,,i_direction]*b,type='l',col='#00618e',lty=2,lwd=1.5)
  # without bootstrap
  # lines(testing_v_input,a+testing_pred_rec[2,,i_direction]*b,type='l',lty=3,lwd=1)
  # lines(testing_v_input,a+testing_pred_rec[3,,i_direction]*b,type='l',lty=3,lwd=1)
  axis(2, at=a+seq(-.03,.03,.005)*b,labels=NA, las=2,col="#00618e")
  axis(2, at=a+seq(-.03,.03,.01)*b,labels=seq(-.03,.03,.01), las=2,col="#00618e")
  count_label_seq = seq(-6000,(ceiling(ylim_hist[2]/6000)+1)*6000,6000)
  axis(4, at=count_label_seq,labels=count_label_seq, las=2,col="pink")
  mtext("z(d)", side = 2, line = 2.6)  
  mtext("cell count", side = 4, line = 3.5)  
  #add lines on top and bottom
  axis(1, at=c(-1,1),labels=c(-1,1),las=0)
  axis(3, at=c(-1,1),labels=c(-1,1),las=0)
  abline(a,b,col="purple")
  dev.off()
  
}




