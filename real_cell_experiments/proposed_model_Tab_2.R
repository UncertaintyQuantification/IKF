library(FastGaSP)
library(RobustGaSP)
source('../functions/functions_particle.R')

particle_data = read.csv('cell_data.csv')

exp = trajectory_data(particle_data)
exp@T_time

exp_train=extract_time_window(exp, 1, 53) # use first 53 time frames to build the model
exp_test=extract_time_window(exp, 54, exp@T_time) # use the last 50 to test

p_bar=2 ##expected number of neighbor
tol=10^{-8}*sum(sapply(exp_train@particle_tracking,nrow))*p_bar

apolar_vicsek=T 
kernel_type='matern_5_2'


### residual bootstrap
param_rec=matrix(NA,2,3,dimnames = list(c("x","y"),c("beta","tau","r")))
pred_save = list()
pred_result=matrix(NA,2,3,dimnames = list(c("x","y"),c("RMSE","L95","P95")))


for(i_direction in 1:2){
  
  if(i_direction == 1){
    direction = "x"
  }else if(i_direction == 2){
    direction = "y"
  }
  
  #testing_input_rec[i_direction,]=testing_input
  #T_test = length(n_record_test)
  
  # estimate parameters
  cut_r_max = 20
  
  param_ini=log(c(1,10,10/(cut_r_max-10)))
  
  est = fit(exp_train,param = param_ini,cut_r_max=cut_r_max, #nx = 30, ny = 20,
            tol = tol, testing_inputs = NULL, apolar_vicsek = apolar_vicsek, direction = direction)
  
  
  param_rec[i_direction,] = est@parameters
  


  test_pred = pred_cell_approx_var(exp_train=exp_train, exp_test=exp_test, est=est, 
                                   cut_r_max=cut_r_max, testing_n = 200,tol=tol, 
                                   apolar_vicsek=apolar_vicsek, direction=direction, 
                                   compute_var_approx = TRUE)
  
  pred_save[[i_direction]] = test_pred
  rmse = sqrt(mean((test_pred$test_output-test_pred$v_neighbor_pred)^2))
  
  LB=test_pred$v_neighbor_pred+sqrt(abs(test_pred$pred_var))*qnorm(0.025)
  UB=test_pred$v_neighbor_pred+sqrt(abs(test_pred$pred_var))*qnorm(0.975)
  
  length95 = mean(UB-LB)
  coverage95 = mean(UB>test_pred$test_output & LB<test_pred$test_output)
  
  pred_result[i_direction,]=c(rmse, length95, coverage95)
  
}


save.image("Train53_Test50_res_max20.RData")



