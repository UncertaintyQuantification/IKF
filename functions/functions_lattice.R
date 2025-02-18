branin <- function(xx, a=1, b=5.1/(4*pi^2), c=5/pi, r=6, s=10, t=1/(8*pi))
{
  ##########################################################################
  #
  # BRANIN FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2)  #x1 in [-5, 10], x2 in [0, 15]
  # a = constant (optional), with default value 1
  # b = constant (optional), with default value 5.1/(4*pi^2)
  # c = constant (optional), with default value 5/pi
  # r = constant (optional), with default value 6
  # s = constant (optional), with default value 10
  # t = constant (optional), with default value 1/(8*pi)
  #
  ##########################################################################
  
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- a * (x2 - b*x1^2 + c*x1 - r)^2
  term2 <- s*(1-t)*cos(x1)
  
  y <- term1 + term2 + s
  return(y)
}


CV_missing_zero_mean = function(param,kernel_type,input1,input2,
                                Y_train,Y_valid,miss_ind,valid_ind,
                                tilde_nu,tol,maxIte,print_par = T){
  #param=log(c(beta1,beta2,tau))
  beta1 = exp(param[1])
  beta2 = exp(param[2])
  tau = exp(param[3])
  
  n1 = length(input1)
  n2 = length(input2)
  
  delta_x1=input1[-1]-input1[-n1]
  delta_x2=input2[-1]-input2[-n2]
  
  # #original CG
  # m_CG=fast_pred_sparse_CG_missing(param, kernel_type, delta_x1, delta_x2, Y_train, n1,
  #                                  n2, miss_ind, tilde_nu, tol = tol, maxIte = maxIte)
  #CG with specified initial x0
  m_CG=fast_pred_sparse_CG_missing_with_ini(param, kernel_type, delta_x1, delta_x2, Y_train, n1,
                                            n2, miss_ind, x0, tilde_nu, tol = tol, maxIte = maxIte)
  x0 <<- m_CG[[1]] # Update global x0
  
  CG_ans=m_CG[[1]]
  CG_ans_tilde=tau*CG_ans
  
  A_t_CG = A_t_times_x_sparse_missing(CG_ans_tilde,miss_ind,N)
  pred_Z_CG = R_times_z_kronecker(A_t_CG,param,n1,n2,delta_x1,delta_x2,tilde_nu,kernel_type)
  
  rmse=sqrt(mean((pred_Z_CG[valid_ind]-Y_valid)^2))
  #rmse
  #sqrt(mean((pred_Z_CG[valid_ind]-Z_full[valid_ind])^2))
  
  # image2D(matrix(Y_obs_NA,n1,n2),x=input1,y=input2,main='observed output',xlab='x1',ylab='x2')
  # image2D(matrix(pred_Z_CG,n1,n2),x=input1,y=input2,main='predicted mean',xlab='x1',ylab='x2')
  
  if(print_par){print(c(m_CG[[3]],beta1,beta2,tau,rmse))}
  return(log(rmse))
}




branin_est_param = function(Y_obs_NA, obs_ind,inside_ind, N, N0, valid_prop = 0.2, 
                            input1, input2, kernel_type, param_ini, 
                            tilde_nu, CGtol, CGmaxIte,maxit=100,print_par){
  
  n = round(N0 * (1-valid_prop))
  train_ind = sort(sample(obs_ind,n))
  
  miss_ind = setdiff(1:N,train_ind)
  valid_ind = setdiff(miss_ind,inside_ind)
  
  Y_train = Y_obs_NA[train_ind]
  Y_valid = Y_obs_NA[valid_ind]
  
  x0<<-rep(0,n)
  
  m = optim(param_ini,CV_missing_zero_mean,control=list(maxit=maxit),#maxit=200
            kernel_type=kernel_type,input1=input1,input2=input2,
            Y_train=Y_train,Y_valid=Y_valid,miss_ind=miss_ind,valid_ind=valid_ind,
            tilde_nu=tilde_nu,tol=CGtol,maxIte=CGmaxIte,print_par = print_par)
  
  param=m$par
  #print(c(param,m$value))
  return(param)
}


missing_prediction = function(param, Y_obs, input1, input2, non_obs_ind, theta = 0,
                              kernel_type, tilde_nu, tol, maxIte){
  n1 = length(input1)
  n2 = length(input2)
  N = n1*n2
  N0 = length(Y_obs)
  
  beta1 = exp(param[1])
  beta2 = exp(param[2])
  tau = exp(param[3])
  
  
  delta_x1=input1[-1]-input1[-n1]
  delta_x2=input2[-1]-input2[-n2]
  
  A_X_theta = rep(theta, N0)
  
  m_CG=fast_pred_sparse_CG_missing(param, kernel_type, delta_x1, delta_x2, Y_obs-A_X_theta, n1, 
                                   n2, non_obs_ind, tilde_nu, tol = tol, maxIte = maxIte)
  
  
  CG_ans=m_CG[[1]]
  CG_ans_tilde=tau*CG_ans
  
  A_t_CG = A_t_times_x_sparse_missing(CG_ans_tilde,non_obs_ind,N)
  pred_Z_CG = rep(theta,N)+R_times_z_kronecker(A_t_CG,param,n1,n2,delta_x1,delta_x2,tilde_nu,kernel_type)
  
  return(pred_Z_CG)
}


satellite_est_param = function(param_ini,valid_prop=0.1,split_method,
                               N,N0,full_input,full_output,x1,x2,obs_ind,kernel_type,zero_mean = TRUE,
                               tilde_nu, CGtol,CGmaxIte,maxit=100,print_par){
  n = round(N0 * (1-valid_prop)) ##use around 10% as cross-validation
  
  #n = round(N0 * (1-(0.1))) 
  
  if(split_method == 'random'){
    train_ind = sort(sample(obs_ind,n))
    non_train_ind = setdiff(1:N,train_ind)
    valid_ind = setdiff(obs_ind,train_ind)
    non_valid_ind = setdiff(1:N,valid_ind)
    
  }else if(split_method == 'disk'){
    ##need to contain signal so it won't be over-smooth
    #center = c(.3*min(x1)+.7*max(x1),
    #           .4*min(x2)+.6*max(x2))
    center = c(min(x1)+0.8*(max(x1)-min(x1)),
               min(x2)+0.8*(max(x2)-min(x1)))
    
    # center = c(min(x1)+0.8*(max(x1)-min(x1)),
    #            min(x2)+0.6*(max(x2)-min(x1)))
    
    #center = c(.5*min(x1)+.5*max(x1),
    #           .5*min(x2)+.5*max(x2))
    
    
    dist_x1 = mean(diff(unique(x1)))
    dist_x2 = mean(diff(unique(x2)))
    
    radius = sqrt((N0-n) / pi *dist_x1*dist_x2)
    distances = apply(full_input, 1, function(x) sqrt(sum((x-center)^2)))
    valid_ind = setdiff(which(distances < radius),non_obs_ind) 
    
    non_valid_ind = setdiff(1:N,valid_ind)
    train_ind = setdiff(obs_ind,valid_ind)
    non_train_ind = setdiff(1:N,train_ind)
    n = length(train_ind)
  }
  
  
  Y_train = full_output[train_ind]
  Y_valid = full_output[valid_ind]
  
  
  
  # # ##plot the training and testing
  # Y_obs_NA = full_output
  # Y_obs_NA[-obs_ind] = NA
  # Y_train_NA = Y_obs_NA
  # Y_train_NA[-train_ind] = NA
  # 
  # image2D(matrix(Y_obs_NA,n1,n2),x=x1,y=x2,main='observed output',xlab='x1',ylab='x2',zlim=zlim)
  # image2D(matrix(Y_train_NA,n1,n2),x=x1,y=x2,main='train output',xlab='x1',ylab='x2',zlim=zlim)
  
  if(zero_mean){
    x0<<-rep(0,n)
    m = optim(param_ini,CV_missing_zero_mean,control=list(maxit=maxit),
              kernel_type=kernel_type,input1=x1,input2=x2,
              Y_train=Y_train,Y_valid=Y_valid,miss_ind=non_train_ind,valid_ind=valid_ind,
              tilde_nu=tilde_nu,tol=CGtol,maxIte=CGmaxIte, print_par = print_par)
    return(m$par)
  }else{
    x0_R_y<<-rep(0,n)
    x0_R_A_X <<- rep(0,n)
    m = optim(param_ini,CV_missing_minimize_l2_mean,control=list(maxit=maxit),
              kernel_type=kernel_type,input1=x1,input2=x2,A_train_X=rep(1,n),a_valid_X = rep(1,N0-n),
              Y_train=Y_train,Y_valid=Y_valid,miss_ind=non_train_ind,valid_ind=valid_ind,non_valid_ind=non_valid_ind,
              tilde_nu=tilde_nu,tol=CGtol,CG_MaxIter=CGmaxIte,return_theta=F, print_par = print_par)
    
    x0_R_y<<-rep(0,n)
    x0_R_A_X <<- rep(0,n)
    theta = CV_missing_minimize_l2_mean(m$par,kernel_type,input1=x1,input2=x2,A_train_X=rep(1,n),a_valid_X = rep(1,N0-n),
                                        Y_train=Y_train,Y_valid=Y_valid,miss_ind=non_train_ind,valid_ind=valid_ind,non_valid_ind=non_valid_ind,
                                        tilde_nu=tilde_nu,tol=CGtol*10^(-5),CG_MaxIter=CGmaxIte*10,return_theta=T,print_par = print_par)
    
    return(list(param=m$par,theta = theta))
    
  }
  
  
}

