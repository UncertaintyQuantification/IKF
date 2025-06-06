#### functions copied from the FastGaSP package ####
IKF = function(beta, tilde_nu, delta_x, output, kernel_type=kernel_type){
  if(kernel_type=='matern_5_2'){
    lambda=sqrt(5)*beta
    
    W0=Construct_W0_matern_5_2(1,lambda)  
    GG=Construct_G_matern_5_2(delta_x,lambda) 
    W=Construct_W_matern_5_2(1.0,delta_x,lambda,W0)
  }else if(kernel_type=='exp'){
    lambda=beta
    W0=Construct_W0_exp(1,lambda)
    GG=Construct_G_exp(delta_x,lambda) 
    W=Construct_W_exp(1,delta_x,lambda,W0)
  }
  Q_K=Get_Q_K(GG,W,W0,tilde_nu)
  res=Get_R_y(GG,Q_K[[1]],Q_K[[2]],output)#-tilde_nu*output
  
  return(res)
}

get_boundary_grid=function(px_min,px_max,py_min,py_max,nx,ny){
  grid_boundary_mat=matrix(NA,nx*ny,4)
  colnames(grid_boundary_mat)=c('pos_x_min','pos_x_max','pos_y_min','pos_y_max')
  
  ##get slightly larger grid 
  len_x_ori=px_max-px_min
  len_y_ori=py_max-py_min
  delta_x=len_x_ori/nx*0.1
  delta_y=len_y_ori/ny*0.1
  
  px_seq=seq(px_min-delta_x,px_max+delta_x,(px_max-px_min+2*delta_x)/(nx))
  py_seq=seq(py_min-delta_y,py_max+delta_y,(py_max-py_min+2*delta_y)/(ny))
  grid_boundary_mat[,1]=rep(px_seq[1:nx],ny)
  grid_boundary_mat[,2]=rep(px_seq[2:(nx+1)],ny)
  grid_boundary_mat[,3]=as.numeric(t(matrix(py_seq[1:ny],ny,nx)))
  grid_boundary_mat[,4]=as.numeric(t(matrix(py_seq[2:(ny+1)],ny,nx)))
  
  my_grid=list()
  
  
  Lx_min = min(grid_boundary_mat[,1:2]);
  Lx_max = max(grid_boundary_mat[,1:2]);
  Ly_min = min(grid_boundary_mat[,3:4]);
  Ly_max = max(grid_boundary_mat[,3:4]);
  
  len_x=  (max(grid_boundary_mat[,1:2])-min(grid_boundary_mat[,1:2]))/nx
  len_y=  (max(grid_boundary_mat[,3:4])-min(grid_boundary_mat[,3:4]))/ny
  
  
  grid_boundary_info=list()
  grid_boundary_info$grid_boundary_mat=grid_boundary_mat
  grid_boundary_info$grid_info=as.matrix(c(Lx_min,Lx_max,Ly_min,Ly_max,nx,ny,len_x,len_y),8,1)
  rownames(  grid_boundary_info$grid_info)=c('Lx_min','Lx_max','Ly_min','Ly_max','nx','ny','len_x','len_y')
  
  return(grid_boundary_info)
  
}

initiate_grid=function(grid_boundary_info){
  #include_theta = !is.null(theta) #if provided theta, then update theta to grid
  
  # Lx_min=unname(grid_boundary_info$grid_info['Lx_min',])
  # Lx_max=unname(grid_boundary_info$grid_info['Lx_max',])
  # Ly_min=unname(grid_boundary_info$grid_info['Ly_min',])
  # Ly_max=unname(grid_boundary_info$grid_info['Ly_max',])
  nx=unname(grid_boundary_info$grid_info['nx',])
  ny=unname(grid_boundary_info$grid_info['ny',])
  # len_x=unname(grid_boundary_info$grid_info['len_x',])
  # len_y=unname(grid_boundary_info$grid_info['len_y',])
  
  
  ##form neighboring particle
  ##from x first
  neighbor_index_list=as.list(1:(nx*ny))
  ##exterior
  for(i in 1:(nx*ny)){
    i_x=(i%%nx)
    if(i_x==0){
      i_x=nx
    }
    i_y=ceiling(i/nx)
    
    
    if((i_x-1)>0&(i_y-1)>0){
      neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y-2)*nx+(i_x-1))
    }
    
    if((i_y-1)>0){
      neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y-2)*nx+(i_x))
    }
    
    if((i_x+1)<=nx&(i_y-1)>0){
      neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y-2)*nx+(i_x+1))
    }
    
    if((i_x-1)>0){
      neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y-1)*nx+(i_x-1))
    }
    
    if((i_x+1)<=nx){
      neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y-1)*nx+(i_x+1))
    }
    
    if((i_x-1)>0&(i_y+1)<=ny){
      neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y)*nx+(i_x-1))
    }
    
    if((i_y+1)<=ny){
      neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y)*nx+(i_x))
    }
    if((i_x+1)<=nx&(i_y+1)<=ny){
      neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y)*nx+(i_x+1))
    }
  }
  
  
  
  
  #return(list(m_grid = m_grid, neighbor_index_list = neighbor_index_list))
  return(neighbor_index_list)
}

create_particle_grid=function(grid_boundary_info, pos_x,pos_y,vel_x,vel_y,neighbor_index_list){ #,return_theta=FALSE
  Lx_min=unname(grid_boundary_info$grid_info['Lx_min',])
  Lx_max=unname(grid_boundary_info$grid_info['Lx_max',])
  Ly_min=unname(grid_boundary_info$grid_info['Ly_min',])
  Ly_max=unname(grid_boundary_info$grid_info['Ly_max',])
  nx=unname(grid_boundary_info$grid_info['nx',])
  ny=unname(grid_boundary_info$grid_info['ny',])
  len_x=unname(grid_boundary_info$grid_info['len_x',])
  len_y=unname(grid_boundary_info$grid_info['len_y',])
  
  m_grid=vector(mode = 'list', nx*ny) #as.list(rep(NA,nx*ny))
  n_t=length(pos_x)
  
  for(i in 1:(nx*ny)){
    m_grid[[i]]=list(particle_pos = NULL, neighbor_pos=NULL, particle_vel=NULL, neighbor_vel=NULL)
    # if(return_theta){
    #   m_grid[[i]]=list(particle_pos = NULL, neighbor_pos=NULL, particle_vel=NULL, neighbor_vel=NULL,
    #                    particle_theta=NULL, neighbor_theta=NULL)
    # }else{
    #   
    # }
  }
  
  
  for(i in 1:n_t){
    i_x=ceiling((pos_x[i]-Lx_min)/len_x)
    i_y=ceiling((pos_y[i]-Ly_min)/len_y)
    
    index_grid=(i_y-1)*nx+i_x
    
    ##update pos and vel
    m_grid[[index_grid]]$particle_pos=cbind((m_grid[[index_grid]]$particle_pos),c(pos_x[i],pos_y[i]))
    m_grid[[index_grid]]$particle_vel=cbind((m_grid[[index_grid]]$particle_vel),c(vel_x[i],vel_y[i]))
    
    #if(return_theta) m_grid[[index_grid]]$particle_theta = c(m_grid[[index_grid]]$particle_theta,atan2(vel_y[i],vel_x[i]))
    
  }
  
  for(i in 1:(nx*ny)){
    #print(i)
    if(!is.null(m_grid[[i]]$particle_pos)){# only care about the grid that has particles
      neighbor_index = neighbor_index_list[[i]]
      for(idx in neighbor_index){
        if(!is.null(m_grid[[idx]]$particle_pos)){
          m_grid[[i]]$neighbor_pos=cbind(m_grid[[i]]$neighbor_pos,m_grid[[idx]]$particle_pos)
          m_grid[[i]]$neighbor_vel=cbind(m_grid[[i]]$neighbor_vel,m_grid[[idx]]$particle_vel)
          
          #if(return_theta) m_grid[[i]]$neighbor_theta = c(m_grid[[i]]$neighbor_theta, m_grid[[idx]]$particle_theta)
        }
      }
    }
    
  }
  
  return(m_grid)
  
}

find_grid_neighbors=function(pos_x_list,pos_y_list, vel_x_list,vel_y_list, time_range, #n_t,T_time,
                             grid_boundary_info){ #
  
  
  # Lx_min=grid_boundary_info$grid_info[1]
  # Lx_max=grid_boundary_info$grid_info[2]
  # Ly_min=grid_boundary_info$grid_info[3]
  # Ly_max=grid_boundary_info$grid_info[4]
  # nx=grid_boundary_info$grid_info[5]
  # ny=grid_boundary_info$grid_info[6]
  # len_x=grid_boundary_info$grid_info[7]
  # len_y=grid_boundary_info$grid_info[8]
  
  if(min(time_range)<1 | max(time_range)>length(pos_x_list)){
    stop('invalid time_range')
  }
  
  
  neighbor_index_list=initiate_grid(grid_boundary_info=grid_boundary_info)
  
  neighbors_info = vector(length(time_range), mode = "list")
  names(neighbors_info) <- paste0("time", time_range)
  for(t in time_range){
    #print(t)
    
    m_grid_here=create_particle_grid(grid_boundary_info=grid_boundary_info,
                                     pos_x=pos_x_list[[t]],pos_y=pos_y_list[[t]],
                                     vel_x=vel_x_list[[t]],vel_y=vel_y_list[[t]],
                                     neighbor_index_list=neighbor_index_list)
    
    neighbors_info[[t]] = m_grid_here
    
    
  }
  
  
  return(neighbors_info)
}

form_neighbors_unnormalized_Vicsek_with_r=function(threshold_r,pos_x_list,pos_y_list,vel_x_list,vel_y_list,
                                                   time_range,grid_boundary_info,neighbors_info,D_y){
  
  T_total = length(time_range)
  n_t = length(pos_x_list[[1]])
  
  Lx_min=unname(grid_boundary_info$grid_info['Lx_min',])
  Lx_max=unname(grid_boundary_info$grid_info['Lx_max',])
  Ly_min=unname(grid_boundary_info$grid_info['Ly_min',])
  Ly_max=unname(grid_boundary_info$grid_info['Ly_max',])
  nx=unname(grid_boundary_info$grid_info['nx',])
  ny=unname(grid_boundary_info$grid_info['ny',])
  len_x=unname(grid_boundary_info$grid_info['len_x',])
  len_y=unname(grid_boundary_info$grid_info['len_y',])
  
  
  A_vec=rep(NA,n_t*T_total*15*D_y)
  v_vec=rep(NA,n_t*T_total*15*D_y)
  num_neighbors_vec=rep(NA,n_t*T_total) ##here N neighbor each particle (obs) has
  
  count=0
  
  for(i_t in 1:length(time_range)){
    t = time_range[i_t]
    pos_x_t=pos_x_list[[t]]
    pos_y_t=pos_y_list[[t]]
    vel_x_t=vel_x_list[[t]]
    vel_y_t=vel_y_list[[t]]
    m_grid_here = neighbors_info[[t]]
    
    for(i in 1:n_t){
      input_pos_i=as.vector(c(pos_x_t[i],pos_y_t[i]))
      
      i_x=ceiling((input_pos_i[1]-Lx_min)/len_x)
      i_y=ceiling((input_pos_i[2]-Ly_min)/len_y)
      
      index_grid=(i_y-1)*nx+i_x
      d_vec_here_all=input_pos_i-as.matrix(m_grid_here[[index_grid]]$neighbor_pos)
      
      d_here=sqrt(colSums(d_vec_here_all^2))
      
      index_neighbor=which(d_here<threshold_r)
      # if(apolar_vicsek==F){
      #   index_neighbor=which(d_here<threshold_r)
      # }else{
      #   index_neighbor=which(d_here<threshold_r)
      #   index_same_v_direction=which(colSums(m_grid_here[[index_grid]]$neighbor_vel*input_vel_i)>=0)
      #   index_neighbor=intersect(index_neighbor,index_same_v_direction)
      # }
      n_neighbor=length(index_neighbor)
      num_neighbors_vec[n_t*(i_t-1)+i]=n_neighbor
      
      v_vec[count*D_y+(1:(n_neighbor*D_y))] = c(m_grid_here[[index_grid]]$neighbor_vel[1,index_neighbor],
                                                m_grid_here[[index_grid]]$neighbor_vel[2,index_neighbor])
      A_vec[2*count*D_y+(1:(2*n_neighbor*D_y))]=c(rep(c(1/n_neighbor,0),n_neighbor),rep(c(0,1/n_neighbor),n_neighbor))
      
      
      
      
      count=count+n_neighbor
      
      
    }
  }
  
  
  A_vec=A_vec[1:(2*count*D_y)]
  v_vec=v_vec[1:(count*D_y)]
  
  return(list(A_vec=A_vec, v_vec=v_vec, num_neighbors_vec=num_neighbors_vec))
}


####### New Functions ########
# using predicted interaction to predict trajectory for unnormalized Vicsek model
unnormalized_Vicsek_pred_trajec = function(est_param,v1,p1,v_abs,n_t,T_test,h,sigma_0,weights,v_train_vec,kernel_type = 'matern_5_2',tilde_nu=0.1,noise_type='Gaussian'){
  beta=est_param[1]
  tau=est_param[2]
  threshold_r=est_param[3]
  
  pos=matrix(NA,2*n_t,T_test+1)
  v=matrix(NA,2*n_t,T_test+1)
  
  v[,1]=v1
  pos[,1]=p1
  
  for(t in 1:T_test){
    input_here=matrix(pos[,t],2,n_t)
    v_vec_here=matrix(v[,t],2,n_t)
    
    ##augment it for fast prediction 
    d_aug=c(as.vector(v_vec_here),(v_train_vec))
    d_aug_sort=sort(d_aug,index.return=T)
    d_aug_sort_x=d_aug_sort$x
    d_aug_sort_rev_ix=sort(d_aug_sort$ix,index.return=T)$ix ###this is to reverse the previous sort 
    
    delta_x_aug=d_aug_sort_x[2:length(d_aug_sort_x)]-d_aug_sort_x[1:(length(d_aug_sort_x)-1)]
    
    w_aug=c(rep(0,n_t*2),weights)
    
    
    # param_tilde=log(c(beta,tilde_nu)) 
    # pred_mean_aug=R_times_z(param_tilde, have_noise=T, delta_x=delta_x_aug, z=w_aug[d_aug_sort$ix],
    #                         kernel_type=kernel_type)-tilde_nu*w_aug[d_aug_sort$ix]
    pred_mean_aug=IKF(beta, tilde_nu, delta_x_aug, w_aug[d_aug_sort$ix], kernel_type=kernel_type)-tilde_nu*w_aug[d_aug_sort$ix]
    
    pred_mean_fast=matrix(pred_mean_aug[d_aug_sort_rev_ix][1:(2*n_t)],2,n_t)
    
    
    
    for(i in 1:n_t){
      input_i=input_here[,i]#as.vector(pos[(i-1)*D+1:D,t])
      d_vec_here_all=(input_i-input_here) ###this is negative d for A
      #d_vec_here_all=(input_i-input_here) ### d_vec for A
      
      d_here=sqrt(colSums(d_vec_here_all^2))
      #v_here=sqrt(colSums(v_vec_here^2))
      
      index_neighbor=which(d_here<threshold_r) #&v_here>0
      
      v[(i-1)*2+1,t+1]=mean(pred_mean_fast[1,index_neighbor])
      v[(i-1)*2+2,t+1]=mean(pred_mean_fast[2,index_neighbor])
      
      ##add noise
      if(noise_type=='Gaussian'){
        v[(i-1)*2+1,t+1]=  v[(i-1)*2+1,t+1]+sigma_0*rnorm(1)  ##
        v[(i-1)*2+2,t+1]= v[(i-1)*2+2,t+1]+sigma_0*rnorm(1)  ##
      }else if(noise_type=='Uniform'){
        v[(i-1)*2+1,t+1]=  v[(i-1)*2+1,t+1]+sigma_0*(runif(1)-0.5)  ##
        v[(i-1)*2+2,t+1]= v[(i-1)*2+2,t+1]+sigma_0*(runif(1)-0.5)
      }
      
      pos[(i-1)*2+1:2,t+1]=pos[(i-1)*2+1:2,t]+v[(i-1)*2+1:2,t]*h
    }
  }
  
  return(list(pos_pred=pos,v_pred=v))
}

# residual bootstrap for unnormalized Vicsek model
residual_bootstrap_unnormalized_Vicsek = function(est, sim, B = 100, 
                                                  param_ini, cut_r_max,nx=NULL, ny=NULL,
                                                  kernel_type = 'matern_5_2', 
                                                  tilde_nu=0.1,tol = 1e-6,maxIte=1000) {
  
  T_sim = sim@T_time
  D_y = sim@D_y
  n_t = sim@n_particles
  N=n_t*T_sim*D_y
  
  # Initialize arrays to store results
  beta_est_rec = numeric(3)
  names(beta_est_rec) = c("est","LB","UB")
  tau_est_rec = numeric(3)
  names(tau_est_rec) = c("est","LB","UB")
  radius_est_rec = numeric(3)
  names(radius_est_rec) = c("est","LB","UB")
  
  beta=est@parameters['beta']
  tau=est@parameters['tau']
  radius=est@parameters['radius']
  
  beta_est_rec[1]=beta
  tau_est_rec[1]=tau
  radius_est_rec[1]=radius
  
  v_all_vec = est@training_data$training_velocity
  A_all_vec = est@training_data$A_v
  num_neighbors_all_vec = est@training_data$num_neighbors
  
  sort_v_all=sort(v_all_vec,index.return=T)
  N_tilde_all=length(v_all_vec)
  delta_x_all=sort_v_all$x[-1]-sort_v_all$x[-N_tilde_all]
  
  d_sort_rev_ix=sort(sort_v_all$ix,index.return=T)$ix ###this is to reverse the previous sort 
  
  
  output_all=as.vector(rbind(
    unlist(sim@vx_list[1+1:T_sim]),
    unlist(sim@vy_list[1+1:T_sim])
  ))
  
  
  m_CG=IKF_CG_particle(param=log(c(beta,tau)), kernel_type=kernel_type, delta_x_all=delta_x_all, output=output_all, 
                       A_all_v = A_all_vec, sort_d_all_ix=sort_v_all$ix,  num_neighbors_vec=2*num_neighbors_all_vec, tilde_nu=tilde_nu,
                       D=D_y,  N=N,   tol=tol,  maxIte = maxIte)
  
  ans_CG=m_CG[[1]] 
  
  ans_CG_tilde=ans_CG*tau
  
  w_CG=A_t_times_x_particle(output=ans_CG_tilde, A_all_v=A_all_vec,  num_neighbors_vec=2*num_neighbors_all_vec,  
                            D_y=D_y, N_tilde=N_tilde_all)
  
  # param_tilde=log(c(beta,tilde_nu)) 
  # 
  # pred_mean_rev=R_times_z(param_tilde, have_noise=T, delta_x=delta_x_all, z=w_CG[sort_theta_all$ix],
  #                         kernel_type=kernel_type)-tilde_nu*w_CG[sort_theta_all$ix]
  
  pred_mean_rev = IKF(beta, tilde_nu, delta_x=delta_x_all, output=w_CG[sort_v_all$ix], kernel_type=kernel_type)-tilde_nu*w_CG[sort_v_all$ix]
  
  pred_mean_fast = pred_mean_rev[d_sort_rev_ix]
  
  pred_mean_output=A_times_x_particle( pred_mean_fast,  A_all_vec,  2*num_neighbors_all_vec,  
                                       D_y,  N)
  
  res=output_all-pred_mean_output
  
  beta_rec=rep(NA,B)
  tau_rec=rep(NA,B)
  radius_rec=rep(NA,B)
  
  for(b in 1:B){
    cat(sprintf("\rDoing bootstrap %d of %d", b, B))
    flush.console()
    
    output_all_new=pred_mean_output+sample(res,replace = T)
    
    # sim_new = sim
    # sim_new@theta_list = split(c(sim@theta_list[[1]], output_all_new), rep(1:((T_sim+1)), each=n_t))
    
    T_index_time = 1:T_sim
    T_index_ho=seq(5,T_sim,5) ##every 5 use the last one as holdout
    T_index_train=(1:T_sim)[-T_index_ho]
    
    output_new= output_all_new[unlist(lapply(T_index_train, function(t) {
      ((t-1)*n_t*D_y + 1):(t*n_t*D_y)
    }))]
    ho_output_new=output_all_new[unlist(lapply(T_index_ho, function(t) {
      ((t-1)*n_t*D_y + 1):(t*n_t*D_y)
    }))]
    
    #length(as.vector(rbind(unlist(sim@vx_list[T_index_train + 1]),unlist(sim@vy_list[T_index_train + 1]))))
    
    est_new = particle_interaction_est_unnormalized_Vicsek(sim, param=param_ini, cut_r_max=cut_r_max, nx=nx, ny=ny,
                                                           kernel_type=kernel_type,tol = tol,
                                                           output = output_new, ho_output = ho_output_new)
    
    
    beta_rec[b]=est_new@parameters['beta']
    tau_rec[b]=est_new@parameters['tau']
    radius_rec[b]=est_new@parameters['radius']
    
  }
  cat("\n")
  
  beta_est_rec[2]=quantile(beta_rec,0.025)
  beta_est_rec[3]=quantile(beta_rec,0.975)
  
  tau_est_rec[2]=quantile(tau_rec,0.025)
  tau_est_rec[3]=quantile(tau_rec,0.975)
  
  radius_est_rec[2]=quantile(radius_rec,0.025)
  radius_est_rec[3]=quantile(radius_rec,0.975)
  
  # Return results
  return(list(
    beta_est_rec = beta_est_rec,
    tau_est_rec = tau_est_rec,
    radius_est_rec = radius_est_rec
  ))
}

model_pred_fixed_param_IKF = function(sim,param,cut_r_max,nx=NULL, ny=NULL,
                                      kernel_type, tol, maxIte=1000,tilde_nu=0.1, testing_input){ 
  
  beta=exp(param[1])
  tau=exp(param[2])  ###sigma_2/sigma_2_0
  threshold_r=exp(param[3])/(1+exp(param[3])) * cut_r_max
  
  
  px_list = sim@px_list
  py_list = sim@py_list
  vx_list = sim@vx_list
  vy_list = sim@vy_list
  n_t = sim@n_particles
  T_sim = sim@T_time
  D_y = sim@D_y
  
  N=n_t*T_sim*D_y 
  
  T_index_time = 1:T_sim
  
  px_min=min(unlist(px_list))
  px_max=max(unlist(px_list))
  py_min=min(unlist(py_list))
  py_max=max(unlist(py_list))
  
  if(is.null(nx)){
    nx=floor((px_max-px_min)/cut_r_max)
  }else{
    if(cut_r_max>(px_max-px_min)/nx) nx=floor((px_max-px_min)/cut_r_max)
  }
  
  if(is.null(ny)){
    ny=floor((py_max-py_min)/cut_r_max)
  }else{
    if(cut_r_max>(py_max-py_min)/ny) ny=floor((py_max-py_min)/cut_r_max)
  }
  
  grid_boundary_info=get_boundary_grid(px_min=px_min,px_max=px_max,
                                       py_min=py_min,py_max=py_max,nx=nx,ny=ny)
  #print(grid_boundary_info)
  
  
  neighbors_info = find_grid_neighbors(pos_x_list=px_list,pos_y_list=py_list,
                                       vel_x_list=vx_list,vel_y_list=vy_list, 
                                       time_range=T_index_time, grid_boundary_info=grid_boundary_info)
  
  
  
  output_all=as.vector(rbind(
    unlist(vx_list[1+T_index_time]),
    unlist(vy_list[1+T_index_time])
  ))
  
  
  time_rec=system.time({
    ans_neighbors_all=form_neighbors_unnormalized_Vicsek_with_r(threshold_r=threshold_r,pos_x_list=px_list,pos_y_list=py_list,
                                                                vel_x_list=vx_list,vel_y_list=vy_list,time_range=T_index_time,
                                                                grid_boundary_info=grid_boundary_info,neighbors_info=neighbors_info,D_y=D_y)
    
    
    A_all_vec=ans_neighbors_all$A_vec
    v_all_vec=ans_neighbors_all$v_vec
    num_neighbors_all_vec=ans_neighbors_all$num_neighbors_vec
    sort_v_all=sort(v_all_vec,index.return=T)
    N_tilde_all=length(v_all_vec) ###this is N_j in the paper
    
    
    delta_x_all=sort_v_all$x[-1]-sort_v_all$x[-N_tilde_all]
    
    
    m_CG=IKF_CG_particle(param=log(c(beta,tau)), kernel_type=kernel_type, delta_x_all=delta_x_all, output=output_all, 
                         A_all_v = A_all_vec, sort_d_all_ix=sort_v_all$ix,  num_neighbors_vec=2*num_neighbors_all_vec, tilde_nu=tilde_nu,
                         D=D_y,  N=N,   tol=tol,  maxIte = maxIte)
    
    ans_CG=m_CG[[1]] 
    
    ans_CG_tilde=ans_CG*tau
    
    sigma_2_0_est = output_all%*%ans_CG/length(output_all) ##sometimes negative? solved
    
    
    w_CG=A_t_times_x_particle(output=ans_CG_tilde, A_all_v=A_all_vec,  num_neighbors_vec=2*num_neighbors_all_vec,  
                              D_y=D_y, N_tilde=N_tilde_all)
    
    
    testing_n = length(testing_input)
    
    
    sigma_2_est=sigma_2_0_est*tau
    
    param=log(c(beta,tau))
    
    d_aug=c(testing_input,(v_all_vec))
    d_aug_sort=sort(d_aug,index.return=T)
    d_aug_sort_x=d_aug_sort$x
    d_aug_sort_rev_ix=sort(d_aug_sort$ix,index.return=T)$ix ###this is to reverse the previous sort 
    
    delta_x_aug=d_aug_sort_x[2:length(d_aug_sort_x)]-d_aug_sort_x[1:(length(d_aug_sort_x)-1)]
    
    
    w_aug=c(rep(0,testing_n),w_CG)
    
    pred_mean_aug = IKF(beta=beta, tilde_nu=tilde_nu, 
                        delta_x=delta_x_aug, output=w_aug[d_aug_sort$ix], kernel_type=kernel_type)-tilde_nu*w_aug[d_aug_sort$ix]
    
    pred_mean_fast=pred_mean_aug[d_aug_sort_rev_ix][1:testing_n]
    
  })
  
  return(list(pred_mean_fast=pred_mean_fast,time_rec=time_rec))
  
  
}


model_pred_fixed_param_CG = function(sim,param,cut_r_max,nx=NULL, ny=NULL,
                                     kernel_type, tol, maxIte=1000,testing_input){
  
  beta=exp(param[1])
  tau=exp(param[2])  ###sigma_2/sigma_2_0
  threshold_r=exp(param[3])/(1+exp(param[3])) * cut_r_max
  
  
  px_list = sim@px_list
  py_list = sim@py_list
  vx_list = sim@vx_list
  vy_list = sim@vy_list
  n_t = sim@n_particles
  T_sim = sim@T_time
  D_y = sim@D_y
  
  N=n_t*T_sim*D_y 
  
  
  T_index_time = 1:T_sim
  
  
  px_min=min(unlist(px_list))
  px_max=max(unlist(px_list))
  py_min=min(unlist(py_list))
  py_max=max(unlist(py_list))
  
  if(is.null(nx)){
    nx=floor((px_max-px_min)/cut_r_max)
  }else{
    if(cut_r_max>(px_max-px_min)/nx) nx=floor((px_max-px_min)/cut_r_max)
  }
  
  if(is.null(ny)){
    ny=floor((py_max-py_min)/cut_r_max)
  }else{
    if(cut_r_max>(py_max-py_min)/ny) ny=floor((py_max-py_min)/cut_r_max)
  }
  
  grid_boundary_info=get_boundary_grid(px_min=px_min,px_max=px_max,
                                       py_min=py_min,py_max=py_max,nx=nx,ny=ny)
  #print(grid_boundary_info)
  
  
  neighbors_info = find_grid_neighbors(pos_x_list=px_list,pos_y_list=py_list,
                                       vel_x_list=vx_list,vel_y_list=vy_list, 
                                       time_range=T_index_time, grid_boundary_info=grid_boundary_info)
  
  output_all=as.vector(rbind(
    unlist(vx_list[1+T_index_time]),
    unlist(vy_list[1+T_index_time])
  ))
  
  time_rec=system.time({
    ans_neighbors_all=form_neighbors_unnormalized_Vicsek_with_r(threshold_r=threshold_r,pos_x_list=px_list,pos_y_list=py_list,
                                                                vel_x_list=vx_list,vel_y_list=vy_list,time_range=T_index_time,
                                                                grid_boundary_info=grid_boundary_info,neighbors_info=neighbors_info,D_y=D_y)
    
    A_all_vec=ans_neighbors_all$A_vec
    v_all_vec=ans_neighbors_all$v_vec
    num_neighbors_all_vec=ans_neighbors_all$num_neighbors_vec
    
    N_tilde_all=length(v_all_vec) 
    
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
    
    
    nu=1/tau
    res=CG_direct(A=A_full, Sigma=R, b=output_all,nu=nu, tol=tol, maxIte=maxIte)
    R_y_inv_y=res[[1]]
    
    A_full_t_R_y_inv_y=t(A_full)%*%R_y_inv_y
    
    
    r0=abs(outer(testing_input,(v_all_vec),'-'))
    if(kernel_type=='exp'){
      r = exp(-beta*r0)
    }else if(kernel_type=='matern_5_2'){
      r = matern_5_2_funct(r0, beta)
    }
    pred_mean_CG=r%*%A_full_t_R_y_inv_y
  })
  
  return(list(pred_mean_CG=pred_mean_CG,time_rec=time_rec))
  
}




model_pred_fixed_param_direct = function(sim,param,cut_r_max,nx=NULL, ny=NULL,
                                         kernel_type, testing_input){
  beta=exp(param[1])
  tau=exp(param[2])  ###sigma_2/sigma_2_0
  threshold_r=exp(param[3])/(1+exp(param[3])) * cut_r_max
  
  
  px_list = sim@px_list
  py_list = sim@py_list
  vx_list = sim@vx_list
  vy_list = sim@vy_list
  n_t = sim@n_particles
  T_sim = sim@T_time
  D_y = sim@D_y
  
  N=n_t*T_sim*D_y 
  
  
  T_index_time = 1:T_sim
  
  
  px_min=min(unlist(px_list))
  px_max=max(unlist(px_list))
  py_min=min(unlist(py_list))
  py_max=max(unlist(py_list))
  
  if(is.null(nx)){
    nx=floor((px_max-px_min)/cut_r_max)
  }else{
    if(cut_r_max>(px_max-px_min)/nx) nx=floor((px_max-px_min)/cut_r_max)
  }
  
  if(is.null(ny)){
    ny=floor((py_max-py_min)/cut_r_max)
  }else{
    if(cut_r_max>(py_max-py_min)/ny) ny=floor((py_max-py_min)/cut_r_max)
  }
  
  grid_boundary_info=get_boundary_grid(px_min=px_min,px_max=px_max,
                                       py_min=py_min,py_max=py_max,nx=nx,ny=ny)
  #print(grid_boundary_info)
  
  
  neighbors_info = find_grid_neighbors(pos_x_list=px_list,pos_y_list=py_list,
                                       vel_x_list=vx_list,vel_y_list=vy_list, 
                                       time_range=T_index_time, grid_boundary_info=grid_boundary_info)
  
  output_all=as.vector(rbind(
    unlist(vx_list[1+T_index_time]),
    unlist(vy_list[1+T_index_time])
  ))
  
  time_rec=system.time({
    ans_neighbors_all=form_neighbors_unnormalized_Vicsek_with_r(threshold_r=threshold_r,pos_x_list=px_list,pos_y_list=py_list,
                                                                vel_x_list=vx_list,vel_y_list=vy_list,time_range=T_index_time,
                                                                grid_boundary_info=grid_boundary_info,neighbors_info=neighbors_info,D_y=D_y)
    
    A_all_vec=ans_neighbors_all$A_vec
    v_all_vec=ans_neighbors_all$v_vec
    num_neighbors_all_vec=ans_neighbors_all$num_neighbors_vec
    
    N_tilde_all=length(v_all_vec) 
    
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
    
    nu=1/tau
    R_y=A_full%*%R%*%t(A_full)+nu*diag(N)
    R_y_inv=solve(R_y)
    R_y_inv_y=R_y_inv%*%output_all
    
    A_full_t_R_y_inv_y=t(A_full)%*%R_y_inv_y
    
    
    r0=abs(outer(testing_input,(v_all_vec),'-'))
    if(kernel_type=='exp'){
      r = exp(-beta*r0)
    }else if(kernel_type=='matern_5_2'){
      r = matern_5_2_funct(r0, beta)
    }
    pred_mean_direct=r%*%A_full_t_R_y_inv_y
    
  })
  return(list(pred_mean_direct=pred_mean_direct,time_rec=time_rec))
  
}




# perform residual bootstrap in cell data
residual_bootstrap_with_prediction_cell <- function(exp_data, est, B, n_record, output_all, sigma_2_record, 
                                                    param_ini, cut_r_max, nx=NULL, ny=NULL,
                                                    kernel_type, tilde_nu=0.1, tol=1e-6, maxIte=1000, 
                                                    testing_inputs, apolar_vicsek = FALSE, direction) {
  D_y = exp_data@D_y
  N=D_y*sum(n_record)
  T_time = exp_data@T_time
  
  beta=est@parameters['beta']
  tau=est@parameters['tau']
  radius=est@parameters['radius']
  
  v_all_vec = est@training_data$training_velocity
  A_v_all_vec = est@training_data$A_v
  num_neighbors_all_vec = est@training_data$num_neighbors
  
  sort_v_all=sort(v_all_vec,index.return=T)
  N_tilde_all=length(v_all_vec)
  delta_x_all=sort_v_all$x[-1]-sort_v_all$x[-N_tilde_all]
  
  d_sort_rev_ix=sort(sort_v_all$ix,index.return=T)$ix ###this is to reverse the previous sort 
  
  m_CG=IKF_CG_particle_cell(param=log(c(beta,tau)), kernel_type=kernel_type, delta_x_all=delta_x_all, output=output_all, 
                            A_all_v=A_v_all_vec, sort_d_all_ix=sort_v_all$ix, 
                            sigma_2_vec=sigma_2_record, num_neighbors_vec=num_neighbors_all_vec, tilde_nu=tilde_nu, 
                            D=D_y, n_t_record=n_record, tol = tol, maxIte = maxIte)
  
  ans_CG=m_CG[[1]] 
  
  ans_CG_tilde=ans_CG*tau
  
  w_CG=A_t_times_x_particle(output=ans_CG_tilde, A_all_v=A_v_all_vec,  num_neighbors_vec=num_neighbors_all_vec,  
                            D_y=D_y, N_tilde=N_tilde_all)
  
  pred_mean_rev = IKF(beta, tilde_nu, delta_x=delta_x_all, output=w_CG[sort_v_all$ix], kernel_type=kernel_type)-tilde_nu*w_CG[sort_v_all$ix]
  
  pred_mean_fast = pred_mean_rev[d_sort_rev_ix]
  
  pred_mean_output=A_times_x_particle(pred_mean_fast, A_v_all_vec, num_neighbors_all_vec,  
                                      D_y, N)
  
  res=output_all-pred_mean_output
  
  boot_est_param = matrix(NA, nrow = B, ncol = 3)
  colnames(boot_est_param) = c("beta", "tau", "radius")
  
  testing_n = nrow(testing_inputs)
  boot_output = array(NA, c(3, testing_n, B), dimnames = list(c("est","LB","UB"), NULL,paste0('b',1:B)))
  
  for(b in 1:B){
    print(paste0("Doing bootstrap ", b))
    output_all_new=pred_mean_output+sample(res,replace = T)
    
    T_index_ho=seq(5,T_time,5) ##every 5 use the last one as holdout
    T_index_train=(1:T_time)[-T_index_ho]
    
    n_record_cum = cumsum(c(0, n_record))
    output_new = unlist(mapply(function(i) output_all_new[(n_record_cum[i]+1):n_record_cum[i+1]], T_index_train))
    ho_output_new = unlist(mapply(function(i) output_all_new[(n_record_cum[i]+1):n_record_cum[i+1]], T_index_ho))
    
    est_new = fit(exp_data,param=param_ini,cut_r_max=cut_r_max, kernel_type = kernel_type, tol = tol, 
                  output = output_new, ho_output = ho_output_new, testing_inputs=testing_inputs, #compute_CI = FALSE,
                  apolar_vicsek = apolar_vicsek, direction = direction)
    
    boot_est_param[b, ] = c(est_new@parameters['beta'], est_new@parameters['tau'], est_new@parameters['radius'])
    
    boot_output[1,,b] = est_new@predictions$mean
    boot_output[2,,b] = est_new@predictions$lower95
    boot_output[3,,b] = est_new@predictions$upper95
  }
  
  return(list(boot_est_param = boot_est_param, boot_output = boot_output))
}


find_grid_neighbors_with_pred_cell=function(px_list_test,py_list_test, vx_list_test,vy_list_test, #time_range, #n_t,T_time,
                                            grid_boundary_info,pred_v_list){ #
  
  T_test = length(px_list_test)
  
  Lx_min=unname(grid_boundary_info$grid_info['Lx_min',])
  Lx_max=unname(grid_boundary_info$grid_info['Lx_max',])
  Ly_min=unname(grid_boundary_info$grid_info['Ly_min',])
  Ly_max=unname(grid_boundary_info$grid_info['Ly_max',])
  nx=unname(grid_boundary_info$grid_info['nx',])
  ny=unname(grid_boundary_info$grid_info['ny',])
  len_x=unname(grid_boundary_info$grid_info['len_x',])
  len_y=unname(grid_boundary_info$grid_info['len_y',])
  
  
  neighbor_index_list=initiate_grid(grid_boundary_info=grid_boundary_info)
  
  
  neighbors_info_with_pred = vector(T_test, mode = "list")
  names(neighbors_info_with_pred) = paste0("time", 1:T_test)
  for(t in 1:T_test){
    pos_x=px_list_test[[t]]
    pos_y=py_list_test[[t]]
    vel_x=vx_list_test[[t]]
    vel_y=vy_list_test[[t]]
    pred=pred_v_list[[t]]
    
    
    m_grid=vector(mode = 'list', nx*ny) #as.list(rep(NA,nx*ny))
    n_t=length(pos_x)
    
    for(i in 1:(nx*ny)){
      m_grid[[i]]=list(neighbor_pos=NULL, neighbor_vel=NULL, neighbor_pred = NULL)
    }
    
    
    for(i in 1:n_t){
      i_x=ceiling((pos_x[i]-Lx_min)/len_x)
      i_y=ceiling((pos_y[i]-Ly_min)/len_y)
      
      index_grid=(i_y-1)*nx+i_x
      
      neighbor_index = neighbor_index_list[[index_grid]]
      
      for(idx in neighbor_index){
        m_grid[[idx]]$neighbor_pos = cbind(m_grid[[idx]]$neighbor_pos, c(pos_x[i],pos_y[i]))
        m_grid[[idx]]$neighbor_vel = cbind(m_grid[[idx]]$neighbor_vel, c(vel_x[i],vel_y[i]))
        m_grid[[idx]]$neighbor_pred = c(m_grid[[idx]]$neighbor_pred, pred[i])
      }
      
    }
    neighbors_info_with_pred[[t]] = m_grid
    
  }
  
  return(neighbors_info_with_pred)
}

pred_cell_approx_var = function(exp_train, exp_test, est, cut_r_max, testing_n = 200, 
                                nx = NULL, ny = NULL, tol, tilde_nu = 0.1, maxIte = 1000,
                                apolar_vicsek, direction, compute_var_approx = TRUE){
  
  # train data extraction (used for variance only)
  vx_list = get_consecutive_data(exp_train, "vx")$start
  vy_list = get_consecutive_data(exp_train, "vy")$start
  
  n_record = sapply(vx_list,length)
  
  # test data extraction
  px_list_test = get_consecutive_data(exp_test, "px")$start
  py_list_test = get_consecutive_data(exp_test, "py")$start
  
  vx_list_test = get_consecutive_data(exp_test, "vx")$start
  vy_list_test = get_consecutive_data(exp_test, "vy")$start
  
  
  
  n_record_test = sapply(px_list_test,length)
  T_test = exp_test@T_time
  D_y = exp_test@D_y
  
  if(direction == "x"){
    test_data_input = unlist(vx_list_test)
    test_output = unlist(get_consecutive_data(exp_test, "vx")$end)
    if(compute_var_approx){
      testing_v_input = seq(min(unlist(vx_list)),max(unlist(vx_list)),length.out=testing_n)
      sigma_2_record=sapply(vx_list, var)
      sigma_2_record_test=sapply(vx_list_test, var)
    } 
  }else if(direction == "y"){
    test_data_input = unlist(vy_list_test)
    test_output = unlist(get_consecutive_data(exp_test, "vy")$end)
    if(compute_var_approx){
      testing_v_input = seq(min(unlist(vy_list)),max(unlist(vy_list)),length.out=testing_n)
      sigma_2_record=sapply(vy_list, var)
      sigma_2_record_test=sapply(vy_list_test, var)
    }
  }
  
  # get predictive mean for each input data
  param = log(c(est@parameters[1], est@parameters[2], est@parameters[3]/(cut_r_max-est@parameters[3])))
  est_test = fit(exp_train,param=param,cut_r_max=cut_r_max,est_param = FALSE,
                 tol = tol,testing_inputs = matrix(test_data_input), compute_CI= FALSE,
                 apolar_vicsek = apolar_vicsek, direction = direction)
  
  
  # convert mean as list
  pred_mean = est_test@predictions$mean
  pred_mean_list=list()
  count_sum=0
  for(t in 1:T_test){
    pred_mean_list[[t]]=(pred_mean[count_sum+1:n_record_test[t]])
    count_sum=count_sum+n_record_test[t]
  }
  
  
  px_min=min(unlist(px_list_test))
  px_max=max(unlist(px_list_test))
  py_min=min(unlist(py_list_test))
  py_max=max(unlist(py_list_test))
  
  
  if(is.null(nx)) nx=floor((px_max-px_min)/cut_r_max)
  if(is.null(ny)) ny=floor((py_max-py_min)/cut_r_max)
  
  grid_boundary_info = get_boundary_grid(px_min,px_max,py_min,py_max,nx,ny)
  
  Lx_min=unname(grid_boundary_info$grid_info['Lx_min',])
  Lx_max=unname(grid_boundary_info$grid_info['Lx_max',])
  Ly_min=unname(grid_boundary_info$grid_info['Ly_min',])
  Ly_max=unname(grid_boundary_info$grid_info['Ly_max',])
  nx=unname(grid_boundary_info$grid_info['nx',])
  ny=unname(grid_boundary_info$grid_info['ny',])
  len_x=unname(grid_boundary_info$grid_info['len_x',])
  len_y=unname(grid_boundary_info$grid_info['len_y',])
  
  
  
  neighbors_info_with_pred = find_grid_neighbors_with_pred_cell(px_list_test=px_list_test,py_list_test=py_list_test,
                                                                vx_list_test=vx_list_test,vy_list_test=vy_list_test, 
                                                                grid_boundary_info=grid_boundary_info,pred_v_list=pred_mean_list)
  
  
  num_neighbors_vec=rep(NA,sum(n_record_test))
  v_neighbor_pred=rep(NA,sum(n_record_test))
  
  beta = unname(est@parameters['beta'])
  tau = unname(est@parameters['tau'])
  threshold_r = unname(est@parameters['radius'])
  
  if(compute_var_approx){
    pred_var=rep(NA,sum(n_record_test))
    
    num_neighbors_all_vec = est_test@training_data$num_neighbors
    A_v_all_vec = est_test@training_data$A_v
    v_all_vec = est_test@training_data$training_velocity
    
    sigma_2_0_prop_est = est_test@sigma_2_0_est
    sigma_2_est=sigma_2_0_prop_est*tau
    
    sort_v_all=sort(v_all_vec,index.return=T)
    N_tilde_all=length(v_all_vec) ###this is N_j in the paper
    delta_x_all=sort_v_all$x[-1]-sort_v_all$x[-N_tilde_all]
    
    
    N=length(num_neighbors_all_vec)
    A_r_v_rec=matrix(NA,N,testing_n) # used to record predictive variance of testing_v_input
    R_inv_A_r_v_rec=matrix(NA,N,testing_n)
    
    
    r0_v=abs(outer(testing_v_input,(v_all_vec),'-'))
    if(kernel_type=='exp'){
      r_v = exp(-beta*r0_v)
    }else if(kernel_type=='matern_5_2'){
      r_v = matern_5_2_funct(r0_v, beta)
    }
    
  }
  
  
  
  for(t in 1:T_test){
    print(t)
    
    px_t=px_list_test[[t]]
    py_t=py_list_test[[t]]
    vx_t=vx_list_test[[t]]
    vy_t=vy_list_test[[t]]
    m_grid_here = neighbors_info_with_pred[[t]]
    n_t = n_record_test[t]
    
    index_start = ifelse(t==1, 0, sum(n_record_test[1:(t-1)]))
    
    for(i in 1:n_t){
      input_pos_i=as.vector(c(px_t[i],py_t[i]))
      input_vel_i=as.vector(c(vx_t[i],vy_t[i]))
      
      i_x=ceiling((input_pos_i[1]-Lx_min)/len_x)
      i_y=ceiling((input_pos_i[2]-Ly_min)/len_y)
      
      index_grid=(i_y-1)*nx+i_x
      d_vec_here_all=input_pos_i-as.matrix(m_grid_here[[index_grid]]$neighbor_pos)
      
      d_here=sqrt(colSums(d_vec_here_all^2))
      if(apolar_vicsek){
        index_neighbor=which(d_here<threshold_r)
        index_same_v_direction=which(colSums(m_grid_here[[index_grid]]$neighbor_vel*input_vel_i)>=0)
        index_neighbor=intersect(index_neighbor,index_same_v_direction)
      }else{
        index_neighbor=which(d_here<threshold_r)
      }
      
      #proposed model
      v_neighbor_pred[index_start+i]=mean(m_grid_here[[index_grid]]$neighbor_pred[index_neighbor])
      num_neighbors_vec[index_start+i]=length(index_neighbor)
      
      
      
      
      
      if(compute_var_approx){
        if(direction=="x"){
          neighbors_velocity=m_grid_here[[index_grid]]$neighbor_vel[1,index_neighbor]
        }else if(direction=="y"){
          neighbors_velocity=m_grid_here[[index_grid]]$neighbor_vel[2,index_neighbor]
        }
        
        cov_approx=matrix(NA,length(index_neighbor),length(index_neighbor))
        d_approx_ind=rep(NA,length(index_neighbor))
        for(p in 1:length(index_neighbor)){
          d_approx_ind[p]=which.min(abs(neighbors_velocity[p]-testing_v_input))
        }
        
        
        for(p_i in 1:length(index_neighbor)){
          for(p_j in 1:length(index_neighbor)){
            if(is.na(A_r_v_rec[1,d_approx_ind[p_i]])){
              A_r_v_i=A_times_x_particle(output=r_v[d_approx_ind[p_i],], A_all_v=A_v_all_vec,  num_neighbors_vec=num_neighbors_all_vec,
                                         D=D_y, N)
              A_r_v_rec[,d_approx_ind[p_i]]=A_r_v_i
              
              
              tol_interval=tol*10^{-14}
              R_inv_r_v_all=IKF_CG_particle_cell(param=log(c(beta,tau)), kernel_type=kernel_type, delta_x_all=delta_x_all, output=A_r_v_i,
                                                 A_all_v=A_v_all_vec, sort_d_all_ix=sort_v_all$ix,
                                                 sigma_2_vec=sigma_2_record, num_neighbors_vec=num_neighbors_all_vec, tilde_nu=tilde_nu,
                                                 D=D_y, n_t_record=n_record, tol = tol_interval, maxIte = maxIte)
              
              R_inv_A_r_v_rec[,d_approx_ind[p_i]]=R_inv_r_v_all[[1]]
              
            }
            
            if(is.na(A_r_v_rec[1,d_approx_ind[p_j]])){
              A_r_v_i=A_times_x_particle(output=r_v[d_approx_ind[p_j],], A_all_v=A_v_all_vec,  num_neighbors_vec=num_neighbors_all_vec,
                                         D=D_y, N)
              A_r_v_rec[,d_approx_ind[p_j]]=A_r_v_i
              
              
              tol_interval=tol*10^{-14}
              R_inv_r_v_all=IKF_CG_particle_cell(param=log(c(beta,tau)), kernel_type=kernel_type, delta_x_all=delta_x_all, output=A_r_v_i,
                                                 A_all_v=A_v_all_vec, sort_d_all_ix=sort_v_all$ix,
                                                 sigma_2_vec=sigma_2_record, num_neighbors_vec=num_neighbors_all_vec, tilde_nu=tilde_nu,
                                                 D=D_y, n_t_record=n_record, tol = tol_interval, maxIte = maxIte)
              
              R_inv_A_r_v_rec[,d_approx_ind[p_j]]=R_inv_r_v_all[[1]]
              
            }
            
            
            
            cov_approx[p_i,p_j]=t(A_r_v_rec[,d_approx_ind[p_i]])%*%R_inv_A_r_v_rec[,d_approx_ind[p_j]]*tau
            
            r0_star_star=abs(outer(neighbors_velocity,(neighbors_velocity),'-'))
            if(kernel_type=='exp'){
              r_star_star = exp(-beta*r0_star_star)
            }else if(kernel_type=='matern_5_2'){
              r_star_star = matern_5_2_funct(r0_star_star, beta)
            }
            
            c_v_star=r_star_star-cov_approx
            c_y_star=sum(c_v_star)/(length(index_neighbor)^2)*sigma_2_est+sigma_2_0_prop_est*sigma_2_record_test[t]
            pred_var[index_start+i]=c_y_star
            
          }
        }
        
      }
      
    }
    
  }
  res = list(v_neighbor_pred=v_neighbor_pred,num_neighbors_vec=num_neighbors_vec,test_output = test_output)
  
  if(compute_var_approx){
    res$pred_var = pred_var
    # res$A_r_v_rec = A_r_v_rec
    # res$R_inv_A_r_v_rec = R_inv_A_r_v_rec
  }
  return(res)
}






