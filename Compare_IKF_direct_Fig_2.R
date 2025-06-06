library(FastGaSP)
library(RobustGaSP)

kernel_type = 'matern_5_2'

#### 1. time comparison ####


M=20 # number of points
n_seq_vec=1000+((1:M))*1000

repeated_times=20 # number of repeated times

record_time_IKF=matrix(NA,M,repeated_times)
record_time_direct=matrix(NA,M,repeated_times)


for(i_M in 1:M){
  print(i_M)
  for(j_repeat in 1:repeated_times){
    set.seed(j_repeat+i_M*repeated_times)
    n=n_seq_vec[i_M]
    
    input=sort(runif(n))
    z=rnorm(n)
    
    beta=10
    
    record_time_direct[i_M,j_repeat]=system.time({
      R0=abs(outer(input,input,'-'))
      if(kernel_type=='matern_5_2'){
        R = matern_5_2_funct(R0,beta_i=beta)
      }else if(kernel_type=='exp'){
        R = exp(-beta*R0)
      }
      y_direct=R%*%z
    })[3]
    
    tilde_nu = 10^{-11}
    record_time_IKF[i_M,j_repeat]=system.time({
      delta_x = input[-1] - input[-n]   
      y_IKF = IKF(beta = beta, tilde_nu = tilde_nu, delta_x = delta_x, 
                  output = z, kernel_type = kernel_type) - tilde_nu * z
    })[3]
    
  }
}
  
max(abs(y_IKF-y_direct))

#### 2. time comparison (log) ####

M=20
n_seq_log_10_spaced=floor(10^{0.15*(1:M)+3})

repeated_times=10

record_time_IKF_log10_spaced=matrix(NA,M,repeated_times)
record_time_direct_log10_spaced=matrix(NA,M,repeated_times)

for(i_M in 1:M){
  print(i_M)
  for(j_repeat in 1:repeated_times){
    set.seed(j_repeat+i_M*repeated_times)
    n=n_seq_log_10_spaced[i_M]
    
    input=sort(runif(n))
    z=rnorm(n)
    
    beta=10
    
    if(n<=25000){
      record_time_direct_log10_spaced[i_M,j_repeat]=system.time({
        R0=abs(outer(input,input,'-'))
        if(kernel_type=='matern_5_2'){
          R = matern_5_2_funct(R0,beta_i=beta)
        }else if(kernel_type=='exp'){
          R = exp(-beta*R0)
        }
        y_direct=R%*%z
      })[3]
    }
    
    
    tilde_nu = 10^{-11}
    record_time_IKF_log10_spaced[i_M,j_repeat]=system.time({
      delta_x = input[-1] - input[-n]   
      y_IKF = IKF(beta = beta, tilde_nu = tilde_nu, delta_x = delta_x, 
                  output = z, kernel_type = kernel_type) - tilde_nu * z
    })[3]
    
  }
}
#max(abs(y_IKF-y_direct))

#### 3. robust & non-robust IKF ####
M=20
n_seq_test_error=1000+((1:M))*1000

repeated_times=10

max_error_IKF=matrix(NA,M,repeated_times)
max_error_non_robust_IKF=matrix(NA,M,repeated_times)

for(i_M in 1:M){
  print(i_M)
  for(j_repeat in 1:repeated_times){
    set.seed(j_repeat+i_M*repeated_times)
    n=n_seq_test_error[i_M]
    
    input=sort(runif(n))
    z=rnorm(n)
    
    beta=10
    
    # direct computation
    R0=abs(outer(input,input,'-'))
    if(kernel_type=='matern_5_2'){
      R = matern_5_2_funct(R0,beta_i=beta)
    }else if(kernel_type=='exp'){
      R = exp(-beta*R0)
    }
    y_direct=R%*%z
    
    # robust IKF
    delta_x = input[-1] - input[-n]   
    y_IKF = IKF(beta = beta, tilde_nu = tilde_nu, delta_x = delta_x, 
                output = z, kernel_type = kernel_type) - tilde_nu * z
    
    # non-robust IKF
    y_non_robust_IKF = IKF(beta = beta, tilde_nu = 0, delta_x = delta_x, 
                           output = z, kernel_type = kernel_type)
    
    max_error_IKF[i_M,j_repeat]=max(abs(y_IKF-y_direct))
    max_error_non_robust_IKF[i_M,j_repeat]=max(abs(y_non_robust_IKF-y_direct))
    
  }
}


#### 4. plot ####
library(tidyverse)
library(scales)


dat1=data.frame(N=rep(n_seq_vec,2),method=rep(c('Direct computation','IKF'),each=length(n_seq_vec)),
                time=c(rowMeans(record_time_direct),rowMeans(record_time_IKF)))

#pdf(file='time_comparison.pdf',height=3.5,width=5)
ggplot(dat1,aes(N,time))+geom_point(size = 3,aes(shape=method,fill=method))+
  xlab(expression(N))+ylab('time (s)')+
  scale_shape_manual(values=c(24,21))+
  scale_fill_manual(values=c("#fbb4ae","#80b1d3"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position ="top",#legend.position = c(.26, .86),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.text=element_text(size=15),
        #legend.box.background = element_rect(),legend.box.margin = margin(0,0,0,0),
        axis.text=element_text(size=15),axis.title=element_text(size=20),
        legend.title=element_blank())
#dev.off() 


dat2=data.frame(N=rep(n_seq_log_10_spaced,2),
                method=rep(c('Direct computation','IKF'),each=length(n_seq_log_10_spaced)),
                time=c(rowMeans(record_time_direct_log10_spaced),rowMeans(record_time_IKF_log10_spaced)))
#pdf(file='time_comparison_log10_spaced.pdf',height=3.5,width=5)
ggplot(dat2,aes(N,time,shape=method,fill=method))+geom_point(size = 3)+
  xlab(expression(N))+ylab('time (s)')+
  scale_shape_manual(values=c(24,21))+
  scale_fill_manual(values=c("#fbb4ae","#80b1d3"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position ="top",#legend.position = c(.72, .15),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.text=element_text(size=15),
        #legend.box.background = element_rect(),legend.box.margin = margin(0,0,0,0),
        axis.text=element_text(size=15),axis.title=element_text(size=20),
        legend.title=element_blank())+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
#dev.off()

dat3=data.frame(N=rep(n_seq_test_error,2),
                method=factor(rep(c('IKF','Non-robust IKF'),each=length(n_seq_test_error)),levels = c('Non-robust IKF','IKF')),
                error=c(rowMeans(max_error_IKF),rowMeans(max_error_non_robust_IKF,na.rm=T)))

#pdf(file='max_abs_error_log10_spaced.pdf',height=3.5,width=5)
ggplot(dat3,aes(N,error,shape=method,fill=method))+geom_point(size = 3)+
  xlab(expression(N))+ylab('maximum abs error')+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=c("#abecec","#80b1d3"))+
  #scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position ="top",#legend.position = c(.75, .5),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.text=element_text(size=15),
        #legend.box.background = element_rect(),legend.box.margin = margin(0,0,0,0),
        axis.text=element_text(size=15),axis.title=element_text(size=20),
        legend.title=element_blank())+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
#dev.off()

