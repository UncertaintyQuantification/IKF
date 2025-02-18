
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*- 

// we only include RcppEigen.h which pulls Rcpp.h in for us 
#include <RcppEigen.h> 
#include <Rcpp.h> 
#include <cmath> 
#include "ctools.h"
//#include <valarray>

// [[Rcpp::depends(RcppEigen)]] 

using namespace Rcpp;
using namespace std;
using namespace Eigen; 


//July 16, 2018
////Construct_W0_matern_5_2 
// [[Rcpp::export]] 
MatrixXd Construct_W0_matern_5_2(const double sigma2, const double lambda){ 
  //int num_dim=sigma2.size(); 
  
  Eigen::MatrixXd W0= Eigen::MatrixXd::Zero(3,3); 
  //Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
  
  W0(0,0)=sigma2; 
  W0(0,2)=W0(2,0)=-sigma2*pow(lambda,2.0)/3.0; 
  W0(1,1)=sigma2*pow(lambda,2.0)/3.0; 
  W0(2,2)=sigma2*pow(lambda,4.0); 
  
  return W0; 
} 

//July 16, 2018
////Construct_W0_exp 
// [[Rcpp::export]] 
MatrixXd Construct_W0_exp(const double sigma2, const double lambda){ 
  //int num_dim=sigma2.size(); 
  
  Eigen::MatrixXd W0= Eigen::MatrixXd::Zero(1,1); 
  
  W0(0,0)=sigma2; 
  
  return W0; 
} 



////Construct_G_matern_5_2 
// [[Rcpp::export]] 
List Construct_G_matern_5_2(Eigen::VectorXd delta_x, double lambda){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  //int num_dim=lambda.size();  
  List GG(num_obs);  
  GG[0]=Eigen::MatrixXd::Zero(3,3); 
  
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
  
  // num_dim list, each is 3(num_obs)\times 3 list 
  // for(int i_GG=0;i_GG<num_dim;i_GG++){ 
  //   Eigen::MatrixXd d= Eigen::MatrixXd::Zero(num_obs,9);  //the first row has all zeros  
  for(int j_GG=0;j_GG<(num_obs-1);j_GG++){ 
    int j_GG_1=j_GG+1;    
    d(0,0)=pow(lambda,2.0)*pow(delta_x[j_GG],2.0)+2*lambda*delta_x[j_GG]+2; 
    d(1,0)=-pow(lambda,3.0)*pow(delta_x[j_GG],2.0); 
    d(2,0)=pow(lambda,4.0)*pow(delta_x[j_GG],2.0)-2*pow(lambda,3.0)*delta_x[j_GG]; 
    d(0,1)=2*(lambda*pow(delta_x[j_GG],2.0)+delta_x[j_GG]); 
    d(1,1)=-2*(pow(lambda,2.0)*pow(delta_x[j_GG],2.0)-lambda*delta_x[j_GG]-1); 
    d(2,1)=2*(pow(lambda,3.0)*pow(delta_x[j_GG],2.0)-3*pow(lambda,2.0)*delta_x[j_GG]); 
    d(0,2)=pow(delta_x[j_GG],2); 
    d(1,2)=2*delta_x[j_GG]-lambda*pow(delta_x[j_GG],2.0); 
    d(2,2)=pow(lambda,2.0)*pow(delta_x[j_GG],2.0)-4*lambda*delta_x[j_GG]+2;     
    d=exp(-lambda*delta_x[j_GG])/2.0*d;
    GG[j_GG_1]=d; 
  } 
  // GG[i_GG]=d; 
  //} 
  return GG; 
} 

////Construct_G_exp
// [[Rcpp::export]] 
List Construct_G_exp(Eigen::VectorXd delta_x, double lambda){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  //int num_dim=lambda.size();  
  List GG(num_obs);  
  GG[0]=Eigen::MatrixXd::Zero(1,1); 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(1,1); 
  
  for(int j_GG=0;j_GG<(num_obs-1);j_GG++){ 
    d(0,0)=exp(-delta_x[j_GG]*lambda);
    GG[j_GG+1]=d; 
  }
  
  return GG;
}

////Construct_W_matern_5_2  
// [[Rcpp::export]] 
List Construct_W_matern_5_2(double sigma2,Eigen::VectorXd delta_x, double lambda, MatrixXd W0){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  //int num_dim=sigma2.size();  
  List Wi(num_obs);  
  Wi[0]=W0; 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
  
  
  // List Wi(num_obs); 
  //for(int i_Wi=0;i_Wi<num_dim;i_Wi++){ 
  //Eigen::MatrixXd d= Eigen::MatrixXd::Zero(num_obs,9);   
  double  lambda_delta_x;
  double exp_neg_2_lambda_delta_x;
  int  j_Wi_1;
  for(int j_Wi=0;j_Wi<(num_obs-1);j_Wi++){ 
    j_Wi_1= j_Wi+1; 
    lambda_delta_x=lambda*delta_x[j_Wi];  //close and jump then it is... 
    exp_neg_2_lambda_delta_x=exp(-2*lambda_delta_x); 
    
    d(0,0)=(exp_neg_2_lambda_delta_x*(3+6*lambda_delta_x+6*pow(lambda_delta_x,2.0)+4*pow(lambda_delta_x,3.0)+2*pow(lambda_delta_x,4.0))-3 )/(-4*pow(lambda,5.0)); 
    d(1, 0)=  d(0, 1)=exp_neg_2_lambda_delta_x*pow(delta_x[j_Wi],4.0)/2.0; 
    d(2, 0)=  d(0, 2)=(exp_neg_2_lambda_delta_x*(1+2*lambda_delta_x+2*pow(lambda_delta_x,2.0)+4*pow(lambda_delta_x,3.0)-2*pow(lambda_delta_x,4.0))-1 )/(4*pow(lambda,3.0)); 
    d(1, 1)= (exp_neg_2_lambda_delta_x*(1+2*lambda_delta_x+2*pow(lambda_delta_x,2.0)-4*pow(lambda_delta_x,3.0)+2*pow(lambda_delta_x,4.0))-1 )/(-4*pow(lambda,3.0)); 
    d(2, 1)=  d(1, 2)=exp_neg_2_lambda_delta_x*pow(delta_x[j_Wi],2.0)*(4-4*lambda_delta_x+pow(lambda_delta_x,2.0) )/2.0; 
    d(2, 2)=(exp_neg_2_lambda_delta_x*(-3+10*lambda_delta_x-22*pow(lambda_delta_x,2.0)+12*pow(lambda_delta_x,3.0)-2*pow(lambda_delta_x,4.0))+3 )/(4*lambda)     ;  
    d=d*(4*sigma2*pow(lambda,5.0)/3.0); 
    Wi[j_Wi_1]=d; 
    
    //} 
  } 
  return Wi; 
}


////Construct_W_matern_5_2  
// [[Rcpp::export]] 
List Construct_W_exp(double sigma2, Eigen::VectorXd delta_x, double lambda, MatrixXd W0){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  //int num_dim=sigma2.size();  
  List Wi(num_obs);  
  Wi[0]=W0; 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(1,1); 
  
  for(int j_Wi=0;j_Wi<(num_obs-1);j_Wi++){ 
    d(0,0)=1-exp(-2*delta_x[j_Wi]*lambda);
    Wi[j_Wi+1]=d;
  }
  
  return Wi;
}

////Get_Q_K  
// [[Rcpp::export]]
List Get_Q_K(const List GG,const List  W,const Eigen::MatrixXd C0,const double VV){

  int n=GG.size();
  int k=C0.rows();

  Eigen::VectorXd Q=Eigen::VectorXd::Zero(n);
  Eigen::MatrixXd K=Eigen::MatrixXd::Zero(n,k);
  Eigen::MatrixXd C=C0;

  Eigen::MatrixXd GG_matrix;
  Eigen::MatrixXd W_matrix;

  Eigen::MatrixXd RR;


  // num_dim list, each is 3(num_obs)\times 3 list
  for(int t=0;t<n;t++){
    GG_matrix=GG[t];
    W_matrix=W[t];
    RR=GG_matrix*C*GG_matrix.transpose()+W_matrix;
    //Q[t]=RR(0,0);
    Q[t]=RR(0,0)+VV;
    K.row(t)=RR.col(0).transpose()/Q[t];
    C=RR-RR.col(0)*RR.row(0)/Q[t];
  }

  List return_list;
  return_list.push_back(Q);
  return_list.push_back(K);


  return return_list;
}


// [[Rcpp::export]]
VectorXd Get_L_t_y(const List GG,const Eigen::VectorXd Q, const Eigen::MatrixXd K,  const Eigen::VectorXd output){
  int n=GG.size();

  Eigen::VectorXd sqrt_Q=Q.array().sqrt();


  Eigen::VectorXd res=Eigen::VectorXd::Zero(n);

  res(n-1)=sqrt_Q(n-1)*output(n-1);

  Eigen::MatrixXd GG_matrix;
  if(n>=2){
    GG_matrix=GG[n-1];
    Eigen::VectorXd u=GG_matrix*K.row(n-2).transpose();
    res(n-2)=sqrt_Q(n-2)*(u(0)*output(n-1)+output(n-2)); // this is for Matern 2.5
    if(n>=3){
      //Nov 2022
      Eigen::MatrixXd GG_matrix_plus_1;
      Eigen::MatrixXd g;
      Eigen::VectorXd GG_K;
      Eigen::VectorXd g_GG_K;
      g=GG_matrix.row(0)*output(n-1);
      for(int t=n-3;t>=0;t--){
        GG_matrix_plus_1=GG[t+1];
        GG_K=GG_matrix_plus_1*K.row(t).transpose();
        g_GG_K=g*GG_K;
        res(t)=sqrt_Q(t)*(g_GG_K(0)+GG_K(0)*output(t+1)+output(t));
        g=g*GG_matrix_plus_1+GG_matrix_plus_1.row(0)*output(t+1);
      }
    }
  }

  return res;
}

// [[Rcpp::export]]
VectorXd Get_L_y(const List GG,const Eigen::VectorXd Q, const Eigen::MatrixXd K,  const Eigen::VectorXd output){
  int n=GG.size();

  Eigen::VectorXd sqrt_Q=Q.array().sqrt();

  Eigen::MatrixXd m=Eigen::VectorXd::Zero(K.cols());

  Eigen::VectorXd a;

  //Eigen::VectorXd Y_minus_a_1_scaled_vec=Eigen::VectorXd::Zero(n);
  //Eigen::VectorXd Y_minus_a_1=Eigen::VectorXd::Zero(n);

  Eigen::VectorXd res=Eigen::VectorXd::Zero(n);
  Eigen::MatrixXd GG_matrix;

  for(int t=0;t<n;t++){
    GG_matrix=GG[t];
    a=GG_matrix*m;
    res[t]=a[0]+output[t]*sqrt_Q[t];

    //Y_minus_a_1[t]=(res[t]-a[0]);

    // Y_minus_a_1_scaled_vec[t]=(res[t]-a[0])/sqrt_Q[t];
    m=a+K.row(t).transpose()*(res[t]-a[0]);
  }
  //return output;

  return res;
}
// 
// [[Rcpp::export]]
VectorXd Get_R_y(const List GG,const Eigen::VectorXd Q, const Eigen::MatrixXd K,  const Eigen::VectorXd output){

  Eigen::VectorXd tilde_z = Get_L_t_y(GG,Q,K,output);
  Eigen::VectorXd res = Get_L_y(GG,Q,K,tilde_z);

  return res;
}



// [[Rcpp::export]]
VectorXd A_t_times_x_sparse_missing(const VectorXd output, const Eigen::VectorXd miss_ind, const int N){

  //miss_ind starts at 1 (use R syntax)

  VectorXd A_t_x = VectorXd::Zero(N);

  int j_x=0; //point to output
  int j_m=0; //point to miss_ind

  for(int i=0; i<N; i++){
    if(j_m >= miss_ind.rows()){
      A_t_x[i]=output[j_x];
      j_x += 1;
    }else{
      if(i==miss_ind[j_m]-1){
        A_t_x[i]=0;
        j_m += 1;
      }else{
        A_t_x[i]=output[j_x];
        j_x += 1;
      }
    }

  }
  return A_t_x;

}


// [[Rcpp::export]]
VectorXd R_times_z_kronecker(Eigen::VectorXd z, const Eigen::VectorXd param, const int n1, const int n2,
                             const Eigen::VectorXd delta_x1,const Eigen::VectorXd delta_x2,
                             const double tilde_nu, const  String kernel_type){

  //preconstructure
  double VV=tilde_nu; // tilde_nu to stablize the computation

  //for R1
  double gamma1=1.0/exp(param[0]);
  Eigen::MatrixXd W01;
  List GG1;
  List W1;
  List Q_K1;
  double lambda1=0;

  //for R2
  double gamma2=1.0/exp(param[1]);
  Eigen::MatrixXd W02;
  List GG2;
  List W2;
  List Q_K2;
  double lambda2=0;

  if(kernel_type=="matern_5_2"){
    //for R1
    lambda1=sqrt(5.0)/gamma1;
    W01=Construct_W0_matern_5_2(1.0,lambda1);
    GG1=Construct_G_matern_5_2(delta_x1,lambda1);
    W1=Construct_W_matern_5_2(1.0,delta_x1,lambda1,W01);

    //for R2
    lambda2=sqrt(5.0)/gamma2;
    W02=Construct_W0_matern_5_2(1.0,lambda2);
    GG2=Construct_G_matern_5_2(delta_x2,lambda2);
    W2=Construct_W_matern_5_2(1.0,delta_x2,lambda2,W02);

  }else if(kernel_type=="exp"){
    //for R1
    lambda1=1.0/gamma1;
    W01=Construct_W0_exp(1.0,lambda1);
    GG1=Construct_G_exp(delta_x1,lambda1);
    W1=Construct_W_exp(1.0,delta_x1,lambda1,W01);

    //for R2
    lambda2=1.0/gamma2;
    W02=Construct_W0_exp(1.0,lambda2);
    GG2=Construct_G_exp(delta_x2,lambda2);
    W2=Construct_W_exp(1.0,delta_x2,lambda2,W02);
  }

  Q_K1=Get_Q_K(GG1,W1,W01,VV);
  Eigen::VectorXd Q1=Q_K1[0];
  Eigen::MatrixXd K1=Q_K1[1];

  Q_K2=Get_Q_K(GG2,W2,W02,VV);
  Eigen::VectorXd Q2=Q_K2[0];
  Eigen::MatrixXd K2=Q_K2[1];

  //Sep 2022
  Q1 = (Q1.array() < 0).select(0, Q1); //truncate them to be zero to avoid NA in singular case, but error could be large if it is too singular
  Q2 = (Q2.array() < 0).select(0, Q2);

  //Eigen::VectorXd sqrt_Q1=Q1.array().sqrt();
  //Eigen::VectorXd sqrt_Q2=Q2.array().sqrt();



  int N = n1*n2;
  Eigen::Map<Eigen::MatrixXd> mat_n2_n1(z.data(), n2, n1);

  //compute R2 %*% mat(z) %*% R1
  MatrixXd R2_Z=MatrixXd::Zero(n2,n1);
  MatrixXd R2_Z_R1=MatrixXd::Zero(n2,n1);
  Eigen::VectorXd Rz_res;
  Eigen::VectorXd temp;


  // First loop
  for (int c = 0; c < n1; c++) {
    //Get_R_y(GG2,Q2,K2,mat_n2_n1.col(c))
    Rz_res = Get_R_y(GG2,Q2,K2,mat_n2_n1.col(c));
    temp = tilde_nu * mat_n2_n1.col(c);
    R2_Z.col(c) = Rz_res - temp;
  }

  // Second loop
  for (int r = 0; r < n2; r++) {
    //R2_Z_R1.row(r)=R_times_z_pre_constructed(R2_Z.row(r),sqrt_Q2,K2,GG2)-tilde_nu*R2_Z.row(r); can't directly do this (potentially due to row and column operations)
    Rz_res = Get_R_y(GG1,Q1,K1,R2_Z.row(r));
    temp = tilde_nu * R2_Z.row(r);
    R2_Z_R1.row(r) = Rz_res - temp;
  }
  //convert R2_Z_R1 back to vector
  Map<VectorXd> y_KF(R2_Z_R1.data(), N);

  return y_KF;
}
// 
// 
// 
// 
// //Feb 2024
// [[Rcpp::export]]
List fast_pred_sparse_CG_missing(VectorXd param, const  String kernel_type, const VectorXd delta_x1, const VectorXd delta_x2,
                                 const VectorXd output, int n1, int n2, const VectorXd miss_ind,
                                 const double tilde_nu, float tol=1e-6, int maxIte = 1000){

  //param=log(c(beta1,beta2,tau))
  int N0=output.size();
  int N=n1*n2;
  double  tau=exp(param[2]);

  VectorXd vec_N = VectorXd::Zero(N);
  //MatrixXd mat_n2_n1;
  Eigen::MatrixXd R2_Z;
  Eigen::MatrixXd R2_Z_R1;
  VectorXd Rz_res;
  VectorXd temp;
  //VectorXd y_KF;
  VectorXd vec_N0 = VectorXd::Zero(N0);

  //variables for CG itself
  VectorXd x=VectorXd::Zero(N0);
  VectorXd r = output;
  VectorXd p = r;
  double rs_old = (r.transpose() * r).value();
  double rs_new=INFINITY;//1000;//1.0;
  double rs_ratio;
  double alpha;

  int ite = 0;
  VectorXd resid = VectorXd::Zero(maxIte);


  //preconstructure
  double VV=tilde_nu; // tilde_nu to stablize the computation

  //for R1
  double gamma1=1.0/exp(param[0]);
  Eigen::MatrixXd W01;
  List GG1;
  List W1;
  List Q_K1;
  double lambda1=0;

  //for R2
  double gamma2=1.0/exp(param[1]);
  Eigen::MatrixXd W02;
  List GG2;
  List W2;
  List Q_K2;
  double lambda2=0;

  if(kernel_type=="matern_5_2"){
    //for R1
    lambda1=sqrt(5.0)/gamma1;
    W01=Construct_W0_matern_5_2(1.0,lambda1);
    GG1=Construct_G_matern_5_2(delta_x1,lambda1);
    W1=Construct_W_matern_5_2(1.0,delta_x1,lambda1,W01);

    //for R2
    lambda2=sqrt(5.0)/gamma2;
    W02=Construct_W0_matern_5_2(1.0,lambda2);
    GG2=Construct_G_matern_5_2(delta_x2,lambda2);
    W2=Construct_W_matern_5_2(1.0,delta_x2,lambda2,W02);

  }else if(kernel_type=="exp"){
    //for R1
    lambda1=1.0/gamma1;
    W01=Construct_W0_exp(1.0,lambda1);
    GG1=Construct_G_exp(delta_x1,lambda1);
    W1=Construct_W_exp(1.0,delta_x1,lambda1,W01);

    //for R2
    lambda2=1.0/gamma2;
    W02=Construct_W0_exp(1.0,lambda2);
    GG2=Construct_G_exp(delta_x2,lambda2);
    W2=Construct_W_exp(1.0,delta_x2,lambda2,W02);
  }

  Q_K1=Get_Q_K(GG1,W1,W01,VV);
  Eigen::VectorXd Q1=Q_K1[0];
  Eigen::MatrixXd K1=Q_K1[1];

  Q_K2=Get_Q_K(GG2,W2,W02,VV);
  Eigen::VectorXd Q2=Q_K2[0];
  Eigen::MatrixXd K2=Q_K2[1];

  //Sep 2022
  Q1 = (Q1.array() < 0).select(0, Q1); //truncate them to be zero to avoid NA in singular case, but error could be large if it is too singular
  Q2 = (Q2.array() < 0).select(0, Q2);

  //Eigen::VectorXd sqrt_Q1=Q1.array().sqrt();
  //Eigen::VectorXd sqrt_Q2=Q2.array().sqrt();

  //end
  while((ite < maxIte) && (rs_new > tol)){

    //first step: t(A) %*% v
    int j_x=0; //point to v
    int j_m=0; //point to miss_ind

    for(int i=0; i<N; i++){
      if(j_m >= miss_ind.rows()){
        vec_N[i]=p[j_x];
        j_x += 1;
      }else{
        if(i==miss_ind[j_m]-1){
          vec_N[i]=0;
          j_m += 1;
        }else{
          vec_N[i]=p[j_x];
          j_x += 1;
        }
      }

    }

    //second step: kronecker(R1, R2) %*% v = vec(R2 %*% mat(v) %*% R1)
    //convert vec_N to matrix
    Map<MatrixXd> mat_n2_n1(vec_N.data(), n2, n1);

    //compute R2 %*% mat(vec_N) %*% R1
    R2_Z=MatrixXd::Zero(n2,n1);
    R2_Z_R1=MatrixXd::Zero(n2,n1);

    // First loop
    for (int c_i = 0; c_i < n1; c_i++) {
      Eigen::VectorXd Rz_res = Get_R_y(GG2,Q2,K2,mat_n2_n1.col(c_i));
      Eigen::VectorXd temp = tilde_nu * mat_n2_n1.col(c_i);
      R2_Z.col(c_i) = Rz_res - temp;
    }

    // Second loop
    for (int r_i = 0; r_i < n2; r_i++) {
      //R2_Z_R1.row(r_i)=R_times_z_pre_constructed(R2_Z.row(r_i),sqrt_Q2,K2,GG2)-tilde_nu*R2_Z.row(r_i); can't directly do this (not sure why)
      Eigen::VectorXd Rz_res = Get_R_y(GG1,Q1,K1,R2_Z.row(r_i));
      Eigen::VectorXd temp = tilde_nu * R2_Z.row(r_i);
      R2_Z_R1.row(r_i) = Rz_res - temp;
    }
    //convert R2_Z_R1 back to vector
    Map<VectorXd> y_KF(R2_Z_R1.data(), N);

    //third step: A %*% v
    j_x=0; //point to A %*% v
    j_m=0; //point to miss_ind

    for(int i=0; i<N; i++){
      if(j_m >= miss_ind.rows()){
        vec_N0[j_x]=y_KF[i];
        j_x += 1;
      }else{
        if(i==miss_ind[j_m]-1){
          j_m += 1;
        }else{
          vec_N0[j_x]=y_KF[i];
          j_x += 1;
        }
      }
    }
    vec_N0=vec_N0*tau+p;

    alpha = rs_old / (p.transpose() * vec_N0).value();
    x += alpha*p;
    r -= alpha*vec_N0;
    rs_new = (r.transpose() * r).value();
    rs_ratio = rs_new / rs_old;
    p = r + rs_ratio * p;
    rs_old = rs_new;
    resid[ite] = rs_new; //mean((res2[[1]]-z)^2)
    ite++;

  }


  List ans_list;

  ans_list.push_back(x);
  ans_list.push_back(resid.head(ite));
  ans_list.push_back(ite);


  return ans_list;
}
// 
// //Feb 2024
// [[Rcpp::export]]
List fast_pred_sparse_CG_missing_with_ini(VectorXd param, const  String kernel_type, const VectorXd delta_x1, const VectorXd delta_x2,
                                          const VectorXd output, int n1, int n2, const VectorXd miss_ind, VectorXd x,
                                          const double tilde_nu, float tol=1e-6, int maxIte = 1000){

  //param=log(c(beta1,beta2,tau))
  int N0=output.size();
  int N=n1*n2;
  double  tau=exp(param[2]);

  VectorXd vec_N = VectorXd::Zero(N);
  //MatrixXd mat_n2_n1;
  Eigen::MatrixXd R2_Z;
  Eigen::MatrixXd R2_Z_R1;
  VectorXd Rz_res;
  VectorXd temp;
  //VectorXd y_KF;
  VectorXd vec_N0 = VectorXd::Zero(N0);


  //preconstructure
  double VV=tilde_nu; // tilde_nu to stablize the computation

  //for R1
  double gamma1=1.0/exp(param[0]);
  Eigen::MatrixXd W01;
  List GG1;
  List W1;
  List Q_K1;
  double lambda1=0;

  //for R2
  double gamma2=1.0/exp(param[1]);
  Eigen::MatrixXd W02;
  List GG2;
  List W2;
  List Q_K2;
  double lambda2=0;

  if(kernel_type=="matern_5_2"){
    //for R1
    lambda1=sqrt(5.0)/gamma1;
    W01=Construct_W0_matern_5_2(1.0,lambda1);
    GG1=Construct_G_matern_5_2(delta_x1,lambda1);
    W1=Construct_W_matern_5_2(1.0,delta_x1,lambda1,W01);

    //for R2
    lambda2=sqrt(5.0)/gamma2;
    W02=Construct_W0_matern_5_2(1.0,lambda2);
    GG2=Construct_G_matern_5_2(delta_x2,lambda2);
    W2=Construct_W_matern_5_2(1.0,delta_x2,lambda2,W02);

  }else if(kernel_type=="exp"){
    //for R1
    lambda1=1.0/gamma1;
    W01=Construct_W0_exp(1.0,lambda1);
    GG1=Construct_G_exp(delta_x1,lambda1);
    W1=Construct_W_exp(1.0,delta_x1,lambda1,W01);

    //for R2
    lambda2=1.0/gamma2;
    W02=Construct_W0_exp(1.0,lambda2);
    GG2=Construct_G_exp(delta_x2,lambda2);
    W2=Construct_W_exp(1.0,delta_x2,lambda2,W02);
  }

  Q_K1=Get_Q_K(GG1,W1,W01,VV);
  Eigen::VectorXd Q1=Q_K1[0];
  Eigen::MatrixXd K1=Q_K1[1];

  Q_K2=Get_Q_K(GG2,W2,W02,VV);
  Eigen::VectorXd Q2=Q_K2[0];
  Eigen::MatrixXd K2=Q_K2[1];

  //Sep 2022
  Q1 = (Q1.array() < 0).select(0, Q1); //truncate them to be zero to avoid NA in singular case, but error could be large if it is too singular
  Q2 = (Q2.array() < 0).select(0, Q2);

  //Eigen::VectorXd sqrt_Q1=Q1.array().sqrt();
  //Eigen::VectorXd sqrt_Q2=Q2.array().sqrt();

  //variables for CG itself
  //VectorXd x=VectorXd::Zero(N0);

  //Compute output - R_y * x
  //first step: t(A) %*% v
  int j_x=0; //point to v
  int j_m=0; //point to miss_ind

  for(int i=0; i<N; i++){
    if(j_m >= miss_ind.rows()){
      vec_N[i]=x[j_x];
      j_x += 1;
    }else{
      if(i==miss_ind[j_m]-1){
        vec_N[i]=0;
        j_m += 1;
      }else{
        vec_N[i]=x[j_x];
        j_x += 1;
      }
    }
  }

  //second step: kronecker(R1, R2) %*% v = vec(R2 %*% mat(v) %*% R1)
  //convert vec_N to matrix
  Map<MatrixXd> mat_n2_n1(vec_N.data(), n2, n1);

  //compute R2 %*% mat(vec_N) %*% R1
  R2_Z=MatrixXd::Zero(n2,n1);
  R2_Z_R1=MatrixXd::Zero(n2,n1);

  // First loop
  for (int c = 0; c < n1; c++) {
    Eigen::VectorXd Rz_res = Get_R_y(GG2,Q2,K2,mat_n2_n1.col(c));
    Eigen::VectorXd temp = tilde_nu * mat_n2_n1.col(c);
    R2_Z.col(c) = Rz_res - temp;
  }

  // Second loop
  for (int r = 0; r < n2; r++) {
    //R2_Z_R1.row(r)=R_times_z_pre_constructed(R2_Z.row(r),sqrt_Q2,K2,GG2)-tilde_nu*R2_Z.row(r); can't directly do this (not sure why)
    Eigen::VectorXd Rz_res = Get_R_y(GG1,Q1,K1,R2_Z.row(r));
    Eigen::VectorXd temp = tilde_nu * R2_Z.row(r);
    R2_Z_R1.row(r) = Rz_res - temp;
  }
  //convert R2_Z_R1 back to vector
  Map<VectorXd> y_KF(R2_Z_R1.data(), N);

  //third step: A %*% v
  j_x=0; //point to A %*% v
  j_m=0; //point to miss_ind

  for(int i=0; i<N; i++){
    if(j_m >= miss_ind.rows()){
      vec_N0[j_x]=y_KF[i];
      j_x += 1;
    }else{
      if(i==miss_ind[j_m]-1){
        j_m += 1;
      }else{
        vec_N0[j_x]=y_KF[i];
        j_x += 1;
      }
    }
  }
  vec_N0=vec_N0*tau+x;

  VectorXd r = output-vec_N0; //output -R_y*x


  VectorXd p = r;
  double rs_old = (r.transpose() * r).value();
  double rs_new=1000;//1.0;
  double rs_ratio;
  double alpha;

  int ite = 0;
  VectorXd resid = VectorXd::Zero(maxIte);

  //end
  while((ite < maxIte) && (rs_new > tol)){

    //first step: t(A) %*% v
    int j_x=0; //point to v
    int j_m=0; //point to miss_ind

    for(int i=0; i<N; i++){
      if(j_m >= miss_ind.rows()){
        vec_N[i]=p[j_x];
        j_x += 1;
      }else{
        if(i==miss_ind[j_m]-1){
          vec_N[i]=0;
          j_m += 1;
        }else{
          vec_N[i]=p[j_x];
          j_x += 1;
        }
      }

    }

    //second step: kronecker(R1, R2) %*% v = vec(R2 %*% mat(v) %*% R1)
    //convert vec_N to matrix
    Map<MatrixXd> mat_n2_n1(vec_N.data(), n2, n1);

    //compute R2 %*% mat(vec_N) %*% R1
    R2_Z=MatrixXd::Zero(n2,n1);
    R2_Z_R1=MatrixXd::Zero(n2,n1);

    // First loop
    for (int c_i = 0; c_i < n1; c_i++) {
      Eigen::VectorXd Rz_res = Get_R_y(GG2,Q2,K2,mat_n2_n1.col(c_i));
      Eigen::VectorXd temp = tilde_nu * mat_n2_n1.col(c_i);
      R2_Z.col(c_i) = Rz_res - temp;
    }

    // Second loop
    for (int r_i = 0; r_i < n2; r_i++) {
      //R2_Z_R1.row(r_i)=R_times_z_pre_constructed(R2_Z.row(r_i),sqrt_Q2,K2,GG2)-tilde_nu*R2_Z.row(r_i); can't directly do this (not sure why)
      Eigen::VectorXd Rz_res = Get_R_y(GG1,Q1,K1,R2_Z.row(r_i));
      Eigen::VectorXd temp = tilde_nu * R2_Z.row(r_i);
      R2_Z_R1.row(r_i) = Rz_res - temp;
    }
    //convert R2_Z_R1 back to vector
    Map<VectorXd> y_KF(R2_Z_R1.data(), N);

    //third step: A %*% v
    j_x=0; //point to A %*% v
    j_m=0; //point to miss_ind

    for(int i=0; i<N; i++){
      if(j_m >= miss_ind.rows()){
        vec_N0[j_x]=y_KF[i];
        j_x += 1;
      }else{
        if(i==miss_ind[j_m]-1){
          j_m += 1;
        }else{
          vec_N0[j_x]=y_KF[i];
          j_x += 1;
        }
      }

    }
    vec_N0=vec_N0*tau+p;

    alpha = rs_old / (p.transpose() * vec_N0).value();
    x += alpha*p;
    r -= alpha*vec_N0;
    rs_new = (r.transpose() * r).value();
    rs_ratio = rs_new / rs_old;
    p = r + rs_ratio * p;
    rs_old = rs_new;
    resid[ite] = rs_new; //mean((res2[[1]]-z)^2)
    ite++;

  }


  List ans_list;

  ans_list.push_back(x);
  ans_list.push_back(resid.head(ite));
  ans_list.push_back(ite);


  return ans_list;
}

