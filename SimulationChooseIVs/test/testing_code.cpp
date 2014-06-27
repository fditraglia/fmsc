//To use this code, you need to paste it at the bottom of the C++
//file contained in the parent directory.

//Testing code - Make some of the member functions available to R
// [[Rcpp::export]]
colvec tsls_est_cpp(mat X, colvec y, mat Z) {
   tsls_fit results(X, y, Z);
   return(results.est());
}

// [[Rcpp::export]]
colvec tsls_SE_textbook_cpp(mat X, colvec y, mat Z) {
   tsls_fit results(X, y, Z);
   return(results.SE_textbook());
}


// [[Rcpp::export]]
colvec tsls_SE_robust_cpp(mat X, colvec y, mat Z) {
   tsls_fit results(X, y, Z);
   return(results.SE_robust());
}

// [[Rcpp::export]]
colvec tsls_SE_center_cpp(mat X, colvec y, mat Z) {
   tsls_fit results(X, y, Z);
   return(results.SE_center());
}


// [[Rcpp::export]]
colvec test_dgp(double g, double r, int n){
  colvec p = ones(3) / 10;
  double b = 1;
  mat Q = eye(3, 3);
  mat V(3,3); 
  V << 1 << 0.5 - g * r << r << endr
    << 0.5  - g * r << 1 << 0 << endr
    << r << 0 << 1 << endr;
  dgp sims(b, p, g, V, Q, n);
  tsls_fit valid(sims.x, sims.y, sims.z1);
  return(valid.est());
}


// [[Rcpp::export]]
List cancor_cpp(mat X, mat Y){
  cancor results(X, Y);
  return List::create(Named("cor") = results.cor,
                      Named("xcoef") = results.xcoef,
                      Named("ycoef") = results.ycoef);
}

// [[Rcpp::export]]
colvec CCIC_test(double g, double r, int n = 500){
  
  //This is the simulation setup from the original 
  //version of the paper (Section 3.4)
  double b = 1;
  colvec p = 0.1 * ones(3);
  mat Q = eye(3, 3);
  mat V(3,3); 
  V << 1 << 0.5 - g * r << r << endr
    << 0.5  - g * r << 1 << 0 << endr
    << r << 0 << 1 << endr;
  
  dgp sims(b, p, g, V, Q, n);
  
  CCIC HallPeixe(sims.x, join_rows(sims.z1, sims.z2));
  colvec out(3);
  out(0) = HallPeixe.BIC();
  out(1) = HallPeixe.AIC();
  out(2) = HallPeixe.HQ();
  return(out);
}

// [[Rcpp::export]]
colvec Andrews_test(double g, double r, int n = 500){
  
  //This is the simulation setup from the original 
  //version of the paper (Section 3.4)
  double b = 1;
  colvec p = 0.1 * ones(3);
  mat Q = eye(3, 3);
  mat V(3,3); 
  V << 1 << 0.5 - g * r << r << endr
    << 0.5  - g * r << 1 << 0 << endr
    << r << 0 << 1 << endr;
  
  dgp sims(b, p, g, V, Q, n);
  
  linearGMM_msc Andrews(sims.x, sims.y,
                        join_rows(sims.z1, sims.z2));
  colvec out(7);
  out(0) = as_scalar(Andrews.est_1step());
  out(1) = as_scalar(Andrews.est_2step());
  out(2) = Andrews.Jstat();
  out(3) = Andrews.pJtest();
  out(4) = Andrews.GMM_AIC();
  out(5) = Andrews.GMM_BIC();
  out(6) = Andrews.GMM_HQ();
  return(out);
}


// [[Rcpp::export]]
List GMMselect_test(double g, double r, int n = 500){
  
  //This is the simulation setup from the original 
  //version of the paper (Section 3.4)
  double b = 1;
  colvec p = 0.1 * ones(3);
  mat Q = eye(3, 3);
  mat V(3,3); 
  V << 1 << 0.5 - g * r << r << endr
    << 0.5  - g * r << 1 << 0 << endr
    << r << 0 << 1 << endr;
  
  //Only two candidate specifications
  umat valid_full(4,2);
  valid_full << 1 << 1 << endr
             << 1 << 1 << endr
             << 1 << 1 << endr
             << 0 << 1 << endr;
  
  dgp sims(b, p, g, V, Q, n);
  
  linearGMM_select results(sims.x, sims.y,
                           join_rows(sims.z1, sims.z2),
                           valid_full);
//  colvec testy = conv_to<colvec>::from(results.moments_BIC());          
  return List::create(
//                      Named("x") = sims.x,
//                      Named("y") = sims.y,
//                      Named("z1") = sims.z1,
//                      Named("z2") = sims.z2,
                      Named("J") = results.J,
                      Named("pJtest") = results.pJtest,
                      Named("AIC") = results.AIC,
                      Named("BIC") = results.BIC,
                      Named("HQ") = results.HQ,
//                      Named("AIC_CCIC") = results.AIC_CCIC,
//                      Named("BIC_CCIC") = results.BIC_CCIC,
//                      Named("HQ_CCIC") = results.HQ_CCIC,
                      Named("onestep") = results.estimates_1step,
//                      Named("twostep") = results.estimates_2step,
                      Named("estAIC") = results.est_AIC(),
                      Named("estBIC") = results.est_BIC(), 
                      Named("estHQ") = results.est_HQ(),
//                      Named("estCCIC_AIC") = results.est_CCIC_AIC(),
//                      Named("estCCIC_HQ") = results.est_CCIC_HQ(),
//                      Named("estCCIC_BIC") = results.est_CCIC_BIC(),
                      Named("momentsAIC") = results.moments_AIC(),
                      Named("momentsBIC") = results.moments_BIC(),
                      Named("momentsHQ") = results.moments_HQ());
//                      Named("momentsCCIC_AIC") = results.moments_CCIC_AIC(),
//                      Named("momentsCCIC_BIC") = results.moments_CCIC_BIC(),
//                      Named("momentsCCIC_HQ") = results.moments_CCIC_HQ());                       
}

// [[Rcpp::export]]
List fmsc_test(mat x, colvec y, mat z1, mat z2, 
                      umat candidates){
  
  fmsc_chooseIV test(x, y, z1, z2, candidates);
  
  //The target parameter is the OLS slope coefficient
  colvec w(2);
  w << 0 << endr
    << 1 << endr;
  
  return List::create(Named("tau") = test.tau,
                      Named("Psi") = test.Psi,
                      Named("tau.outer") = test.tau_outer_est,
                      Named("Bias.mat") = test.Bias_mat,
                      Named("Estimates") = test.estimates,
                      Named("K") = test.K,
                      Named("Omega") = test.Omega,
                      Named("sq.bias.inner") = test.sqbias_inner,
                      Named("avar.inner") = test.avar_inner,
                      Named("mu") = test.mu(w),
                      Named("mu.valid") = test.mu_valid(w),
                      Named("mu.full") = test.mu_full(w),
                      Named("abias.sq") = test.abias_sq(w),
                      Named("avar") = test.avar(w),
                      Named("fmsc") = test.fmsc(w),
                      Named("fmsc.pos") = test.fmsc_pos(w),
                      Named("mu.fmsc") = test.mu_fmsc(w),
                      Named("mu.fmsc.pos") = test.mu_fmsc_pos(w));
}