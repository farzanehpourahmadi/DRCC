function violation=violation(alpha_uni,epsi,eta,a,b,mode_error_np,cov_w_corr_c_np)

mu_zero=((alpha_uni+1)/alpha_uni)*(-mode_error_np)'*a;
sigma_zero=((alpha_uni+2)/alpha_uni)*a'*(cov_w_corr_c_np+mode_error_np*mode_error_np')*a;
violation=-((eta*(b-a'*mode_error_np)-mu_zero)^2-...
((1-epsi-eta^(-alpha_uni))/epsi)*(sigma_zero-(mu_zero)^2));       



end