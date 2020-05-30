function [etak]=eta_selection(alpha_uni,epsi,a,b,mode_error_n,cov_w_corr_c_n)
for d=1:n_d_c
for t=1:n_t
for j=1:n_g    
mode_error_np=mode_error_n{d,t};
cov_w_corr_c_np=cov_w_corr_c_n{d,t};    
xl=(1/(1-epsi))^(1/alpha_uni);
mu_zero=((alpha_uni+1)/alpha_uni)*(-mode_error_np)'*a;
sigma_zero=((alpha_uni+2)/alpha_uni)*a'*(cov_w_corr_c_np+mode_error_np*mode_error_np')*a;
xu=(alpha_uni*(1-epsi)^((alpha_uni+1)/alpha_uni)*(sigma_zero-(mu_zero)^2))/(2*epsi*(b-a'*mode_error_np)^2);
RR=0.5*(sqrt(5)-1);
dd=RR*(xu-xl);
x1=xu-dd;
x2=xl+dd;
eta=x1;
violation=violate(alpha_uni,epsi,eta,a,b,mode_error_np,cov_w_corr_c_np);
f1=violation;
eta=x2;
violation=violate(alpha_uni,epsi,eta,a,b,mode_error_np,cov_w_corr_c_np);
f2=violation;
tol=1e-4;
err=inf;
while err>tol
    if f1>f2
        xu=x2;
        fu=f2;
        x2=x1;
        f2=f1;
        dd=RR*(xu-xl);
        x1=xu-dd;
        violation=violate(alpha_uni,epsi,eta,a,b,mode_error_np,cov_w_corr_c_np);
        f1=violation;
    elseif f1<f2
        xl=x1;
        fl=f1;
        x1=x2;
        f1=f2;
        dd=RR*(xu-xl);
        x2=xl+dd;
        violation=violate(alpha_uni,epsi,eta,a,b,mode_error_np,cov_w_corr_c_np);
        f2=violation;
    else
        xl=(x1+x2)/2;
        xu=xl;
    end
    err=2*abs(xu-xl)/(xu+xl);
end

violation=violate1up(alpha_uni,epsi,eta,a,b,mode_error_np,cov_w_corr_c_np);
fu1up(j,t,d)=violation;
% if fu1up(j,t,d)<=0 || xu>10^5
if fu1up(j,t,d)<=0
    etak(j,t,d)=NaN; 

else
    etak(j,t,d)=xu;

end
end
end
end
end