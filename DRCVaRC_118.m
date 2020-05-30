%load('parameters_118bus.mat')
tic
n_w=2; %number of m plants
n_ge=19; %number of generators
n_gc=22; %number of candidate generators
n_g=n_gc+n_ge;
n_l=99; %number of loads
n_line=186;
n_d_c=3;
n_t=24;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      00;
mu(n_t,n_w)=0;

alpha_uni=2;

% calculating mean and mode vecotrs and covariance matrix 
 for d=1:n_d_c
 for t=1:n_t  
     %scenario 1
     m_c_n{d,t}=3*m_c{d,t};
     mode_error_n{d,t}=3*mode_error{d,t};
     mode_error_c_n{d,t}=3*mode_error_c{d,t};
     cov_w_corr_c_n{d,t}=3*3*cov_w_corr_c{d,t};
     cov_corr_c_n{d,t}=3*3*cov_corr_c{d,t};
     sigma_uni_n{d,t}=(((alpha_uni+2)/alpha_uni)*cov_w_corr_c_n{d,t})-((1/(alpha_uni)^2)*mode_error_n{d,t}*mode_error_n{d,t}');
     sigma_uni_c_n{d,t}=(((alpha_uni+2)/alpha_uni)*cov_corr_c_n{d,t})-((1/(alpha_uni)^2)*mode_error_n_c{d,t}*mode_error_n_c{d,t}');
     % scenario 2   
     m_c_n2{d,t}=3*0.9*m_c{d,t};
     mode_error_n2{d,t}=3*0.9*mode_error{d,t};
     mode_error_c_n2{d,t}=3*0.9*mode_error_c{d,t};
     cov_w_corr_c_n2{d,t}=3*3*0.9*0.9*cov_w_corr_c{d,t};
     cov_corr_c_n2{d,t}=3*3*0.9*0.9*cov_corr_c{d,t};
     sigma_uni_n2{d,t}=(((alpha_uni+2)/alpha_uni)*cov_w_corr_c_n2{d,t})-((1/(alpha_uni)^2)*mode_error_n2{d,t}*mode_error_n2{d,t}');
     sigma_uni_c_n2{d,t}=(((alpha_uni+2)/alpha_uni)*cov_corr_c_n2{d,t})-((1/(alpha_uni)^2)*mode_error_c_n2{d,t}*mode_error_c_n2{d,t}');
 end
 end


epsi=0.1;
Hg=PTDF*busgen; % sensitivity of each line to each generation unit
Hw=PTDF*buswind; % sensitivity of each line to each m plant
Hl=PTDF*busload;

% definig variables
v=binvar(n_gc,1,'full');
% scenario 1
p1=sdpvar(n_g,n_t,n_d_c,'full');
x1=sdpvar(n_g,n_t,n_d_c,'full');
u1=sdpvar(n_g,n_t,n_d_c,'full');
alpha1=sdpvar(n_g,n_t,n_d_c,'full');
beta1up=sdpvar(n_g,n_t,n_d_c,'full');
beta1dn=sdpvar(n_g,n_t,n_d_c,'full');
delta1up=sdpvar(n_g,n_t,n_d_c,'full');
delta1dn=sdpvar(n_g,n_t,n_d_c,'full');
gamma1up=sdpvar(n_g,n_t,n_d_c,'full');
gamma1dn=sdpvar(n_g,n_t,n_d_c,'full');
phi1up=sdpvar(n_line,n_t,n_d_c,'full');
phi1dn=sdpvar(n_line,n_t,n_d_c,'full');


%scenario 2
p2=sdpvar(n_g,n_t,n_d_c,'full');
x2=sdpvar(n_g,n_t,n_d_c,'full');
u2=sdpvar(n_g,n_t,n_d_c,'full');
alpha2=sdpvar(n_g,n_t,n_d_c,'full');
beta2_up=sdpvar(n_g,n_t,n_d_c,'full');
beta2_dn=sdpvar(n_g,n_t,n_d_c,'full');
beta2up=sdpvar(n_g,n_t,n_d_c,'full');
beta2dn=sdpvar(n_g,n_t,n_d_c,'full');
delta2up=sdpvar(n_g,n_t,n_d_c,'full');
delta2dn=sdpvar(n_g,n_t,n_d_c,'full');
gamma2up=sdpvar(n_g,n_t,n_d_c,'full');
gamma2dn=sdpvar(n_g,n_t,n_d_c,'full');
phi2up=sdpvar(n_line,n_t,n_d_c,'full');
phi2dn=sdpvar(n_line,n_t,n_d_c,'full');


Constraints=[];
Constraints=[Constraints,u1>=zeros(n_g,n_t,n_d_c),x1>=zeros(n_g,n_t,n_d_c),x1<=ones(n_g,n_t,n_d_c),u2>=zeros(n_g,n_t,n_d_c),x2>=zeros(n_g,n_t,n_d_c),x2<=ones(n_g,n_t,n_d_c),ls1>=0,ls2>=0,ws1>=0,ws2>=0];
% power balance
for d=1:n_d_c
 for t=1:n_t
     % scenario 1
    Constraints=[Constraints,ones(1,n_g)*p1(:,t,d)+ones(1,n_w)*(m_c_n{d,t})==ones(1,n_l)*4*(pd(1,:)')];
    Constraints=[Constraints,ones(1,n_g)*alpha1(:,t,d)==-1];
    %scenario 2
    Constraints=[Constraints,ones(1,n_g)*p2(:,t,d)+ones(1,n_w)*(m_c_n2{d,t})==ones(1,n_l)*0.9*4*(pd(1,:)')];
    Constraints=[Constraints,ones(1,n_g)*alpha2(:,t,d)==-1];
 end
end
%All generating units
for d=1:n_d_c
for t=1:n_t
    for j=1:n_g
         % scenario 1
        % capacity constraints
        Constraints=[Constraints,norm([beta1up(j,t,d);sigma_uni_n{d,t}^(1/2)*alpha1(j,t,d)*ones(n_w,1)],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(pgmax(j)*x1(j,t,d)-p1(j,t,d))-(((2*epsi*(alpha_uni+1))/alpha_uni)-1)*beta1up(j,t,d)];
        Constraints=[Constraints,norm([beta1up(j,t,d);sigma_uni_n{d,t}^(1/2)*alpha1(j,t,d)*ones(n_w,1)],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(pgmax(j)*x1(j,t,d)-p1(j,t,d))-(((2*epsi-1)*(alpha_uni+1)-1)/alpha_uni)*beta1up(j,t,d)];
        Constraints=[Constraints,norm([((alpha_uni+1)/alpha_uni)*beta1up(j,t,d);sigma_uni_n{d,t}^(1/2)^(1/2)*alpha1(j,t,d)*ones(n_w,1)],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(pgmax(j)*x1(j,t,d)-p1(j,t,d))-(((2*epsi-1)*(alpha_uni+1))/alpha_uni)*beta1up(j,t,d)];
         
        Constraints=[Constraints,norm([beta1dn(j,t,d);-sigma_uni_n{d,t}^(1/2)*alpha1(j,t,d)*ones(n_w,1)],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p1(j,t,d)-pgmin(j)*x1(j,t,d))-(((2*epsi*(alpha_uni+1))/alpha_uni)-1)*beta1dn(j,t,d)];
        Constraints=[Constraints,norm([beta1dn(j,t,d);-sigma_uni_n{d,t}^(1/2)*alpha1(j,t,d)*ones(n_w,1)],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p1(j,t,d)-pgmin(j)*x1(j,t,d))-(((2*epsi-1)*(alpha_uni+1)-1)/alpha_uni)*beta1dn(j,t,d)];
        Constraints=[Constraints,norm([((alpha_uni+1)/alpha_uni)*beta1dn(j,t,d);-sigma_uni_n{d,t}^(1/2)^(1/2)*alpha1(j,t,d)*ones(n_w,1)],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p1(j,t,d)-pgmin(j)*x1(j,t,d))-(((2*epsi-1)*(alpha_uni+1))/alpha_uni)*beta1dn(j,t,d)];
  

        % ramp-rate constraints
        Constraints=[Constraints,norm([delta1up(j,t,d);sigma_uni_c_n{d,t}^(1/2)*[alpha1(j,t,d);alpha1(j,t,d);-alpha1(j,t-1,d);-alpha1(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p1(j,t-1,d)-p1(j,t,d)+x1(j,t-1,d)*ru(j)+(1-x1(j,t-1,d))*QS(j))-(((2*epsi*(alpha_uni+1))/alpha_uni)-1)*delta1up(j,t,d)];
        Constraints=[Constraints,norm([delta1up(j,t,d);sigma_uni_c_n{d,t}^(1/2)*[alpha1(j,t,d);alpha1(j,t,d);-alpha1(j,t-1,d);-alpha1(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p1(j,t-1,d)-p1(j,t,d)+x1(j,t-1,d)*ru(j)+(1-x1(j,t-1,d))*QS(j))-(((2*epsi-1)*(alpha_uni+1)-1)/alpha_uni)*delta1up(j,t,d)];
        Constraints=[Constraints,norm([((alpha_uni+1)/alpha_uni)*delta1up(j,t,d);sigma_uni_c_n{d,t}^(1/2)*[alpha1(j,t,d);alpha1(j,t,d);-alpha1(j,t-1,d);-alpha1(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p1(j,t-1,d)-p1(j,t,d)+x1(j,t-1,d)*ru(j)+(1-x1(j,t-1,d))*QS(j))-(((2*epsi-1)*(alpha_uni+1))/alpha_uni)*delta1up(j,t,d)];
        
        Constraints=[Constraints,norm([delta1dn(j,t,d);-sigma_uni_c_n{d,t}^(1/2)*[alpha1(j,t,d);alpha1(j,t,d);-alpha1(j,t-1,d);-alpha1(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p1(j,t,d)-p1(j,t-1,d)+x1(j,t,d)*rd(j)+(1-x1(j,t,d))*QS(j))-(((2*epsi*(alpha_uni+1))/alpha_uni)-1)*delta1dn(j,t,d)];
        Constraints=[Constraints,norm([delta1dn(j,t,d);-sigma_uni_c_n{d,t}^(1/2)*[alpha1(j,t,d);alpha1(j,t,d);-alpha1(j,t-1,d);-alpha1(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p1(j,t,d)-p1(j,t-1,d)+x1(j,t,d)*rd(j)+(1-x1(j,t,d))*QS(j))-(((2*epsi-1)*(alpha_uni+1)-1)/alpha_uni)*delta1dn(j,t,d)];
        Constraints=[Constraints,norm([((alpha_uni+1)/alpha_uni)*delta1dn(j,t,d);-sigma_uni_c_n{d,t}^(1/2)*[alpha1(j,t,d);alpha1(j,t,d);-alpha1(j,t-1,d);-alpha1(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p1(j,t,d)-p1(j,t-1,d)+x1(j,t,d)*rd(j)+(1-x1(j,t,d))*QS(j))-(((2*epsi-1)*(alpha_uni+1))/alpha_uni)*delta1dn(j,t,d)];
       
        Constraints=[Constraints,norm([gamma1up(j,t,d);sigma_uni_c_n{d,t}^(1/2)*[alpha1(j,t,d);alpha1(j,t,d);-alpha1(j,t-1,d);-alpha1(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p1(j,t-1,d)-p1(j,t,d)+(pgmin(j)+ru(j))*x1(j,t,d)-pgmin(j)*x1(j,t-1,d)-(pgmin(j)+ru(j)-QS(j))*u1(j,t,d))-(((2*epsi*(alpha_uni+1))/alpha_uni)-1)*gamma1up(j,t,d)];
        Constraints=[Constraints,norm([gamma1up(j,t,d);sigma_uni_c_n{d,t}^(1/2)*[alpha1(j,t,d);alpha1(j,t,d);-alpha1(j,t-1,d);-alpha1(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p1(j,t-1,d)-p1(j,t,d)+(pgmin(j)+ru(j))*x1(j,t,d)-pgmin(j)*x1(j,t-1,d)-(pgmin(j)+ru(j)-QS(j))*u1(j,t,d))-(((2*epsi-1)*(alpha_uni+1)-1)/alpha_uni)*gamma1up(j,t,d)];
        Constraints=[Constraints,norm([((alpha_uni+1)/alpha_uni)*gamma1up(j,t,d);sigma_uni_c_n{d,t}^(1/2)*[alpha1(j,t,d);alpha1(j,t,d);-alpha1(j,t-1,d);-alpha1(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p1(j,t-1,d)-p1(j,t,d)+(pgmin(j)+ru(j))*x1(j,t,d)-pgmin(j)*x1(j,t-1,d)-(pgmin(j)+ru(j)-QS(j))*u1(j,t,d))-(((2*epsi-1)*(alpha_uni+1))/alpha_uni)*gamma1up(j,t,d)];

        Constraints=[Constraints,norm([gamma1dn(j,t,d);-sigma_uni_c_n{d,t}^(1/2)*[alpha1(j,t,d);alpha1(j,t,d);-alpha1(j,t-1,d);-alpha1(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p1(j,t,d)-p1(j,t-1,d)+QS(j)*x1(j,t-1,d)-(pgmin(j)+rd(j)-QS(j))*u1(j,t,d)-(QS(j)-rd(j))*x1(j,t,d))-(((2*epsi*(alpha_uni+1))/alpha_uni)-1)*gamma1dn(j,t,d)];
        Constraints=[Constraints,norm([gamma1dn(j,t,d);-sigma_uni_c_n{d,t}^(1/2)*[alpha1(j,t,d);alpha1(j,t,d);-alpha1(j,t-1,d);-alpha1(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p1(j,t,d)-p1(j,t-1,d)+QS(j)*x1(j,t-1,d)-(pgmin(j)+rd(j)-QS(j))*u1(j,t,d)-(QS(j)-rd(j))*x1(j,t,d))-(((2*epsi-1)*(alpha_uni+1)-1)/alpha_uni)*gamma1dn(j,t,d)];
        Constraints=[Constraints,norm([((alpha_uni+1)/alpha_uni)*gamma1dn(j,t,d);-sigma_uni_c_n{d,t}^(1/2)*[alpha1(j,t,d);alpha1(j,t,d);-alpha1(j,t-1,d);-alpha1(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p1(j,t,d)-p1(j,t-1,d)+QS(j)*x1(j,t-1,d)-(pgmin(j)+rd(j)-QS(j))*u1(j,t,d)-(QS(j)-rd(j))*x1(j,t,d))-(((2*epsi-1)*(alpha_uni+1))/alpha_uni)*gamma1dn(j,t,d)];
 
            
        Constraints=[Constraints,p1(j,t,d)<=pgmax(j)*x1(j,t,d)];
        Constraints=[Constraints,p1(j,t,d)>=pgmin(j)*x1(j,t,d)]; 

         %scenario 2
        % capacity constraints
        Constraints=[Constraints,norm([beta2up(j,t,d);sigma_uni_n2{d,t}^(1/2)*alpha2(j,t,d)*ones(n_w,1)],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(pgmax(j)*x2(j,t,d)-p2(j,t,d))-(((2*epsi*(alpha_uni+1))/alpha_uni)-1)*beta2up(j,t,d)];
        Constraints=[Constraints,norm([beta2up(j,t,d);sigma_uni_n2{d,t}^(1/2)*alpha2(j,t,d)*ones(n_w,1)],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(pgmax(j)*x2(j,t,d)-p2(j,t,d))-(((2*epsi-1)*(alpha_uni+1)-1)/alpha_uni)*beta2up(j,t,d)];
        Constraints=[Constraints,norm([((alpha_uni+1)/alpha_uni)*beta2up(j,t,d);sigma_uni_n2{d,t}^(1/2)^(1/2)*alpha2(j,t,d)*ones(n_w,1)],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(pgmax(j)*x2(j,t,d)-p2(j,t,d))-(((2*epsi-1)*(alpha_uni+1))/alpha_uni)*beta2up(j,t,d)];
         
        Constraints=[Constraints,norm([beta2dn(j,t,d);-sigma_uni_n2{d,t}^(1/2)*alpha2(j,t,d)*ones(n_w,1)],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p2(j,t,d)-pgmin(j)*x2(j,t,d))-(((2*epsi*(alpha_uni+1))/alpha_uni)-1)*beta2dn(j,t,d)];
        Constraints=[Constraints,norm([beta2dn(j,t,d);-sigma_uni_n2{d,t}^(1/2)*alpha2(j,t,d)*ones(n_w,1)],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p2(j,t,d)-pgmin(j)*x2(j,t,d))-(((2*epsi-1)*(alpha_uni+1)-1)/alpha_uni)*beta2dn(j,t,d)];
        Constraints=[Constraints,norm([((alpha_uni+1)/alpha_uni)*beta2dn(j,t,d);-sigma_uni_n2{d,t}^(1/2)^(1/2)*alpha2(j,t,d)*ones(n_w,1)],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p2(j,t,d)-pgmin(j)*x2(j,t,d))-(((2*epsi-1)*(alpha_uni+1))/alpha_uni)*beta2dn(j,t,d)];
  

        % ramp-rate constraints
        Constraints=[Constraints,norm([delta2up(j,t,d);sigma_uni_c_n2{d,t}^(1/2)*[alpha2(j,t,d);alpha2(j,t,d);-alpha2(j,t-1,d);-alpha2(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p2(j,t-1,d)-p2(j,t,d)+x2(j,t-1,d)*ru(j)+(1-x2(j,t-1,d))*QS(j))-(((2*epsi*(alpha_uni+1))/alpha_uni)-1)*delta2up(j,t,d)];
        Constraints=[Constraints,norm([delta2up(j,t,d);sigma_uni_c_n2{d,t}^(1/2)*[alpha2(j,t,d);alpha2(j,t,d);-alpha2(j,t-1,d);-alpha2(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p2(j,t-1,d)-p2(j,t,d)+x2(j,t-1,d)*ru(j)+(1-x2(j,t-1,d))*QS(j))-(((2*epsi-1)*(alpha_uni+1)-1)/alpha_uni)*delta2up(j,t,d)];
        Constraints=[Constraints,norm([((alpha_uni+1)/alpha_uni)*delta2up(j,t,d);sigma_uni_c_n2{d,t}^(1/2)*[alpha2(j,t,d);alpha2(j,t,d);-alpha2(j,t-1,d);-alpha2(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p2(j,t-1,d)-p2(j,t,d)+x2(j,t-1,d)*ru(j)+(1-x2(j,t-1,d))*QS(j))-(((2*epsi-1)*(alpha_uni+1))/alpha_uni)*delta2up(j,t,d)];
        
        Constraints=[Constraints,norm([delta2dn(j,t,d);-sigma_uni_c_n2{d,t}^(1/2)*[alpha2(j,t,d);alpha2(j,t,d);-alpha2(j,t-1,d);-alpha2(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p2(j,t,d)-p2(j,t-1,d)+x2(j,t,d)*rd(j)+(1-x2(j,t,d))*QS(j))-(((2*epsi*(alpha_uni+1))/alpha_uni)-1)*delta2dn(j,t,d)];
        Constraints=[Constraints,norm([delta2dn(j,t,d);-sigma_uni_c_n2{d,t}^(1/2)*[alpha2(j,t,d);alpha2(j,t,d);-alpha2(j,t-1,d);-alpha2(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p2(j,t,d)-p2(j,t-1,d)+x2(j,t,d)*rd(j)+(1-x2(j,t,d))*QS(j))-(((2*epsi-1)*(alpha_uni+1)-1)/alpha_uni)*delta2dn(j,t,d)];
        Constraints=[Constraints,norm([((alpha_uni+1)/alpha_uni)*delta2dn(j,t,d);-sigma_uni_c_n2{d,t}^(1/2)*[alpha2(j,t,d);alpha2(j,t,d);-alpha2(j,t-1,d);-alpha2(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p2(j,t,d)-p2(j,t-1,d)+x2(j,t,d)*rd(j)+(1-x2(j,t,d))*QS(j))-(((2*epsi-1)*(alpha_uni+1))/alpha_uni)*delta2dn(j,t,d)];
       
        Constraints=[Constraints,norm([gamma2up(j,t,d);sigma_uni_c_n2{d,t}^(1/2)*[alpha2(j,t,d);alpha2(j,t,d);-alpha2(j,t-1,d);-alpha2(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p2(j,t-1,d)-p2(j,t,d)+(pgmin(j)+ru(j))*x2(j,t,d)-pgmin(j)*x2(j,t-1,d)-(pgmin(j)+ru(j)-QS(j))*u2(j,t,d))-(((2*epsi*(alpha_uni+1))/alpha_uni)-1)*gamma2up(j,t,d)];
        Constraints=[Constraints,norm([gamma2up(j,t,d);sigma_uni_c_n2{d,t}^(1/2)*[alpha2(j,t,d);alpha2(j,t,d);-alpha2(j,t-1,d);-alpha2(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p2(j,t-1,d)-p2(j,t,d)+(pgmin(j)+ru(j))*x2(j,t,d)-pgmin(j)*x2(j,t-1,d)-(pgmin(j)+ru(j)-QS(j))*u2(j,t,d))-(((2*epsi-1)*(alpha_uni+1)-1)/alpha_uni)*gamma2up(j,t,d)];
        Constraints=[Constraints,norm([((alpha_uni+1)/alpha_uni)*gamma2up(j,t,d);sigma_uni_c_n2{d,t}^(1/2)*[alpha2(j,t,d);alpha2(j,t,d);-alpha2(j,t-1,d);-alpha2(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p2(j,t-1,d)-p2(j,t,d)+(pgmin(j)+ru(j))*x2(j,t,d)-pgmin(j)*x2(j,t-1,d)-(pgmin(j)+ru(j)-QS(j))*u2(j,t,d))-(((2*epsi-1)*(alpha_uni+1))/alpha_uni)*gamma2up(j,t,d)];

        Constraints=[Constraints,norm([gamma2dn(j,t,d);-sigma_uni_c_n2{d,t}^(1/2)*[alpha2(j,t,d);alpha2(j,t,d);-alpha2(j,t-1,d);-alpha2(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p2(j,t,d)-p2(j,t-1,d)+QS(j)*x2(j,t-1,d)-(pgmin(j)+rd(j)-QS(j))*u2(j,t,d)-(QS(j)-rd(j))*x2(j,t,d))-(((2*epsi*(alpha_uni+1))/alpha_uni)-1)*gamma2dn(j,t,d)];
        Constraints=[Constraints,norm([gamma2dn(j,t,d);-sigma_uni_c_n2{d,t}^(1/2)*[alpha2(j,t,d);alpha2(j,t,d);-alpha2(j,t-1,d);-alpha2(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p2(j,t,d)-p2(j,t-1,d)+QS(j)*x2(j,t-1,d)-(pgmin(j)+rd(j)-QS(j))*u2(j,t,d)-(QS(j)-rd(j))*x2(j,t,d))-(((2*epsi-1)*(alpha_uni+1)-1)/alpha_uni)*gamma2dn(j,t,d)];
        Constraints=[Constraints,norm([((alpha_uni+1)/alpha_uni)*gamma2dn(j,t,d);-sigma_uni_c_n2{d,t}^(1/2)*[alpha2(j,t,d);alpha2(j,t,d);-alpha2(j,t-1,d);-alpha2(j,t-1,d)]],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(p2(j,t,d)-p2(j,t-1,d)+QS(j)*x2(j,t-1,d)-(pgmin(j)+rd(j)-QS(j))*u2(j,t,d)-(QS(j)-rd(j))*x2(j,t,d))-(((2*epsi-1)*(alpha_uni+1))/alpha_uni)*gamma2dn(j,t,d)];


        
        Constraints=[Constraints,p2(j,t,d)<=pgmax(j)*x2(j,t,d)];
        Constraints=[Constraints,p2(j,t,d)>=pgmin(j)*x2(j,t,d)];                 

%        
    end
    for j=1:n_line
    %transmission lines
     %scenario 1
        Constraints=[Constraints,norm([phi1up(j,t,d);sigma_uni_n{d,t}^(1/2)*(Hw(j,:)+Hg(j,1:n_g)*alpha1(:,t,d)*ones(1,n_w))'],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(10*fmax(j)-Hg(j,1:n_g)*p1(:,t,d)-Hw(j,:)*m_c_n{d,t}+Hl(j,:)*(4*pd'))-(((2*epsi*(alpha_uni+1))/alpha_uni)-1)*phi1up(j,t,d)];
        Constraints=[Constraints,norm([phi1up(j,t,d);sigma_uni_n{d,t}^(1/2)*(Hw(j,:)+Hg(j,1:n_g)*alpha1(:,t,d)*ones(1,n_w))'],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(10*fmax(j)-Hg(j,1:n_g)*p1(:,t,d)-Hw(j,:)*m_c_n{d,t}+Hl(j,:)*(4*pd'))-(((2*epsi-1)*(alpha_uni+1)-1)/alpha_uni)*phi1up(j,t,d)];
        Constraints=[Constraints,norm([((alpha_uni+1)/alpha_uni)*phi1up(j,t,d);sigma_uni_n{d,t}^(1/2)*(Hw(j,:)+Hg(j,1:n_g)*alpha1(:,t,d)*ones(1,n_w))'],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(10*fmax(j)-Hg(j,1:n_g)*p1(:,t,d)-Hw(j,:)*m_c_n{d,t}+Hl(j,:)*(4*pd'))-(((2*epsi-1)*(alpha_uni+1))/alpha_uni)*phi1up(j,t,d)];
       
        Constraints=[Constraints,norm([phi1dn(j,t,d);-sigma_uni_n{d,t}^(1/2)*(Hw(j,:)+Hg(j,1:n_g)*alpha1(:,t,d)*ones(1,n_w))'],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(10*fmax(j)+Hg(j,1:n_g)*p1(:,t,d)+Hw(j,:)*m_c_n{d,t}-Hl(j,:)*(4*pd'))-(((2*epsi*(alpha_uni+1))/alpha_uni)-1)*phi1dn(j,t,d)];
        Constraints=[Constraints,norm([phi1dn(j,t,d);-sigma_uni_n{d,t}^(1/2)*(Hw(j,:)+Hg(j,1:n_g)*alpha1(:,t,d)*ones(1,n_w))'],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(10*fmax(j)+Hg(j,1:n_g)*p1(:,t,d)+Hw(j,:)*m_c_n{d,t}-Hl(j,:)*(4*pd'))-(((2*epsi-1)*(alpha_uni+1)-1)/alpha_uni)*phi1dn(j,t,d)];
        Constraints=[Constraints,norm([((alpha_uni+1)/alpha_uni)*phi1dn(j,t,d);-sigma_uni_n{d,t}^(1/2)*(Hw(j,:)+Hg(j,1:n_g)*alpha1(:,t,d)*ones(1,n_w))'],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(10*fmax(j)+Hg(j,1:n_g)*p1(:,t,d)+Hw(j,:)*m_c_n{d,t}-Hl(j,:)*(4*pd'))-(((2*epsi-1)*(alpha_uni+1))/alpha_uni)*phi1dn(j,t,d)];
    
     %scenario 2
        Constraints=[Constraints,norm([phi2up(j,t,d);sigma_uni_n2{d,t}^(1/2)*(Hw(j,:)+Hg(j,1:n_g)*alpha2(:,t,d)*ones(1,n_w))'],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(10*fmax(j)-Hg(j,1:n_g)*p2(:,t,d)-Hw(j,:)*m_c_n{d,t}+Hl(j,:)*(4*0.9*pd'))-(((2*epsi*(alpha_uni+1))/alpha_uni)-1)*phi2up(j,t,d)];
        Constraints=[Constraints,norm([phi2up(j,t,d);sigma_uni_n2{d,t}^(1/2)*(Hw(j,:)+Hg(j,1:n_g)*alpha2(:,t,d)*ones(1,n_w))'],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(10*fmax(j)-Hg(j,1:n_g)*p2(:,t,d)-Hw(j,:)*m_c_n{d,t}+Hl(j,:)*(4*0.9*pd'))-(((2*epsi-1)*(alpha_uni+1)-1)/alpha_uni)*phi2up(j,t,d)];
        Constraints=[Constraints,norm([((alpha_uni+1)/alpha_uni)*phi2up(j,t,d);sigma_uni_n2{d,t}^(1/2)*(Hw(j,:)+Hg(j,1:n_g)*alpha2(:,t,d)*ones(1,n_w))'],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(10*fmax(j)-Hg(j,1:n_g)*p2(:,t,d)-Hw(j,:)*m_c_n{d,t}+Hl(j,:)*(4*0.9*pd'))-(((2*epsi-1)*(alpha_uni+1))/alpha_uni)*phi2up(j,t,d)];
       
        Constraints=[Constraints,norm([phi2dn(j,t,d);-sigma_uni_n2{d,t}^(1/2)*(Hw(j,:)+Hg(j,1:n_g)*alpha2(:,t,d)*ones(1,n_w))'],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(10*fmax(j)+Hg(j,1:n_g)*p2(:,t,d)+Hw(j,:)*m_c_n{d,t}-Hl(j,:)*(4*0.9*pd'))-(((2*epsi*(alpha_uni+1))/alpha_uni)-1)*phi2dn(j,t,d)];
        Constraints=[Constraints,norm([phi2dn(j,t,d);-sigma_uni_n2{d,t}^(1/2)*(Hw(j,:)+Hg(j,1:n_g)*alpha2(:,t,d)*ones(1,n_w))'],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(10*fmax(j)+Hg(j,1:n_g)*p2(:,t,d)+Hw(j,:)*m_c_n{d,t}-Hl(j,:)*(4*0.9*pd'))-(((2*epsi-1)*(alpha_uni+1)-1)/alpha_uni)*phi2dn(j,t,d)];
        Constraints=[Constraints,norm([((alpha_uni+1)/alpha_uni)*phi2dn(j,t,d);-sigma_uni_n2{d,t}^(1/2)*(Hw(j,:)+Hg(j,1:n_g)*alpha2(:,t,d)*ones(1,n_w))'],2)<=((2*epsi*(alpha_uni+1))/alpha_uni)*(10*fmax(j)+Hg(j,1:n_g)*p2(:,t,d)+Hw(j,:)*m_c_n{d,t}-Hl(j,:)*(4*0.9*pd'))-(((2*epsi-1)*(alpha_uni+1))/alpha_uni)*phi2dn(j,t,d)];    
    end
end
end

%investment constraints
for d=1:n_d_c
for t=1:n_t
    for j=n_ge+1:n_g
        % scenario 1
        Constraints=[Constraints,x1(j,t,d)<=v(j-n_ge)];
        %scenario 2
        Constraints=[Constraints,x2(j,t,d)<=v(j-n_ge)];
    end
end
end

for d=1:n_d_c
for t=1+1:n_t
    for j=1:n_g
        %scenario 1
        Constraints=[Constraints,u1(j,t,d)>=x1(j,t,d)-x1(j,t-1,d)];
        %scneario 2
        Constraints=[Constraints,u2(j,t,d)>=x2(j,t,d)-x2(j,t-1,d)];
    end
end
end

% min. up/dn time
for d=1:n_d_c
for j=1:n_g
    for t=lu(j)+1:n_t
        Constraints=[Constraints,u1(j,t-lu(j)+1:t,d)*ones(length(t-lu(j)+1:t),1)<=x1(j,t,d)];
        Constraints=[Constraints,u2(j,t-lu(j)+1:t,d)*ones(length(t-lu(j)+1:t),1)<=x2(j,t,d)];

    end
    for t=ld(j)+1:n_t
        Constraints=[Constraints,u1(j,t-ld(j)+1:t,d)*ones(length(t-ld(j)+1:t),1)<=1-x1(j,t-ld(j),d)];
        Constraints=[Constraints,u2(j,t-ld(j)+1:t,d)*ones(length(t-ld(j)+1:t),1)<=1-x2(j,t-ld(j),d)];
    end
end
end

        


% objective function

Objective=0;
cost_base_operation=0;
for j=n_ge+1:n_g
Objective=(cost_inv(j-n_ge))*v(j-n_ge)+Objective;
end
for d=1:n_d_c
for t=1:n_t

    Objective=Objective+weigth_c(d,1)*(0.5*(cost_op(1:n_g)'*p1(1:n_g,t,d)+(cost_st(1:n_g))'*u1(1:n_g,t,d))+0.5*(cost_op(1:n_g)'*p2(1:n_g,t,d)+(cost_st(1:n_g))'*u2(1:n_g,t,d)));
    cost_base_operation=cost_base_operation+weigth_c(d,1)*(0.5*(cost_op(1:n_g)'*p1(1:n_g,t,d)+(cost_st(1:n_g))'*u1(1:n_g,t,d))+0.5*(cost_op(1:n_g)'*p2(1:n_g,t,d)+(cost_st(1:n_g))'*u2(1:n_g,t,d)));

end
end

ops = sdpsettings('solver','gurobi','verbose',0);
sol=optimize(Constraints,Objective,ops);
% 
v_s=value(v);
p1_s=value(p1);
x1_s=value(x1);
u1_s=value(u1);
alpha1_s=value(alpha1);
ls1_s=value(ls1);
ws1_s=value(ws1);
p2_s=value(p2);
x2_s=value(x2);
u2_s=value(u2);
alpha2_s=value(alpha2);
ls2_s=value(ls2);
ws2_s=value(ws2);





time=toc

