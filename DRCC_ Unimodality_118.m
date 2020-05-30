function [infeas,alpha1v,p1v,x1v,alpha2v,p2v,x2v,u1v,u2v,v_v,result_ct,result_co]=OPF_118_unimodality(eta1up,eta1dn,pi1up,delta1up,delta1dn,gamma1up,gamma1dn,phi1up,phi1dn,eta2up,eta2dn,pi2up,delta2up,delta2dn,gamma2up,gamma2dn,phi2up,phi2dn)

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

 
epsi=0.25;

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
%scenario 2
p2=sdpvar(n_g,n_t,n_d_c,'full');
x2=sdpvar(n_g,n_t,n_d_c,'full');
u2=sdpvar(n_g,n_t,n_d_c,'full');
alpha2=sdpvar(n_g,n_t,n_d_c,'full');


Constraints=[];



Constraints=[Constraints,u1>=zeros(n_g,n_t,n_d_c),x1>=zeros(n_g,n_t,n_d_c),x1<=ones(n_g,n_t,n_d_c),u2>=zeros(n_g,n_t,n_d_c),x2>=zeros(n_g,n_t,n_d_c),x2<=ones(n_g,n_t,n_d_c),ls1>=0,ls2>=0,ws1>=0,ws2>=0];
% Ppower balance
for d=1:n_d_c
 for t=1:n_t
     % scenario 1
    Constraints=[Constraints,ones(1,n_g)*p1(:,t,d)+ones(1,n_w)*(m_c_n{d,t})+ls1-ws1==ones(1,n_l)*4*(pd(1,:)')];
    Constraints=[Constraints,ones(1,n_g)*alpha1(:,t,d)==-1];
%     %scenario 2
    Constraints=[Constraints,ones(1,n_g)*p2(:,t,d)+ones(1,n_w)*(m_c_n2{d,t})+ls2-ws2==ones(1,n_l)*0.9*4*(pd(1,:)')];
    Constraints=[Constraints,ones(1,n_g)*alpha2(:,t,d)==-1];
 end
end

for d=1:n_d_c
for t=1:n_t
   
    for j=1:n_g

        % scenario 1
        for k=1:length(eta1up)
        if eta1up(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(eta1up(k)^(-alpha_uni)))/epsi)*norm(alpha1(j,t,d)*ones(1,n_w)*sigma_uni_n{d,t}^(1/2),2)<=eta1up(k)*(pgmax(j)*x1(j,t,d)-p1(j,t,d)-alpha1(j,t,d)*ones(1,n_w)*mode_error_n{d,t})-((alpha_uni+1)/alpha_uni)*(alpha1(j,t,d)*ones(1,n_w)*(0-mode_error_n{d,t}))];
        end
        end
        for k=1:length(eta1dn)
        if eta1dn(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(eta1dn(k)^(-alpha_uni)))/epsi)*norm(-alpha1(j,t,d)*ones(1,n_w)*sigma_uni_n{d,t}^(1/2),2)<=eta1dn(k)*(p1(j,t,d)-pgmin(j)*x1(j,t,d)+alpha1(j,t,d)*ones(1,n_w)*mode_error_n{d,t})+((alpha_uni+1)/alpha_uni)*(alpha1(j,t,d)*ones(1,n_w)*(0-mode_error_n{d,t}))];   
        end
        end
        for k=1:length(pi1up)
        if pi1up(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(pi1up(k)^(-alpha_uni)))/epsi)*norm(alpha1(j,t-1,d)*ones(1,n_w)*sigma_uni_n{d,t}^(1/2),2)<=pi1up(k)*(QS(j)*x1(j,t-1,d)+(pgmax(j)-QS(j))*(x1(j,t,d)-u1(j,t,d))-p1(j,t-1,d)-alpha1(j,t-1,d)*ones(1,n_w)*mode_error_n{d,t})-((alpha_uni+1)/alpha_uni)*(alpha1(j,t-1,d)*ones(1,n_w)*(0-mode_error_n{d,t}))];
        end
        end
        % ramp-rate constraints
         for k=1:length(delta1up)
        if delta1up(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(delta1up(k)^(-alpha_uni)))/epsi)*norm([alpha1(j,t,d) alpha1(j,t,d) -alpha1(j,t-1,d) -alpha1(j,t-1,d)]*sigma_uni_c_n{d,t}^(1/2),2)<=delta1up(k)*(p1(j,t-1,d)-p1(j,t,d)+x1(j,t-1,d)*ru(j)+(1-x1(j,t-1,d))*QS(j)-[alpha1(j,t,d) alpha1(j,t,d) -alpha1(j,t-1,d) -alpha1(j,t-1,d)]*mode_error_c_n{d,t})-((alpha_uni+1)/alpha_uni)*([alpha1(j,t,d) alpha1(j,t,d) -alpha1(j,t-1,d) -alpha1(j,t-1,d)]*(0-mode_error_n{d,t}))];
        end
        end
        for k=1:length(delta1dn)
        if delta1dn(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(delta1dn(k)^(-alpha_uni)))/epsi)*norm([alpha1(j,t,d) alpha1(j,t,d) -alpha1(j,t-1,d) -alpha1(j,t-1,d)]*sigma_uni_c_n{d,t}^(1/2),2)<=delta1dn(k)*(p1(j,t,d)-p1(j,t-1,d)+x1(j,t,d)*rd(j)+(1-x1(j,t,d))*QS(j)+[alpha1(j,t,d) alpha1(j,t,d) -alpha1(j,t-1,d) -alpha1(j,t-1,d)]*mode_error_c_n{d,t})+((alpha_uni+1)/alpha_uni)*([alpha1(j,t,d) alpha1(j,t,d) -alpha1(j,t-1,d) -alpha1(j,t-1,d)]*(0-mode_error_n{d,t}))];     
        end
        end
        for k=1:length(gamma1up)
        if gamma1up(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(gamma1up(k)^(-alpha_uni)))/epsi)*norm([alpha1(j,t,d) alpha1(j,t,d) -alpha1(j,t-1,d) -alpha1(j,t-1,d)]*sigma_uni_c_n{d,t}^(1/2),2)<=gamma1up(k)*(p1(j,t-1,d)-p1(j,t,d)+(pgmin(j)+ru(j))*x1(j,t,d)-pgmin(j)*x1(j,t-1,d)-(pgmin(j)+ru(j)-QS(j))*u1(j,t,d)-[alpha1(j,t,d) alpha1(j,t,d) -alpha1(j,t-1,d) -alpha1(j,t-1,d)]*mode_error_c_n{d,t})-((alpha_uni+1)/alpha_uni)*([alpha1(j,t,d) alpha1(j,t,d) -alpha1(j,t-1,d) -alpha1(j,t-1,d)]*(0-mode_error_n{d,t}))];     
        end
        end
        for k=1:length(gamma1dn)
        if gamma1dn(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(gamma1dn(k)^(-alpha_uni)))/epsi)*norm([alpha1(j,t,d) alpha1(j,t,d) -alpha1(j,t-1,d) -alpha1(j,t-1,d)]*sigma_uni_c_n{d,t}^(1/2),2)<=gamma1dn(k)*(p1(j,t,d)-p1(j,t-1,d)+QS(j)*x1(j,t-1,d)-(pgmin(j)+rd(j)-QS(j))*u1(j,t,d)-(QS(j)-rd(j))*x1(j,t,d)+[alpha1(j,t,d) alpha1(j,t,d) -alpha1(j,t-1,d) -alpha1(j,t-1,d)]*mode_error_c_n{d,t})+((alpha_uni+1)/alpha_uni)*([alpha1(j,t,d) alpha1(j,t,d) -alpha1(j,t-1,d) -alpha1(j,t-1,d)]*(0-mode_error_n{d,t}))];     
        end
        end
    
        Constraints=[Constraints,p1(j,t,d)<=pgmax(j)*x1(j,t,d)];
        Constraints=[Constraints,p1(j,t,d)>=pgmin(j)*x1(j,t,d)]; 
    
        % scenario 2
        for k=1:length(eta2up)
        if eta2up(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(eta2up(k)^(-alpha_uni)))/epsi)*norm(alpha2(j,t,d)*ones(1,n_w)*sigma_uni_n2{d,t}^(1/2),2)<=eta2up(k)*(pgmax(j)*x2(j,t,d)-p2(j,t,d)-alpha2(j,t,d)*ones(1,n_w)*mode_error_n2{d,t})-((alpha_uni+1)/alpha_uni)*(alpha2(j,t,d)*ones(1,n_w)*(0-mode_error_n2{d,t}))];
        end
        end
        for k=1:length(eta2dn)
        if eta2dn(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(eta2dn(k)^(-alpha_uni)))/epsi)*norm(-alpha2(j,t,d)*ones(1,n_w)*sigma_uni_n2{d,t}^(1/2),2)<=eta2dn(k)*(p2(j,t,d)-pgmin(j)*x2(j,t,d)+alpha2(j,t,d)*ones(1,n_w)*mode_error_n2{d,t})+((alpha_uni+1)/alpha_uni)*(alpha2(j,t,d)*ones(1,n_w)*(0-mode_error_n2{d,t}))];   
        end
        end
        for k=1:length(pi2up)
        if pi2up(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(pi2up(k)^(-alpha_uni)))/epsi)*norm(alpha2(j,t-1,d)*ones(1,n_w)*sigma_uni_n2{d,t}^(1/2),2)<=pi2up(k)*(QS(j)*x2(j,t-1,d)+(pgmax(j)-QS(j))*(x2(j,t,d)-u2(j,t,d))-p2(j,t-1,d)-alpha2(j,t-1,d)*ones(1,n_w)*mode_error_n2{d,t})-((alpha_uni+1)/alpha_uni)*(alpha2(j,t-1,d)*ones(1,n_w)*(0-mode_error_n2{d,t}))];
        end
        end
        % ramp-rate constraints
         for k=1:length(delta2up)
        if delta2up(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(delta2up(k)^(-alpha_uni)))/epsi)*norm([alpha2(j,t,d) alpha2(j,t,d) -alpha2(j,t-1,d) -alpha2(j,t-1,d)]*sigma_uni_c_n2{d,t}^(1/2),2)<=delta2up(k)*(p2(j,t-1,d)-p2(j,t,d)+x2(j,t-1,d)*ru(j)+(1-x2(j,t-1,d))*QS(j)-[alpha2(j,t,d) alpha2(j,t,d) -alpha2(j,t-1,d) -alpha2(j,t-1,d)]*mode_error_c_n{d,t})-((alpha_uni+1)/alpha_uni)*([alpha2(j,t,d) alpha2(j,t,d) -alpha2(j,t-1,d) -alpha2(j,t-1,d)]*(0-mode_error_n2{d,t}))];
        end
        end
        for k=1:length(delta2dn)
        if delta2dn(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(delta2dn(k)^(-alpha_uni)))/epsi)*norm([alpha2(j,t,d) alpha2(j,t,d) -alpha2(j,t-1,d) -alpha2(j,t-1,d)]*sigma_uni_c_n2{d,t}^(1/2),2)<=delta2dn(k)*(p2(j,t,d)-p2(j,t-1,d)+x2(j,t,d)*rd(j)+(1-x2(j,t,d))*QS(j)+[alpha2(j,t,d) alpha2(j,t,d) -alpha2(j,t-1,d) -alpha2(j,t-1,d)]*mode_error_c_n{d,t})+((alpha_uni+1)/alpha_uni)*([alpha2(j,t,d) alpha2(j,t,d) -alpha2(j,t-1,d) -alpha2(j,t-1,d)]*(0-mode_error_n2{d,t}))];     
        end
        end
        for k=1:length(gamma2up)
        if gamma2up(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(gamma2up(k)^(-alpha_uni)))/epsi)*norm([alpha2(j,t,d) alpha2(j,t,d) -alpha2(j,t-1,d) -alpha2(j,t-1,d)]*sigma_uni_c_n2{d,t}^(1/2),2)<=gamma2up(k)*(p2(j,t-1,d)-p2(j,t,d)+(pgmin(j)+ru(j))*x2(j,t,d)-pgmin(j)*x2(j,t-1,d)-(pgmin(j)+ru(j)-QS(j))*u2(j,t,d)-[alpha2(j,t,d) alpha2(j,t,d) -alpha2(j,t-1,d) -alpha2(j,t-1,d)]*mode_error_c_n{d,t})-((alpha_uni+1)/alpha_uni)*([alpha2(j,t,d) alpha2(j,t,d) -alpha2(j,t-1,d) -alpha2(j,t-1,d)]*(0-mode_error_n2{d,t}))];     
        end
        end
        for k=1:length(gamma2dn)
        if gamma2dn(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(gamma2dn(k)^(-alpha_uni)))/epsi)*norm([alpha2(j,t,d) alpha2(j,t,d) -alpha2(j,t-1,d) -alpha2(j,t-1,d)]*sigma_uni_c_n2{d,t}^(1/2),2)<=gamma2dn(k)*(p2(j,t,d)-p2(j,t-1,d)+QS(j)*x2(j,t-1,d)-(pgmin(j)+rd(j)-QS(j))*u2(j,t,d)-(QS(j)-rd(j))*x2(j,t,d)+[alpha2(j,t,d) alpha2(j,t,d) -alpha2(j,t-1,d) -alpha2(j,t-1,d)]*mode_error_c_n{d,t})+((alpha_uni+1)/alpha_uni)*([alpha2(j,t,d) alpha2(j,t,d) -alpha2(j,t-1,d) -alpha2(j,t-1,d)]*(0-mode_error_n2{d,t}))];     
        end
        end
      
        Constraints=[Constraints,p2(j,t,d)<=pgmax(j)*x2(j,t,d)];
        Constraints=[Constraints,p2(j,t,d)>=pgmin(j)*x2(j,t,d)];   
        
    end
    
    
    for j=1:n_line
        % transmission lines
        % scenario 1
        for k=1:length(phi1up)
        if phi1up(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(phi1up(k)^(-alpha_uni)))/epsi)*norm((Hw(j,:)+Hg(j,1:n_g)*alpha1(:,t,d)*ones(1,n_w))*sigma_uni_n{d,t}^(1/2),2)<=phi1up(k)*(10*fmax(j)-Hg(j,1:n_g)*p1(:,t,d)-Hw(j,:)*m_c_n{d,t}+Hl(j,:)*(4*pd')-(Hw(j,:)+Hg(j,1:n_g)*alpha1(:,t,d)*ones(1,n_w))*mode_error_n{d,t})-((alpha_uni+1)/alpha_uni)*((Hw(j,:)+Hg(j,1:n_g)*alpha1(:,t,d)*ones(1,n_w))*(0-mode_error_n{d,t}))];
        end
        end
        for k=1:length(phi1dn)
        if phi1dn(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(phi1dn(k)^(-alpha_uni)))/epsi)*norm((Hw(j,:)+Hg(j,1:n_g)*alpha1(:,t,d)*ones(1,n_w))*sigma_uni_n{d,t}^(1/2),2)<=phi1dn(k)*(10*fmax(j)+Hg(j,1:n_g)*p1(:,t,d)+Hw(j,:)*m_c_n{d,t}-Hl(j,:)*(4*pd')+(Hw(j,:)+Hg(j,1:n_g)*alpha1(:,t,d)*ones(1,n_w))*mode_error_n{d,t})+((alpha_uni+1)/alpha_uni)*((Hw(j,:)+Hg(j,1:n_g)*alpha1(:,t,d)*ones(1,n_w))*(0-mode_error_n{d,t}))];
        end
        end
        % scenario 2
        for k=1:length(phi2up)
        if phi12up(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(phi2up(k)^(-alpha_uni)))/epsi)*norm((Hw(j,:)+Hg(j,1:n_g)*alpha2(:,t,d)*ones(1,n_w))*sigma_uni_n2{d,t}^(1/2),2)<=phi2up(k)*(10*fmax(j)-Hg(j,1:n_g)*p2(:,t,d)-Hw(j,:)*m_c_n{d,t}+Hl(j,:)*(4*pd')-(Hw(j,:)+Hg(j,1:n_g)*alpha2(:,t,d)*ones(1,n_w))*mode_error_n2{d,t})-((alpha_uni+1)/alpha_uni)*((Hw(j,:)+Hg(j,1:n_g)*alpha2(:,t,d)*ones(1,n_w))*(0-mode_error_n2{d,t}))];
        end
        end
        for k=1:length(phi2dn)
        if phi12dn(k)>0
        Constraints=[Constraints,sqrt((1-epsi-(phi2dn(k)^(-alpha_uni)))/epsi)*norm((Hw(j,:)+Hg(j,1:n_g)*alpha2(:,t,d)*ones(1,n_w))*sigma_uni_n2{d,t}^(1/2),2)<=phi2dn(k)*(10*fmax(j)+Hg(j,1:n_g)*p2(:,t,d)+Hw(j,:)*m_c_n{d,t}-Hl(j,:)*(4*pd')+(Hw(j,:)+Hg(j,1:n_g)*alpha2(:,t,d)*ones(1,n_w))*mode_error_n2{d,t})+((alpha_uni+1)/alpha_uni)*((Hw(j,:)+Hg(j,1:n_g)*alpha2(:,t,d)*ones(1,n_w))*(0-mode_error_n2{d,t}))];
        end
        end
    end
        
end
end


%investment constraints
for d=1:n_d_c
for t=1:n_t
    for j=n_ge+1:n_g
        %scenario 1
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
infeas=sol.problem;
% 
v_v=value(v);
p1v=value(p1);
x1v=value(x1);
u1v=value(u1);
alpha1v=value(alpha1);

p2v=value(p2);
x2v=value(x2);
u2v=value(u2);
alpha2v=value(alpha2);


result_ct=value(Objective);
result_co=value(cost_base_operation);




time=toc

end