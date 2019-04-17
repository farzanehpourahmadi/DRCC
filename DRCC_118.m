%load('parameters_118bus.mat')
tic
n_w=2; %number of m plants
n_ge=19; %number of existing generators
n_gc=22; %number of candidate generators
n_g=n_gc+n_ge; %number of total generators
n_l=99; %number of loads
n_line=186; %number of transmission lines
n_d_c=10; %number of representative days
n_t=24;   %number of operating hours                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  00;
mu(n_t,n_w)=0;
epsi=0.1;
c_r=1;
% Hg=PTDF*busgen; % sensitivity of each line to each generation unit
% Hw=PTDF*buswind; % sensitivity of each line to each wind farm
% Hl=PTDF*busload;  % sensitivity of each line to each load

 k=1;
 for d=1:n_d_c
 for t=1:n_t
     m_c_n{d,t}=3*m_c{d,t};
     cov_w_corr_c_n{d,t}=3*3*cov_w_corr_c{d,t};
     cov_corr_c_n{d,t}=3*3*cov_corr_c{d,t};

     m_c_n2{d,t}=3*0.9*m_c{d,t};
     cov_w_corr_c_n2{d,t}=3*3*0.9*0.9*cov_w_corr_c{d,t};
     cov_corr_c_n2{d,t}=3*3*0.9*0.9*cov_corr_c{d,t};
 end
end
% definig variables
p1=sdpvar(n_g,n_t,n_d_c,'full');
x1=sdpvar(n_g,n_t,n_d_c,'full');
u1=sdpvar(n_g,n_t,n_d_c,'full');
alpha1=sdpvar(n_g,n_t,n_d_c,'full');
v=binvar(n_gc,1,'full');
%scenario 2
p2=sdpvar(n_g,n_t,n_d_c,'full');
x2=sdpvar(n_g,n_t,n_d_c,'full');
u2=sdpvar(n_g,n_t,n_d_c,'full');
alpha2=sdpvar(n_g,n_t,n_d_c,'full');


for c_epsi=0.2:0.2
Constraints=[];
Constraints=[Constraints,u1>=zeros(n_g,n_t,n_d_c),x1>=zeros(n_g,n_t,n_d_c),x1<=ones(n_g,n_t,n_d_c),u2>=zeros(n_g,n_t,n_d_c),x2>=zeros(n_g,n_t,n_d_c),x2<=ones(n_g,n_t,n_d_c)];
for d=1:n_d_c
 for t=1:n_t
    Constraints=[Constraints,ones(1,n_g)*p1(:,t,d)+ones(1,n_w)*(m_c_n{d,t})==ones(1,n_l)*4*(pd(1,:)')];
    Constraints=[Constraints,ones(1,n_g)*alpha1(:,t,d)==-1];
%     %scenario 2
    Constraints=[Constraints,ones(1,n_g)*p2(:,t,d)+ones(1,n_w)*(m_c_n2{d,t})==ones(1,n_l)*0.9*4*(pd(1,:)')];
    Constraints=[Constraints,ones(1,n_g)*alpha2(:,t,d)==-1];
 end
end
%exisiting generating units
for d=1:n_d_c
for t=1:n_t
    for j=1:n_g
        Constraints=[Constraints,p1(j,t,d)<=pgmax(j)*x1(j,t,d)-sqrt((1-epsi)/epsi)*norm(alpha1(j,t,d)*ones(1,n_w)*cov_w_corr_c_n{d,t}^(1/2),2)];
        Constraints=[Constraints,p1(j,t,d)>=pgmin(j)*x1(j,t,d)+sqrt((1-epsi)/epsi)*norm(alpha1(j,t,d)*ones(1,n_w)*cov_w_corr_c_n{d,t}^(1/2),2)];   
        Constraints=[Constraints,p1(j,t,d)<=pgmax(j)*x1(j,t,d)-(pgmax(j)-QS(j))*u1(j,t,d)-sqrt((1-epsi)/epsi)*norm(alpha1(j,t,d)*ones(1,n_w)*cov_w_corr_c_n{d,t}^(1/2),2)];
        Constraints=[Constraints,p1(j,t,d)<=pgmax(j)*x1(j,t,d)];
        Constraints=[Constraints,p1(j,t,d)>=pgmin(j)*x1(j,t,d)];   
%         %scenario 2
        Constraints=[Constraints,p2(j,t,d)<=pgmax(j)*x2(j,t,d)-sqrt((1-epsi)/epsi)*norm(alpha2(j,t,d)*ones(1,n_w)*cov_w_corr_c_n2{d,t}^(1/2),2)];
        Constraints=[Constraints,p2(j,t,d)>=pgmin(j)*x2(j,t,d)+sqrt((1-epsi)/epsi)*norm(alpha2(j,t,d)*ones(1,n_w)*cov_w_corr_c_n2{d,t}^(1/2),2)];   
        Constraints=[Constraints,p2(j,t,d)<=pgmax(j)*x2(j,t,d)-(pgmax(j)-QS(j))*u2(j,t,d)-sqrt((1-epsi)/epsi)*norm(alpha2(j,t,d)*ones(1,n_w)*cov_w_corr_c_n2{d,t}^(1/2),2)];
        Constraints=[Constraints,p2(j,t,d)<=pgmax(j)*x2(j,t,d)];
        Constraints=[Constraints,p2(j,t,d)>=pgmin(j)*x2(j,t,d)];   
        
    end
end
end
 for d=1:n_d_c
  for t=1+1:n_t
     for j=1:n_g
     Constraints=[Constraints,p1(j,t,d)-p1(j,t-1,d)-x1(j,t-1,d)*c_r*ru(j)-(1-x1(j,t-1,d))*c_r*QS(j)<=-sqrt((1-epsi)/epsi)*norm([alpha1(j,t-1,d) alpha1(j,t-1,d) alpha1(j,t,d) alpha1(j,t,d)]*c_u*cov_corr_c_n{d,t}^(1/2),2)];
     Constraints=[Constraints,p1(j,t,d)-p1(j,t-1,d)+x1(j,t,d)*c_r*rd(j)+(1-x1(j,t,d))*c_r*QS(j)>=sqrt((1-epsi)/epsi)*norm([alpha1(j,t-1,d) alpha1(j,t-1,d) alpha1(j,t,d) alpha1(j,t,d)]*cov_corr_c_n{d,t}^(1/2),2)];
     Constraints=[Constraints,p1(j,t,d)-p1(j,t-1,d)-(pgmin(j)+c_r*ru(j))*x1(j,t,d)+pgmin(j)*x1(j,t-1,d)+(pgmin(j)+c_r*ru(j)-QS(j))*u1(j,t,d)<=-sqrt((1-epsi)/epsi)*norm([alpha1(j,t-1,d) alpha1(j,t-1,d) alpha1(j,t,d) alpha1(j,t,d)]*cov_corr_c_n{d,t}^(1/2),2)];
     Constraints=[Constraints,p1(j,t,d)-p1(j,t-1,d)+QS(j)*x1(j,t-1,d)-(pgmin(j)+c_r*rd(j)-QS(j))*u1(j,t,d)-(QS(j)-c_r*rd(j))*x1(j,t,d)>=sqrt((1-epsi)/epsi)*norm([alpha1(j,t-1,d) alpha1(j,t-1,d) alpha1(j,t,d) alpha1(j,t,d)]*cov_corr_c_n{d,t}^(1/2),2)];
     Constraints=[Constraints,p1(j,t-1,d)<=QS(j)*x1(j,t-1,d)+(pgmax(j)-QS(j))*(x1(j,t,d)-u1(j,t,d))-sqrt((1-epsi)/epsi)*norm(alpha1(j,t-1,d)*ones(1,n_w)*cov_w_corr_c_n{d,t-1}^(1/2),2)];
% Scenario 2

     Constraints=[Constraints,p2(j,t,d)-p2(j,t-1,d)-x2(j,t-1,d)*c_r*ru(j)-(1-x2(j,t-1,d))*c_r*QS(j)<=-sqrt((1-epsi)/epsi)*norm([alpha2(j,t-1,d) alpha2(j,t-1,d) alpha2(j,t,d) alpha2(j,t,d)]*cov_corr_c_n{d,t}^(1/2),2)];
     Constraints=[Constraints,p2(j,t,d)-p2(j,t-1,d)+x2(j,t,d)*c_r*rd(j)+(1-x2(j,t,d))*c_r*QS(j)>=sqrt((1-epsi)/epsi)*norm([alpha2(j,t-1,d) alpha2(j,t-1,d) alpha2(j,t,d) alpha2(j,t,d)]*cov_corr_c_n{d,t}^(1/2),2)];
     Constraints=[Constraints,p2(j,t,d)-p2(j,t-1,d)-(pgmin(j)+c_r*ru(j))*x2(j,t,d)+pgmin(j)*x2(j,t-1,d)+(pgmin(j)+c_r*ru(j)-QS(j))*u2(j,t,d)<=-sqrt((1-epsi)/epsi)*norm([alpha2(j,t-1,d) alpha2(j,t-1,d) alpha2(j,t,d) alpha2(j,t,d)]*cov_corr_c_n{d,t}^(1/2),2)];
     Constraints=[Constraints,p2(j,t,d)-p2(j,t-1,d)+QS(j)*x2(j,t-1,d)-(pgmin(j)+c_r*rd(j)-QS(j))*u2(j,t,d)-(QS(j)-c_r*rd(j))*x2(j,t,d)>=sqrt((1-epsi)/epsi)*norm([alpha2(j,t-1,d) alpha2(j,t-1,d) alpha2(j,t,d) alpha2(j,t,d)]*cov_corr_c_n{d,t}^(1/2),2)];
     Constraints=[Constraints,p2(j,t-1,d)<=QS(j)*x2(j,t-1,d)+(pgmax(j)-QS(j))*(x2(j,t,d)-u2(j,t,d))-sqrt((1-epsi)/epsi)*norm(alpha2(j,t-1,d)*ones(1,n_w)*cov_w_corr_c_n{d,t-1}^(1/2),2)];

     end
  end
end
%investment constraints
for d=1:n_d_c
for t=1:n_t
    for j=n_ge+1:n_g
        Constraints=[Constraints,x1(j,t,d)<=v(j-n_ge)];
% %         %scenario 2
        Constraints=[Constraints,x2(j,t,d)<=v(j-n_ge)];
    end
end
end

for d=1:n_d_c
for t=1+1:n_t
    for j=1:n_g
        Constraints=[Constraints,u1(j,t,d)>=x1(j,t,d)-x1(j,t-1,d)];
% %         %scneario 2
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

        

% %transmission constraints
for d=1:n_d_c
for t=1:n_t
    for j=1:n_line
        Constraints=[Constraints,Hg(j,1:n_g)*p1(:,t,d)+Hw(j,:)*m_c_n{d,t}-Hl(j,:)*(3*pd')<=10*fmax(j)-sqrt((1-epsi)/epsi)*norm(Hw(j,:)*cov_w_corr_c_n{d,t}^(1/2)+Hg(j,1:n_g)*alpha1(:,t,d)*ones(1,n_w)*cov_w_corr_c_n{d,t}^(1/2),2)];
        Constraints=[Constraints,Hg(j,1:n_g)*p1(:,t,d)+Hw(j,:)*m_c_n{d,t}-Hl(j,:)*(3*pd')>=-10*fmax(j)+sqrt((1-epsi)/epsi)*norm(Hw(j,:)*cov_w_corr_c_n{d,t}^(1/2)+Hg(j,1:n_g)*alpha1(:,t,d)*ones(1,n_w)*cov_w_corr_c_n{d,t}^(1/2),2)];
% Scenario 2
        Constraints=[Constraints,Hg(j,1:n_g)*p2(:,t,d)+Hw(j,:)*m_c_n2{d,t}-Hl(j,:)*(3*pd')<=10*fmax(j)-sqrt((1-epsi)/epsi)*norm(Hw(j,:)*cov_w_corr_c_n{d,t}^(1/2)+Hg(j,1:n_g)*alpha2(:,t,d)*ones(1,n_w)*cov_w_corr_c_n2{d,t}^(1/2),2)];
        Constraints=[Constraints,Hg(j,1:n_g)*p2(:,t,d)+Hw(j,:)*m_c_n2{d,t}-Hl(j,:)*(3*pd')>=-10*fmax(j)+sqrt((1-epsi)/epsi)*norm(Hw(j,:)*cov_w_corr_c_n{d,t}^(1/2)+Hg(j,1:n_g)*alpha2(:,t,d)*ones(1,n_w)*cov_w_corr_c_n2{d,t}^(1/2),2)];

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
    cost_base_operation=cost_base_operation+weigth_c(d,1)*(0.25*(cost_op(1:n_g)'*p1(1:n_g,t,d)+(cost_st(1:n_g))'*u1(1:n_g,t,d))+0.25*(cost_op(1:n_g)'*p2(1:n_g,t,d)+(cost_st(1:n_g))'*u2(1:n_g,t,d)));

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

result_ct{k}=Objective;
result_co{k}=cost_base_operation;
k=k+1;

end

time=toc

