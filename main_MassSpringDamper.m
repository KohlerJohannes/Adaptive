%The code offline analysis regarding stability with a given horizon is
%adapted from the code in Köhler, Johannes, Melanie N. Zeilinger, and Lars Grüne. 
%"Stability and performance analysis of NMPC: Detectable stage costs and general terminal costs." 
% IEEE Transactions on Automatic Control 68.10 (2023): 6114-6129.
%%
clear all
close all
clc
%% Model
%generate (A,B,C) for different model parameters, and make sure nominal MPC
%is stabilizing for all models (offline)
number_systems=10;
n_i=2;
n_x=n_i*number_systems;
n_u=1;n_y=1;
C=[1,0,zeros(1,n_i*(number_systems-1))]; 
perc=0.5;%parmaetric uncertainty
w_max=1e-2;v_max=1e-2;%magnitude noise
%nominal model parameters
mass_nom=1; 
spring_nom=40;damp_nom=13;%2*sqrt(mass_nom*spring_nom);
discretization=0.5;
%create mass/min parameters with given uncertainty
mass_min=mass_nom*(1-perc); mass_max=mass_nom*(1+perc);
spring_min=spring_nom*(1-perc); spring_max=spring_nom*(1+perc);
damping_min=damp_nom*(1-perc); damping_max=damp_nom*(1+perc);
R=1e-3;Q=eye(n_i*number_systems);%state/input penalty
%state constraints
Ax=[C];bx=0.7*ones(1,1);
%input constraint
u_max=25;
Au=[1;-1];bu=u_max*ones(2,1);
soft_weight=2;%soft penalty constraints
Q_soft=Q+soft_weight*Ax'*Ax;  %needed for analysis soft-constraint
grid=10;%check stability based on a grid of parameter values
M_des=11/discretization;%horizon to ensure stability (error will be thrown if too short)
omega=5;%set penalty terminal
epsilon_f=0;%initialze epsilon_f (computed below)
% to verify condition for all dynamics, we grid
for mass=linspace(mass_min,mass_max,grid)
for spring=linspace(spring_min,spring_max,grid)
for damp=linspace(damping_min,damping_max,grid)
    [A,B]=systemmatrices(mass,spring,damp,number_systems,discretization);    
     %This computes epsilon_f using a generalized eigenvalue 
     temp=0;
     for k=0:M_des-1
     temp=temp+A'^k*Q*A^k;
     end
     epsilon_temp=max(eig(A'^(M_des)*Q_soft*A^(M_des)+(1-omega)/omega*Q,temp));
     epsilon_f=max(epsilon_temp,epsilon_f);
end
end
end
if epsilon_f>0
    error('stability not guaranteed for all N')
else
    disp('horizon $M$ suitable for stability, even with N=1')
end
M=M_des;
% show that no common Lyap exists (by changing perc, one can see when robust control methods beging to fail)
P_common=sdpvar(n_x,n_x,'symmetric');
con_LMI=[P_common>=eye(n_x)];%ensures positive definite, can be scaled arbitraril
for mass=linspace(mass_min,mass_max,2)
for spring=linspace(spring_min,spring_max,2)
for damp=linspace(damping_min,damping_max,2)
    [A,B]=systemmatrices(mass,spring,damp,number_systems,discretization);
    %check if a common Lyapunov function exists
     con_LMI=[con_LMI,A'*P_common*A-P_common<=0];
end 
end 
end 
ops = sdpsettings('solver', 'mosek','verbose',0); 
res_LMI=optimize(con_LMI,0,ops);
if res_LMI.problem==1
    disp('no common Lyapunov function exists')
elseif res_LMI.problem==0
    error('common Lyapunov function exists')
else
    error('unexpected error')
end
%% LMS
%compute gain Gamma
%Phi= [x;u]\krocker I_nx
%Gamma=gamma\kronecker I_nx;
%--> Required condition:
% [x;u]*gamma*[x;u]'<=I
x_max=10*ones(n_x,1);%set conservative bound on maximal state we expect
gamma=diag([1./(x_max.^2);1./(u_max.^2)])/(n_x+n_u);
Gamma=kron(gamma,eye(n_x));
Phi=kron([x_max;u_max]',eye(n_x));
if eig(Phi*Gamma*Phi')>1
    error('gain wrong')%sanity check if gain was computed correctly
end
%% formulate proposed MPC
% cost
N=3/discretization; %division by discretization means horizon of 3 [s]
ell = @(x,u) x' * Q * x + u' * R * u;
T=1e3;
Vo = @(xr,yr) (C*xr-yr)' * T * (C*xr-yr);
Vo_NoRollout = @(xr,yr) (C*xr-yr)' * T * (C*xr-yr)*(N+0)/(N+M);%make it scaled to have similar weight with shorter M horizon in no rollout
ell_xi = @(x,u,xi) ell(x,u)+soft_weight*(xi)'*(xi);
ell_xi0= @(x,u,xi) ell(x,u);%no soft penalty 

yalmip_proposed = create_mpc(N, M,n_x,n_u, Ax, bx, Au, bu,ell_xi,Vo,omega);
%no terminal
yalmip_no_rollout = create_mpc(N,0,n_x,n_u, Ax, bx, Au, bu,ell_xi,Vo_NoRollout,omega);
%no state -> remove soft penalty using ell_xi0 instead of ell_xi (still has state constraint on artifical;
%could also change bx otherwise)
yalmip_no_state_con = create_mpc(N, M,n_x,n_u, Ax, bx, Au, bu,ell_xi0,Vo,omega);
%no adaptation (doesn't need separate yalimp, just set A const)
%% closed-loop simulation
Tsim=300/discretization;
x0=0*ones(n_x,1);%rand(n_x,1);
n_steps=20/discretization;
yT=[0.5*ones(n_steps,1);-0.5*ones(n_steps,1);1*ones(n_steps,1);-1*ones(n_steps,1);1*ones(n_steps,1);-1*ones(n_steps*2,1);bx*ones(Tsim-7*n_steps,1)];
yT_feasible=min(yT,bx);%best feasible output reference y_rd
[A_tr,B_tr]=systemmatrices(mass_min,spring_min,damping_min,number_systems,discretization);
[A_0,B_0]=systemmatrices(mass_nom,spring_nom,damp_nom,number_systems,discretization);
theta_tr=[A_tr,B_tr]; theta_tr=theta_tr(:);
theta_0=[A_0,B_0];theta_0=theta_0(:);
n_theta=length(theta_0);
w=w_max*(rand(n_x,Tsim)-0.5);
v=v_max*(rand(n_x,Tsim)-0.5);
%check if true system is minimum phase (since short horizon controllers can
%fail, even with a perfect model in this case)
sys=ss(A_tr,B_tr,C,0,discretization);
zeros_TF=tzero(tf(sys));
if max(abs(zeros_TF))>0
    disp('unstable zero dynamics')
    max(abs(zeros_TF))
end
%%
% simulation of proposed MPC
%initialize
x_prop=zeros(n_x,Tsim);x_NoRollout=zeros(n_x,Tsim);x_NoState=zeros(n_x,Tsim);x_NoAdapt=zeros(n_x,Tsim);
x_prop(:,1)=x0;
x_NoRollout(:,1)=x0;
x_NoState(:,1)=x0;
x_NoAdapt(:,1)=x0;
x_noisy_prop(:,1)=x0+v(:,1);
x_noisy_NoRollout(:,1)=x0+v(:,1);
x_noisy_NoState(:,1)=x0+v(:,1);
x_noisy_NoAdapt(:,1)=x0+v(:,1);
x_hat_1_prop=zeros(n_x,Tsim);x_hat_1_NoRollout=zeros(n_x,Tsim);x_hat_1_NoState=zeros(n_x,Tsim);x_hat_1_NoAdapt=zeros(n_x,Tsim);
u_prop=zeros(n_u,Tsim);u_NoRollout=zeros(n_u,Tsim);  u_NoState=zeros(n_u,Tsim); u_NoAdapt=zeros(n_u,Tsim);

theta_hat_prop=zeros(n_theta,Tsim);theta_hat_NoRollout=zeros(n_theta,Tsim);theta_hat_NoState=zeros(n_theta,Tsim);
theta_hat_prop(:,1)=theta_0;
theta_hat_NoRollout(:,1)=theta_0;
theta_hat_NoState(:,1)=theta_0;
solvetime=inf(Tsim-1,1);
for k=1:Tsim-1
    AB_prop=reshape(theta_hat_prop(:,k),n_x,n_x+n_u);
    A_prop=AB_prop(1:n_x,1:n_x);B_prop=AB_prop(1:n_x,n_x+1:n_x+n_u);
    AB_NoRollout=reshape(theta_hat_NoRollout(:,k),n_x,n_x+n_u);
    A_NoRollout=AB_NoRollout(1:n_x,1:n_x);B_NoRollout=AB_NoRollout(1:n_x,n_x+1:n_x+n_u);
    AB_NoState=reshape(theta_hat_NoState(:,k),n_x,n_x+n_u);
    A_NoState=AB_NoState(1:n_x,1:n_x);B_NoState=AB_NoState(1:n_x,n_x+1:n_x+n_u);
    t=tic;
    %proposed MPC 
    [res_prop,flag] = yalmip_proposed({x_noisy_prop(:,k),A_prop ,B_prop,yT(k)});
    assert(flag==0||flag==3||flag==5);%0:sucess, 3:max-iteration, 5=lack of progress, 
    solvetime(k)=toc(t);
    %no rollout
    [res_NoRollout,flag] = yalmip_no_rollout({x_noisy_NoRollout(:,k),A_NoRollout,B_NoRollout,yT(k)});
    assert(flag==0||flag==3||flag==5);%0:sucess, 3:max-iteration, 5=lack of progress,
    %no state constraints
    [res_NoState,flag] = yalmip_no_state_con({x_noisy_NoState(:,k),A_NoState,B_NoState,yT(k)});
    assert(flag==0||flag==3||flag==5);%0:sucess, 3:max-iteration, 5=lack of progress,
    %no model adptation
    [res_NoAdapt,flag] = yalmip_proposed({x_noisy_NoAdapt(:,k),A_0 ,B_0,yT(k)});%don't update parameter estimate
    assert(flag==0||flag==3||flag==5);%0:sucess, 3:max-iteration, 5=lack of progress, 
    %retrive input + one-step prediction
    U_prop=res_prop{2};X_prop=res_prop{3};Xr_prop=res_prop{4};
    u_prop(:,k)=U_prop(:,1);x_hat_1_prop(:,k)=X_prop(:,2);%one-step prediction
    U_NoRollout=res_NoRollout{2};X_NoRollout=res_NoRollout{3};Xr_NoRollout=res_NoRollout{4};
    u_NoRollout(:,k)=U_NoRollout(:,1);x_hat_1_NoRollout(:,k)=X_NoRollout(:,2);%one-step prediction
    U_NoState=res_NoState{2};X_NoState=res_NoState{3};Xr_NoState=res_NoState{4};
    u_NoState(:,k)=U_NoState(:,1);x_hat_1_NoState(:,k)=X_NoState(:,2);%one-step prediction
    U_NoAdapt=res_NoAdapt{2};X_NoAdapt=res_NoAdapt{3};Xr_NoAdapt=res_NoAdapt{4};
    u_NoAdapt(:,k)=U_NoAdapt(:,1);x_hat_1_NoAdapt(:,k)=X_NoAdapt(:,2);%one-step prediction
    %1. closed-loop without adaptaiton, noise, etc.
    x_prop(:,k+1)=A_tr*x_prop(:,k)+B_tr*u_prop(:,k)+w(:,k);
    x_NoRollout(:,k+1)=A_tr*x_NoRollout(:,k)+B_tr*u_NoRollout(:,k)+w(:,k);
    x_NoState(:,k+1)=A_tr*x_NoState(:,k)+B_tr*u_NoState(:,k)+w(:,k);
    x_NoAdapt(:,k+1)=A_tr*x_NoAdapt(:,k)+B_tr*u_NoAdapt(:,k)+w(:,k);
    %noisy measurement
    x_noisy_prop(:,k+1)=x_prop(:,k+1)+v(:,k+1);
    x_noisy_NoRollout(:,k+1)=x_NoRollout(:,k+1)+v(:,k+1);
    x_noisy_NoState(:,k+1)=x_NoState(:,k+1)+v(:,k+1);
    x_noisy_NoAdapt(:,k+1)=x_NoAdapt(:,k+1)+v(:,k+1);
    %parameter update
    Phi_prop=kron([x_noisy_prop(:,k);u_prop(:,k)]',eye(n_x));
    theta_hat_prop(:,k+1)=theta_hat_prop(:,k)+Gamma*Phi_prop'*(x_noisy_prop(:,k+1)-x_hat_1_prop(:,k));
    Phi_NoRollout=kron([x_noisy_NoRollout(:,k);u_NoRollout(:,k)]',eye(n_x));
    theta_hat_NoRollout(:,k+1)=theta_hat_NoRollout(:,k)+Gamma*Phi_NoRollout'*(x_noisy_NoRollout(:,k+1)-x_hat_1_NoRollout(:,k));
    Phi_NoState=kron([x_noisy_NoState(:,k);u_NoState(:,k)]',eye(n_x));
    theta_hat_NoState(:,k+1)=theta_hat_NoState(:,k)+Gamma*Phi_NoState'*(x_noisy_NoState(:,k+1)-x_hat_1_NoState(:,k));
    if max(eig(Phi_prop*Gamma*Phi_prop'))>1 || max(eig(Phi_NoRollout*Gamma*Phi_NoRollout'))>1|| max(eig(Phi_NoState*Gamma*Phi_NoState'))>1
        error('parameter gain wrong')
    end
    %parameter error
    V_theta_prop(k)=(theta_hat_prop(:,k)-theta_tr)'*inv(Gamma)*(theta_hat_prop(:,k)-theta_tr);
    V_theta_NoRollout(k)=(theta_hat_NoRollout(:,k)-theta_tr)'*inv(Gamma)*(theta_hat_NoRollout(:,k)-theta_tr);
    V_theta_NoState(k)=(theta_hat_NoState(:,k)-theta_tr)'*inv(Gamma)*(theta_hat_NoState(:,k)-theta_tr);
    %constraint vioation
    con_prop(k)=norm(max(Ax*x_prop(:,k)-bx,0),2)^2;
    con_NoRollout(k)=norm(max(Ax*x_NoRollout(:,k)-bx,0),2)^2;
    con_NoState(k)=norm(max(Ax*x_NoState(:,k)-bx,0),2)^2;
    con_NoAdapt(k)=norm(max(Ax*x_NoAdapt(:,k)-bx,0),2)^2;
    %tracking error
    track_prop(k)=norm(yT(k)-C*x_prop(:,k),2)^2;
    track_NoRollout(k)=norm(yT(k)-C*x_NoRollout(:,k),2)^2;
    track_NoState(k)=norm(yT(k)-C*x_NoState(:,k),2)^2;
    track_NoAdapt(k)=norm(yT(k)-C*x_NoAdapt(:,k),2)^2;
    %tracking error to feasible reference
    track_feas_prop(k)=norm(yT_feasible(k)-C*x_prop(:,k),2)^2;
    track_feas_NoRollout(k)=norm(yT_feasible(k)-C*x_NoRollout(:,k),2)^2;
    track_feas_NoState(k)=norm(yT_feasible(k)-C*x_NoState(:,k),2)^2;
    track_feas_NoAdapt(k)=norm(yT_feasible(k)-C*x_NoAdapt(:,k),2)^2; 
 end
%% save metrics
K=200/discretization %onyl consider first 200 seconds
disp('output error - optimal reference')%sum(track_feas_NoState),
round([sum(track_feas_prop(1:K)),sum(track_feas_NoRollout(1:K)),sum(track_feas_NoAdapt(1:K))]./sum(track_feas_prop(1:K)),2)
disp('constraints')%sum(con_NoState)
round([sum(con_prop(1:K)),sum(con_NoRollout(1:K)),sum(con_NoAdapt(1:K))]/sum(con_prop(1:K)),2)
%disp('output error - target')%,sum(track_NoState)
%[sum(track_prop),sum(track_NoRollout),sum(track_NoAdapt)]
disp('decrease $V_theta$ [%]')
100*round(1-sqrt(V_theta_prop(K))/sqrt(V_theta_prop(1)),2)
disp('solvetimes [ms]')
round(mean(solvetime(1:K))*1e3,0)
round(std(solvetime(1:K))*1e3,0)
%% Plots
close all
t=(0:Tsim-1)*discretization;
plot(t,C*x_prop(:,:),'b','linewidth',1)
hold on
plot(t,C*x_NoRollout(:,:),'green','linewidth',1)
%plot(C*x_NoState(:,:),'y')
plot(t,C*x_NoAdapt(:,:),'red','linewidth',1)
plot(t,yT,'black--','linewidth',1)
plot(t,bx*ones(Tsim,1),'black:','linewidth',1.5)
plot(t,C*x_prop(:,:),'b','linewidth',1)
axis([0,230,-1.1,1.1])
legend({'Proposed', 'No adaptation','Target $y^d$', 'Constraints ${X}$'},'location','southeast','interpreter','latex')
legend({'Proposed','No-term', 'No-adapt','Target $y^d$', 'Constraint'},'location','southeast','interpreter','latex')
xlabel('Time [s]','interpreter','latex')
ylabel('Position $y$ [m]','interpreter','latex')
set(gca, 'fontname','Arial','fontsize',16)
print('MassSpring_Adaptive','-depsc')
%%
figure
plot(u_prop)
hold on
plot(u_NoRollout)
plot(u_NoState)
plot(u_NoAdapt)
legend('Prop','NoRollout','NoState','NoAdapt')
%clear temp sys A b C ans i k Phi
%save model %A,B,C, 
%%
%%
function [A,B]=systemmatrices(m,k,d,M,h)
n_i=2;
if M<3
   error('not programmed') 
end
n_x=n_i*M;
A_c=zeros(n_x);
%i=1:also contected to ground
A_c(1,:)=[0,1,zeros(1,n_x-n_i)];
A_c(2,:)=1/m*[-k-k,-d-d,k,d,zeros(1,n_x-2*n_i)];
%1<i<M
for i=2:M-1
A_c((i-1)*n_i+1,:)=[zeros(1,(i-1)*n_i+1),1,zeros(1,n_x-i*n_i)];
A_c(i*n_i,:)=1/m*[zeros(1,n_i*(i-2)),k,d,-2*k,-2*d,k,d,zeros(1,n_x-n_i*(i+1))];
end
%i=M, measured output
A_c(M*n_i-1,:)=[zeros(1,M*n_i-1),1];
A_c(M*n_i,:)=1/m*[zeros(1,n_x-2*n_i),k,d,-k,-d]; 
%acutation on last mass
n_u=1;
B_c=1/m*[zeros(n_x-n_i,1);0;1;];
[A,B]=c2d(A_c,B_c,h);
%output: position last mass
%A=A_c*h+eye(n_x);
%B=B_c*h;
end
function yalmip_optimizer = create_mpc(N, M,n,m, Ax, bx, Au, bu,ell_xi,Vo,omega)
%nominal MPC
 %to code:
  %terminal horizon
    %nx = size(A,1);
    %nu = size(B,2);
    Xr=sdpvar(n,1,'full');
    Ur=sdpvar(m,1,'full');
    yT=sdpvar(m,1,'full');

    % define optimization variables
    U = sdpvar(m,N+M-1,'full');
    X = sdpvar(n,N+M,'full');
    Xi=sdpvar(length(bx),N+M,'full');
    x0 = sdpvar(n,1,'full');
    A = sdpvar(n,n,'full');
    B = sdpvar(n,m,'full');
    C = sdpvar(m,n,'full');
    objective = 0;
    constraints = [X(:,1) == x0]; % initial condition
    for k = 1:N+M-1
        constraints = [constraints, X(:,k+1) == A * X(:,k) + B * U(:,k)]; % dynamics constraint
        constraints = [constraints, Au * U(:,k) <= bu]; % hard input constraint
        constraints = [constraints, Ax * X(:,k) <= bx+Xi(:,k)]; % hard state constraint
        constraints = [constraints, Xi(k)>=0];
        if k>=N
        objective = objective + omega*ell_xi(X(:,k)-Xr,U(:,k)-Ur,Xi(:,k));
        else
        objective = objective + ell_xi(X(:,k)-Xr,U(:,k)-Ur,Xi(:,k));
        end
    end
    for k= N:N+M-1
        constraints = [constraints, U(:,k) == Ur]; %
    end
    objective=objective + Vo(Xr,yT);
    constraints = [constraints, A*Xr+B*Ur==Xr];
    constraints = [constraints, Au * Ur <= bu];
    constraints = [constraints, Ax * Xr <= bx];
    % terminal objective
    %objective = objective + ell_f_LQR(X(:,N)); 
    %constraints = [constraints, G * X(:,N) <= g];
    % setup optimizer object
    ops = sdpsettings('verbose',0,'solver','quadprog');%'mosek' %'quadprog'  %,'warmstart',1
    parameters = {x0, A, B,yT};   
    yalmip_optimizer = optimizer(constraints,objective,ops,parameters,{objective, U, X,Xr,Ur});
end