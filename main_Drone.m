clear all
close all
clc
%% define dynamics
%define model constants
params.g = 9.81; % (m/s^2) gravity
l = 0.25;  % (m) half-width of quad rotor
m = 0.486;  % (kg) mass of the quad rotor. Used as the first nominal 
J = 0.00383; % (kgm^2), moment of inertia
lJ=l/J;m_inv=1/m;%make it linear in param
theta_tr=[lJ;m_inv];
theta_nom=[theta_tr(1)/2;theta_tr(2)*2];%Initial estimate
params.nx=6;
params.nu=2;
params.nw=1;
params.ny=2;params.C=[eye(params.ny),zeros(params.ny,params.nx-params.ny)];
params.ntheta=length(theta_tr);
params.discretization=0.025;

params.w_max=1e-1;%disturbance
params.v_max=1*[1e-3;1e-3;1e-2;1e-2;1e-2;1e-1]; %measurement noise
%noise: position very accurate, angle/velocity less so, angular velocity
T_sim=1e2/params.discretization;%simulation can be a lot shorter if only proposed is considered
%% noise
rng('default');%fix random seed for reproducability
noise_cl=diag(params.v_max)*(rand(params.nx,T_sim)-0.5);
w_cl=params.w_max*(rand(params.nw,T_sim)-0.5);%zeros(params.nw,T_sim);
%% MPC
params.N =5;% horizon length
params.M=10;
params.N_ext=params.N+params.M;
params.Q=eye(params.nx);params.R=1e-4*eye(params.nu);
params.T=4e3*eye(params.ny);
params.omega=10;
params.q_xi=5e3;
params.obstacles=true;%active/de-active obstacles
params.skip_noadapt=false;%turn of simulation of non-adaptive, which goes unstable (takes significantly longer...)
%input constraints
params.u_min = -1*ones(2,1);params.u_max=4*ones(2,1);
%state constraints
params.x_max=[10*ones(2,1);pi/3;1.5;1.5;pi];
params.x_min=-params.x_max;

%define obstacles  %x-y coordiates + radius (I always added one when the
%controller found a trivial path)
params.obs = [ ...
   -1.8  -1.1  0.20;  
   -1.3  -1.5  0.20;
   -1.3  -0.8  0.20;
   -1.5  -0.9  0.20;
   -0.8  -1.2  0.20;
   -0.4  -0.9  0.20;
   -0.4  -0.1  0.20;
   -0.8  -0.2  0.20;
   -1.7  -2.0  0.20;
   -2.0  -1.7  0.20]';  

% MPC horizon and dimensions
[solver,lb_con,ub_con]=MPC(params);
params_noTerm=params; params_noTerm.M=0;params_noTerm.N_ext=params.N;
[solver_noTerm,lb_con_noTerm,ub_con_noTerm]=MPC(params_noTerm);
%set LMS gain 
params.R_bounds=[inf(2,1);pi/2;params.x_max(3:5)*1.5];
Gamma=LMS_gain(params);
[K_nom,~]=feedback(theta_nom(:,1),params);%offline computed LQR K for rollout terminal penalty

%% closed loop
x_cl=zeros(params.nx,T_sim);
u_cl=zeros(params.nu,T_sim);
x_cl_noAdapt=zeros(params.nx,T_sim);
u_cl_noAdapt=zeros(params.nu,T_sim);
theta_hat=repmat(theta_nom,1,T_sim); 
theta_hat_noTerm=theta_hat;
theta_cl=repmat(theta_tr,1,T_sim);
comp=inf(T_sim,1);
x_0=[-2;-2;zeros(params.nx-2,1)];
x_cl(:,1)=x_0;x_cl_noAdapt(:,1)=x_0;x_cl_noTerm(:,1)=x_0;
y_init=zeros(params.nx*(2*params.N_ext+2*1+1)+params.nu*(params.N_ext+1),1);
y_init_noAdapt=y_init;y_init_noTerm=zeros(length(y_init)-params.M*(params.nx*2+params.nu),1);
for k=1:T_sim
fprintf('Simulation done %.2f\n', k/T_sim);

%initial guess
y_target(:,k)=zeros(2,1);
x_measure(:,k)=x_cl(:,k)+noise_cl(:,k);
if k>1 %LMS update
  theta_hat(:,k)=theta_hat(:,k-1)+Gamma*Phi_noisy'*(x_measure(:,k)-x_prediction);
  if max(eig(Phi_noisy*Gamma*Phi_noisy'))>1+1e-4
      error('LMS gain too large')
  end
end
[K_cur,~]=feedback(theta_hat(:,k),params);
%p_cur = [x_measure(:,k); theta_hat(:,k);y_target(:,k);vec(K_cur)];
p_cur = [x_measure(:,k); theta_hat(:,k);y_target(:,k);K_cur(:)];%vec no longer supported

t=tic;
sol = solver('x0',y_init , 'p', p_cur, 'lbg', lb_con, 'ubg', ub_con);
sol_X = reshape(full(sol.x(1:params.nx*(params.N_ext+1))), params.nx, params.N_ext+1);
sol_U = reshape(full(sol.x(params.nx*(params.N_ext+1)+1:params.nx*(params.N_ext+1)+params.N_ext*params.nu)), params.nu, params.N_ext);
sol_xs = full(sol.x(params.nx*(params.N_ext+1)+params.N_ext*params.nu+1:params.nx*(params.N_ext+1)+params.N_ext*params.nu+params.nx));
sol_us = full(sol.x(params.nx*(params.N_ext+1)+params.N_ext*params.nu+params.nx+1:params.nx*(params.N_ext+1)+params.N_ext*params.nu+params.nx+params.nu));
sol_Xi = reshape(full(sol.x(params.nx*(params.N_ext+1)+params.N_ext*params.nu+params.nx+params.nu+1:end)),params.nx,params.N_ext+1);
comp(k)=toc(t);%comptation time
y_error(k)=norm(params.C*x_cl(:,k)-y_target(:,k),2)^2;%tracking error
xi_cl(:,k)=sol_Xi(:,1);
viol(k)=norm(xi_cl(:,k),2)^2;%constraint violation
y_init=sol.x;%just take previous optimal solution
u_cl(:,k)=sol_U(:,1);    
x_cl(:,k+1)=dynamics(x_cl(:,k),u_cl(:,k),w_cl(:,k),params,theta_cl(:,k));
%LMs updat
[x_prediction,Phi_noisy] = dynamics(x_measure(:,k),u_cl(:,k),zeros(params.nw,1),params,theta_hat(:,k));

if params.skip_noadapt==false
%also implement without model adaptation 
x_measure_noAdapt(:,k)=x_cl_noAdapt(:,k)+noise_cl(:,k);
p_cur_noAdapt = [x_measure_noAdapt(:,k); theta_nom;y_target(:,k);K_nom(:)];
sol_noAdapt = solver('x0',y_init_noAdapt , 'p', p_cur_noAdapt, 'lbg', lb_con, 'ubg', ub_con);
sol_U_noAdapt = reshape(full(sol_noAdapt.x(params.nx*(params.N_ext+1)+1:params.nx*(params.N_ext+1)+params.N_ext*params.nu)), params.nu, params.N_ext);
sol_Xi_noAdapt = reshape(full(sol_noAdapt.x(params.nx*(params.N_ext+1)+params.N_ext*params.nu+params.nx+params.nu+1:end)),params.nx,params.N_ext+1);
y_error_noAdapt(k)=norm(params.C*x_cl_noAdapt(:,k)-y_target(:,k),2)^2;%tracking error
xi_cl_noAdapt(:,k)=sol_Xi_noAdapt(:,1);
viol_noAdapt(k)=norm(xi_cl_noAdapt(:,k),2)^2;%constraint violation
y_init_noAdapt=sol_noAdapt.x;%just take previous optimal solution
u_cl_noAdapt(:,k)=sol_U_noAdapt(:,1);    
x_cl_noAdapt(:,k+1)=dynamics(x_cl_noAdapt(:,k),u_cl_noAdapt(:,k),w_cl(:,k),params,theta_cl(:,k));
end

%also implement without terminal
x_measure_noTerm(:,k)=x_cl_noTerm(:,k)+noise_cl(:,k);
if k>1 %LMS update
  if max(eig(Phi_noisy_noTerm*Gamma*Phi_noisy_noTerm'))>1+1e-4
      error('LMS gain too large')
  end
  theta_hat_noTerm(:,k)=theta_hat_noTerm(:,k-1)+Gamma*Phi_noisy_noTerm'*(x_measure_noTerm(:,k)-x_prediction_noTerm);
end
[K_cur_noTerm,~]=feedback(theta_hat_noTerm(:,k),params);
p_cur_noTerm = [x_measure_noTerm(:,k); theta_hat_noTerm(:,k);y_target(:,k);K_cur_noTerm(:)];
sol_noTerm = solver_noTerm('x0',y_init_noTerm , 'p', p_cur_noTerm, 'lbg', lb_con_noTerm, 'ubg', ub_con_noTerm);
sol_X_noTerm = reshape(full(sol_noTerm.x(1:params.nx*(params_noTerm.N_ext+1))), params.nx, params_noTerm.N_ext+1);
sol_U_noTerm = reshape(full(sol_noTerm.x(params.nx*(params_noTerm.N_ext+1)+1:params.nx*(params_noTerm.N_ext+1)+params_noTerm.N_ext*params.nu)), params.nu, params_noTerm.N_ext);
sol_Xi_noTerm = reshape(full(sol_noTerm.x(params.nx*(params_noTerm.N_ext+1)+params_noTerm.N_ext*params.nu+params.nx+params.nu+1:end)),params.nx,params_noTerm.N_ext+1);
y_error_noTerm(k)=norm(params.C*x_cl_noTerm(:,k)-y_target(:,k),2)^2;%tracking error
xi_cl_noTerm(:,k)=sol_Xi_noTerm(:,1);
viol_noTerm(k)=norm(xi_cl_noTerm(:,k),2)^2;%constraint violation
y_init_noTerm=sol_noTerm.x;%just take previous optimal solution
u_cl_noTerm(:,k)=sol_U_noTerm(:,1);    
x_cl_noTerm(:,k+1)=dynamics(x_cl_noTerm(:,k),u_cl_noTerm(:,k),w_cl(:,k),params,theta_cl(:,k));
%LMs updat
[x_prediction_noTerm,Phi_noisy_noTerm] = dynamics(x_measure_noTerm(:,k),u_cl_noTerm(:,k),zeros(params.nw,1),params,theta_hat_noTerm(:,k));
end
%%
CompTime_av=mean(comp);
CompTime_st=std(comp);
fprintf('Computation time [ms] %.1f +- %.1f\n', CompTime_av*1e3,CompTime_st*1e3);
y_error_av=mean(y_error);
viol_av=mean(viol);

y_error_av_noTerm=mean(y_error_noTerm);
viol_av_noTerm=mean(viol_noTerm);


fprintf('Average tracking error (Proposed) [m] %.1f\n', y_error_av);
fprintf('Average tracking error (NoTerm) [m] %.1f\n', y_error_av_noTerm);

fprintf('Constraint violation (Proposed) %f\n', viol_av);
fprintf('Constraint violation (NoTerm) %f\n', viol_av_noTerm);

time_reach_proposed=round(find(y_error<1e-4,1)*params.discretization,1);
time_reach_noTerm=round(find(y_error_noTerm<1e-4,1)*params.discretization,1);
fprintf('Time to target (Proposed) [s] %f\n',time_reach_proposed);
fprintf('Time to target (NoTerm) [s] %f\n',time_reach_noTerm);

if params.skip_noadapt==false
y_error_av_noAdapt=mean(y_error_noAdapt);
viol_av_noAdapt=mean(viol_noAdapt);
    
fprintf('Average tracking error (NoAdapt) [m] %.1f\n', y_error_av_noAdapt);
fprintf('Constraint violation (NoAdapt) %f\n', viol_av_noAdapt);
end
%% plot
figure(1)
plot(x_cl(1,:),x_cl(2,:),'b-','linewidth',2)
hold on
if params.skip_noadapt==false
plot(x_cl_noAdapt(1,:),x_cl_noAdapt(2,:),'g--','linewidth',2)
end
plot(x_cl_noTerm(1,:),x_cl_noTerm(2,:),'r:','linewidth',2)
plot(0,0,'black+','linewidth',2)
for j=1:size(params.obs,2)
plot(params.obs(1,j),params.obs(2,j),'ko')
 ellipse(params.obs(3,j),params.obs(3,j),zeros(size(params.obs,1),1),params.obs(1,j),params.obs(2,j),'k',[],0.0);
%plot(params.obs(1,j)*[1,1],params.obs(2,j)+[0,params.obs(3,j)],'--')
end
plot(x_cl(1,:),x_cl(2,:),'b-','linewidth',2)
if params.skip_noadapt==false
plot(x_cl_noAdapt(1,:),x_cl_noAdapt(2,:),'g')
end
plot(x_cl_noTerm(1,:),x_cl_noTerm(2,:),'r:','linewidth',2)
plot(0,0,'black+','linewidth',2)
if params.skip_noadapt
legend({'Proposed','No-term','Goal','Obstacles'},'location','southeast','interpreter','latex')
else
legend({'Proposed','No-adapt','No-term','Goal','Obstacles'},'location','southeast','interpreter','latex')
end
axis([-2.05,0.1,-2.25,0.1])
xlabel('Position $p_1$ [m]','interpreter','latex')
ylabel('Position $p_2$ [m]','interpreter','latex')
set(gca, 'fontname','Arial','fontsize',16)
%print('Drone','-depsc')
print('Drone_compare','-depsc')
%%
figure
plot(u_cl')
hold on
plot(u_cl_noAdapt')
%% 
figure
plot(xi_cl','b')
hold on
plot(xi_cl_noTerm','r')
%%
figure
plot(1:T_sim,theta_hat','-b')
hold on
plot(1:T_sim,theta_hat_noTerm',':r')
plot([1,T_sim],[theta_tr,theta_tr]','--')
%legend('Prop','NoTerm','True')
%(theta_hat-theta_tr)')
%% tracking error
width=10;%for smoothing, nicer to understand plots
error_prop=sqrt(sum(x_cl(1:2,:).^2));
error_noTerm=sqrt(sum(x_cl_noTerm(1:2,:).^2));
error_prop_filt = filter(1/width*ones(width,1),1,error_prop);
error_noTerm_filt = filter(1/width*ones(width,1),1,error_noTerm);
figure
plot(error_prop_filt,'b')
hold on
plot(error_noTerm_filt,'r') 
%%
function [solver,lb_con,ub_con]=MPC(params)
import casadi.*
N=params.N;
M=params.M;
N_ext=params.N_ext;
R=params.R; Q=params.Q;
nx=params.nx; nu = params.nu; nw =params.nw;
%cost
% Define symbolic variables
x = SX.sym('x', params.nx);
u = SX.sym('u', params.nu);
% Decision variables
X = SX.sym('X', params.nx, N_ext+1);
U = SX.sym('U', params.nu, N_ext);
Xi = SX.sym('Xi', params.nx, N_ext+1);%slack variables
xs = SX.sym('xs',params.nx,1);
us = SX.sym('us',params.nu,1);
% Parameters (x0,theta,K)
X0 = SX.sym('X0', params.nx,1); 
theta = SX.sym('theta',params.ntheta);
y_target=SX.sym('y_target',params.ny);
K_loc = SX.sym('K_loc',params.nu,params.nx);
% Cost and constraints
cost = 0;
con=[];
con_eq=[]; 
con_eq=[con_eq; X(:,1)-X0];%initial constraint
con_eq=[con_eq;dynamics(xs,us,zeros(nw,1),params,theta)-xs];
for k = 1:N
    % dynamics 
    con_eq = [con_eq; X(:,k+1) - dynamics(X(:,k),U(:,k),zeros(params.nw,1),params,theta)];
    % stage cost
    cost = cost + (X(:,k)-xs)'*Q*(X(:,k)-xs) + (U(:,k)-us)'*R*(U(:,k)-us);
    %slack
    cost = cost + params.q_xi*Xi(:,k)'*Xi(:,k);
end
for k=N+1:N+M
    %extend horizon with terminal penalty
    % dynamics with feedback k
    con_eq = [con_eq;U(:,k)-us+K_loc*(X(:,k)-xs)];
    con_eq = [con_eq; X(:,k+1) - dynamics(X(:,k),U(:,k),zeros(params.nw,1),params,theta)];
    % stage cost
    cost = cost + params.omega*( (X(:,k)-xs)'*Q*(X(:,k)-xs) + (U(:,k)-us)'*R*(U(:,k)-us));
    cost = cost + params.omega*params.q_xi*Xi(:,k)'*Xi(:,k);
    
end
cost = cost + params.omega*( (X(:,N+M+1)-xs)'*Q*(X(:,N+M+1)-xs));
cost = cost + params.omega*params.q_xi*Xi(:,N+M+1)'*Xi(:,N+M+1);

%obstacle avoidance
con_obs=[];
for j=1:size(params.obs,2)
    con_obs=[con_obs;norm(params.C*xs-params.obs(1:2,j),2)-params.obs(3,j)];%steady-state
        for k=1:N+M+1
            con_obs=[con_obs; norm(params.C*(X(:,k)+Xi(:,k))-params.obs(1:2,j),2)-params.obs(3,j)];%should be >>=0 for avoidance
        end
end
con_obs_lb=zeros(size(con_obs));
con_obs_ub=inf(size(con_obs));

con_input=[vec(U);us];
con_input_ub=repmat(params.u_max,N_ext+1,1);
con_input_lb=repmat(params.u_min,N_ext+1,1);
%params.u_max*ones(length(con_input),1);
%con_input_lb=params.u_min*ones(length(con_input),1);
con_state=[vec(X+Xi);xs];%this is not the most computational efficent implementation, but the quickest
con_state_ub=repmat(params.x_max,N_ext+2,1);
con_state_lb=repmat(params.x_min,N_ext+2,1);
cost = cost + (params.C*xs-y_target)'*params.T*(params.C*xs-y_target);
con=[con_eq;con_input;con_state;con_obs];
lb_con = [zeros(size(con_eq));con_input_lb;con_state_lb;con_obs_lb]; % lower bound
ub_con = [zeros(size(con_eq));con_input_ub;con_state_ub;con_obs_ub];% upper bound
%
%add state/input constraints
% Build NLP
vars = [reshape(X, nx*(N_ext+1), 1); reshape(U, nu*N_ext, 1);xs;us;reshape(Xi,nx*(N_ext+1),1)];
nlp = struct('x', vars, 'f', cost, 'g', con, 'p', [X0; theta; y_target;vec(K_loc)]);
% Solver options
opts = struct;
opts.ipopt.print_level = 0;
opts.print_time = false;
solver = nlpsol('solver', 'ipopt', nlp, opts);
end
%%
function [x_new,Phi]=dynamics(x,u,w,params,theta)
%cont-time dynamics
%x=[p_1,p_2,pgi,v_1,v_2,dphi];
g=params.g;lJ=theta(1);m_inv=theta(2);
dot_x = [x(4)*cos(x(3)) - x(5)*sin(x(3));     %px
            x(4)*sin(x(3)) + x(5)*cos(x(3)); %pz
            x(6);                               %phi
            x(6).*x(5)-g*sin(x(3));             %vx
            -x(6)*x(4)-g*cos(x(3));            %vz
            zeros(1,size(x,2))]+...
            [zeros(4,2);m_inv, m_inv; lJ -lJ]*u+...
            [zeros(3,1);cos(x(3));-sin(x(3));0]*w; 
%discrete-time with Euler
%Phi=[zeros(5,2); u' -u'];%lJ;m
x_new=x+params.discretization*dot_x;
%x_new=f(x,u)+Phi(x,u)*theta
Phi=params.discretization*[zeros(4,2); 0, u(1)+u(2); u(1)-u(2), 0];
end
%%
function [K,u_s]=feedback(theta,params)
%compute local feedback around hovering
%1. compute hovering input 
%hovering: phi=0, v=0
g=params.g;lJ=theta(1);m_inv=theta(2);
u_s=g/(2*m_inv)*ones(2,1);
%Jacobian, cont.-time (hard coded)
A=[0,0,0,1,0,0;...
   0,0,0,0,1,0;...
   0,0,0,0,0,1;...
   0,0,-g,0,0,0;...
   0,0,0,0,0,0;...
   0,0,0,0,0,0];
B=[zeros(4,2);m_inv, m_inv; lJ -lJ];
A=eye(params.nx)+params.discretization*A;
B=params.discretization*B;
K=dlqr(A,B,params.Q,params.R);
K=-K;
end
%%
function Gamma=LMS_gain(params)
Gamma=sdpvar(params.ntheta,params.ntheta,'symmetric');
con_LMI=[Gamma>=0];
for phi=linspace(-params.R_bounds(3),-params.R_bounds(3),3)
for v1=linspace(-params.R_bounds(4),-params.R_bounds(4),3)
for v2=linspace(-params.R_bounds(5),-params.R_bounds(5),3)
for dphi=linspace(-params.R_bounds(6),-params.R_bounds(6),3)
for u1=linspace(params.u_min(1),params.u_max(1),3)
for u2=linspace(params.u_min(2),params.u_max(2),3)
   x=[0;0;phi;v1;v2;dphi];
   u=[u1;u2];
   [~,Phi]=dynamics(x,u,zeros(params.nw,1),params,zeros(params.ntheta,1));
   con_LMI=[con_LMI, Phi*Gamma*Phi' <=eye(params.nx)];
end
end
end
end
end
end
options = sdpsettings('verbose', 0);  % suppress solver output
optimize(con_LMI,-trace(Gamma),options);
Gamma=value(Gamma);
end