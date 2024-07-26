%Matlab simulation for paper submitted to MECC2024
%Author: Michele Ambrosino 

clc, close all,clear all

%folder with some utilities functions
addpath('extraFunctionV2_new') %
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))

%% * Set simulation parameters

% * Time vector
dt = 0.1; %samplig time
N = 150; %predition horizon
TEnd = N*dt; %end time simulaiton
t = linspace(0,TEnd,N+1); %time vector

% * Initial Conditions of the system
x0(:,1) = [0.5;0.5]; 
x0(:,2) = [0.4;0.7]; 
x0(:,3) = [-0.8;-0.5]; 
x0(:,4) = [-1;-0.5];
x0(:,5) = [0.8;0.65];
x0(:,6) = [-1.1;-1];
x0(:,7) = [1;-1];
x0(:,8) = [0.6;-0.5];
x0(:,9) = [0.2;0.8];
x0(:,10) = [0.7;0.7];


for i_x0 = 1:size(x0,2)
    i_x0
    X0 = x0(:,i_x0);
% * State space matrices   
global A B Ad Bd

% State Vector (X): omega1
%                   omega2


% Continuous-time dynamics
n_spin = 1; % rad/s, base spin rate
J_ratio = 0.05; % J3/J 
lambda_const = n_spin*(1 - J_ratio); % rad/s^2, lambda = n*(J-J3)/J
Jt = 1; % Nm
J3 = J_ratio*Jt;
A = [0,lambda_const;
     -lambda_const, 0];
B = [0; 1/Jt];

% Get discrete state-space matrices
[Ad,Bd] = c2d(A,B,dt);
n = size(A,1);
m = size(B,2);


%% * Set control parameters and run the simulation

% Set horizon (same as the intial setting)
T_MPC = dt;
N_MPC = N; 

% * Set P, Q, R
Q = diag([1 1]); %  State error cost
R = 20; % Control effort cost
[K,P] = dlqr(Ad,Bd,Q,R);

% * Set constraints 
umax = 0.2*ones(m,1); %input constraints
umin = -umax;
xmax = ones(n,1)*1000; % just a placeholder, doesn't get used since we use PGM
xmin = -xmax;
MaxIter = 10000; 
xTol_MPC = 1e-8;  % ||x||_inf optimization tolerance for MPC optimizer


%% Calculate alpha, d, and lambda

% -------------- Calculate constant d and alpha --------------
% Compute alpha according to Marco's formula - taking the minimum
% OmegaBar = { x | -Kx in U }
ACon = [-K; K];
bCon = [umax; -umax];
alpha = 1e12; % initialize
invP = inv(P);
for i = 1:length(bCon)
    b_i = bCon(i);
    a_i = ACon(i,:)';
    VNew = b_i^2/(a_i'*invP*a_i); % assumes x_eq = 0
    if VNew < alpha
        alpha = VNew;
    end
end
sqrtQ  = sqrtm(Q);
invsqrtQ = inv(sqrtQ);
eigVec = eig(invsqrtQ*P*invsqrtQ);
gamma = alpha/max(eigVec);

InitialNorm = X0'*P*X0;
alpha;



%% Calculate lambda and beta

% Get initial cost given lambda
lambda = 10;
[H,G,W,ACon,FCon,LCon,S,M] = generateQPMatrices_compressed(N_MPC,Ad,Bd,lambda*P,Q,R,xmax,xmin,umax,umin);
H = 0.5*(H+H');
% 
% Get n*m length vectors u_low and u_up
u_low = zeros(m*(N-1),1);
u_up = zeros(m*(N-1),1);
for j = 1:N
    u_low(1 + m*(j-1) : m + m*(j-1)) = umin;
    u_up(1 + m*(j-1) : m + m*(j-1)) = umax;
end

% Calculate true value function 
H_QP = H;
f_QP = G*X0;
[U_MPC,iterCount] = accelProjGradSolver(H_QP,f_QP,zeros(m*N,1),u_low,u_up,700,1e-12,0);
 
% Calculate cost
COST = U_MPC'*H*U_MPC + 2*U_MPC'*G*X0 + X0'*W*X0;
BOUND = N*gamma + lambda*alpha;
beta = BOUND - COST;


% Get initial controls and states
u_plot = reorderMPCVector(U_MPC,m,N);
XVec = S*X0 + M*U_MPC;
x_plot = reorderMPCVector(XVec,n,N);

if InitialNorm < alpha
    error('Initial state is in terminal region')
end
if beta < 0
    error('Negative beta')
end

%% Get shrinking horizon constants

% Define the sequence beta0,...,betaN. Implicitly just know betaN = 0, no
% need to actually store this.
betaVec = zeros(N,1);
betaVec(1) = beta;

% LINEAR BETA SEQUENCE
for jj = 2:N
    kk = jj - 1;
    betaVec(jj) = min([beta*(N - kk)/N, 0.99*(N-kk)*gamma + lambda*alpha]);
end

% * Add controller arguments to a structure
% controlArgs.P = Sf;
controlArgs.P = P;
controlArgs.Q = Q;
controlArgs.R = R; 
controlArgs.K = K;
controlArgs.N = N_MPC; % Horzion lengt
controlArgs.T = T_MPC; % Sampling period
controlArgs.MaxIter = MaxIter;
controlArgs.xTol = xTol_MPC;
controlArgs.umax = umax; % Control constraints
controlArgs.umin = umin;
controlArgs.alpha = alpha;
controlArgs.gamma = gamma;
controlArgs.lambda = lambda;
controlArgs.dMax = 0;
controlArgs.betaVec = betaVec;
controlArgs.X0 = X0;

% Run function for the first sequence
output1 = calcParamsSHMPC(controlArgs);
init1 = output1.ellBarVec(1);


%% Run the analysis with disturbances

% Run offline analysis 
dBar  = (min(output1.wBar));
d1Bar = i_x0/10*(min(output1.wBar));
controlArgs.dMax = d1Bar;
output1_dist = calcParamsSHMPC(controlArgs);

% Online disturbance
dOnline = d1Bar; % 0.008; % 4*(min(output1.wBar))
controlArgs.dMax = dOnline;

controlArgs.V0 = COST;
controlArgs.U_Init = U_MPC;

% Run the simulation open loop
[X_ol,U_ol] = runShrinkingHorizonMPC_openLoop(X0,U_MPC,controlArgs);
COST_opt = U_ol*H*U_ol' + 2*U_ol*G*X0 + X0'*W*X0;

iter = output1_dist.ellBarVec; 

%loop of simulation for the three cases

output1_dist.ellBarVec = iter;
[X_fixed,U_fixed,iter_out] = runShrinkingHorizonMPC_fixedIter(X0,controlArgs,output1_dist);
COST_sub = U_fixed*H*U_fixed' + 2*U_fixed*G*X0 + X0'*W*X0;
%sum(iter_out)
%eval gap
[Gap,R_U,R_X,delta,L,psi,psi_bar,bound] = eval_R_v2(Ad,Bd,controlArgs,output1_dist,dBar,N);

res.X = X_fixed;
res.U = U_fixed;
res.cost = COST_sub;
res.nIter = iter;
res.Gap = Gap; 

res_fin.s(i_x0) = res;



end


close all
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))

set(0,'defaultLineLineWidth',2)
set(0,'defaultAxesFontSize',25)
labelsize = 35;
legendsize = 22;
colorMatrix = [0 0 0;  135 135 135; 60 60 60; 55 126 184; 228 26 28; 77 175 74]./255;



%% Plotting


figure('units','normalized','outerposition',[0 0 0.6 0.8])
%subplot(2,1,1);
% Plot the sublevel set
NGrid = 250;
x1 = linspace(-1,1,NGrid);
x2 = linspace(-1,1,NGrid);
[XX,YY]=meshgrid(x1,x2);
Z= P(1,1)*XX.^2 + P(1,2).*XX.*YY + P(2,1).*XX.*YY+ P(2,2)*YY.^2;
[~,hContour] = contourf(XX,YY,-Z,-[alpha alpha]);
setColor = [255 153 1]/255;
set(hContour,'facecolor',setColor,'edgecolor',setColor,'linewidth',1)
drawnow;  % this is important, to ensure that FacePrims is ready in the next line!
hFills = hContour.FacePrims;  % array of TriangleStrip objects
[hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
for idx = 1 : numel(hFills)
   hFills(idx).ColorData(4) = 0.4*255;   % default=255
end


hold on; box on; grid on
xlabel('$\dot{\theta}_1 \; [{{rad}\over{s}}]$','interpreter','latex','fontsize',labelsize)
ylabel('$\dot{\theta}_2 \; [{{rad}\over{s}}]$','interpreter','latex','fontsize',labelsize)

for i = 1:size(x0,2)
plot(res_fin.s(i).X(1,:),res_fin.s(i).X(2,:),'k',res_fin.s(i).X(1,1),res_fin.s(i).X(2,1),'or',res_fin.s(i).X(1,end),res_fin.s(i).X(2,end),'ob')
hold on
end

