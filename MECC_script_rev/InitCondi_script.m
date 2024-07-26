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
X0 = [0.2;0.2]; 

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

BOUND = N*gamma + lambda*alpha;

if InitialNorm < alpha
    error('Initial state is in terminal region')
end
