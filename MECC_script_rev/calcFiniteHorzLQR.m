function [U,X,K,P] = calcFiniteHorzLQR(A,B,x0,Q,R,P_N,N,umax,umin)
% Gets the saturated DLQR sequence over the horizon
m = size(R,1);
n = size(Q,2);

% Set bounds to infinity if none given (i.e. regular LQR)
if nargin < 8
    umin = -Inf*ones(m,1);
    umax = Inf*ones(m,1);
end

% Initialize
U = zeros(m*N,1);
X = zeros(n,N+1);
P = zeros(n,n,N+1);
K = zeros(m,n,N);

% Run backward recursion to get P, K
P(:,:,end) = P_N;
PPlus = P_N;
for k = 1:N
    i = N+1 - k;
    K_i =  (R+B'*PPlus*B)\B'*PPlus*A;
    P_i = Q + A'*PPlus*A - A'*PPlus*B*K_i;
    PPlus = P_i;
    P(:,:,i) = P_i;
    K(:,:,i) = K_i;
end


% Run forward recursion to get U, X
x_i = x0;
X(:,1) = x0;
for i = 1:N
   u_i = -K(:,:,i)*x_i; % Get LQR
   
   % Saturate
   for j = 1:m
       if u_i(j) > umax(j)
           u_i(j) = umax(j);
       elseif u_i(j) < umin(j)
           u_i(j) = umin(j);
       end
   end
   
   % Store and increment x_i
   U(1+m*(i-1):m+m*(i-1)) = u_i;
   x_i = A*x_i + B*u_i;
   X(:,i+1) = x_i;
end
end

