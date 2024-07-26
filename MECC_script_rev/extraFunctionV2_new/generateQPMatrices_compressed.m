function [H,G,W,ACon,FCon,LCon,AHat,BHat] = generateQPMatrices_compressed(N,Ad,Bd,P,Q,R,xmax,xmin,umax,umin)
% * Generates the matrices H, g, Ain, bin to use in a quadratic program
% * from linear state-space MPC parameters. Unconstrained output and
% * control elements can be specified as 'NaN', these do NOT need to be
% * identical between ymax and ymin or umax and umin.

% Size variables
n = size(Ad,1);
m = size(Bd,2);

% First, construct the intermediate cost function matrices AA, BB, QQ, RR
AHat = zeros(n*N,n);
for i = 1:N
    AHat(1+(i-1)*n : n+(i-1)*n,:) = Ad^i;
end

BHat = zeros(n*N,m*N);
for i = 1:N % Skip the first row since its zeros
    for j = 1:i
        BHat(1+(i-1)*n : n+(i-1)*n, 1+(j-1)*m : m+(j-1)*m) = Ad^(i-j)*Bd;
    end
end

% COST MATRICES
QQ = zeros(n*N);
QQ(1:n*(N-1),1:n*(N-1)) = kron(eye(N-1),Q);
QQ(n*(N-1)+1:end,n*(N-1)+1:end) = P;
RR = kron(eye(N),R);

% Now, construct the final Hessian matrix H and vector g
H = (RR + BHat'*QQ*BHat); 
G = BHat'*QQ*AHat; % st f = G*x
W = Q +  AHat'*QQ*AHat; % st. extra cost = x'*W*x
% W =  M'*QQ*M; % st. extra cost = x'*W*x

% Next, generate the associated  constraint matrices 
 % st the constraints are represted by ACon*U <= FCon + LCon*xk
ACon = [BHat; -BHat; eye(m*N); -eye(m*N)];
    
xMaxVec = kron(ones(N,1),xmax);
xMinVec = kron(ones(N,1),xmin);
uMaxVec = kron(ones(N,1),umax);
uMinVec = kron(ones(N,1),umin);
FCon = [xMaxVec; -xMinVec; uMaxVec; -uMinVec];

LCon = [-AHat; AHat; zeros(m*N,n); zeros(m*N,n)];

end

