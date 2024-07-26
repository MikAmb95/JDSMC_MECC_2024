function [L,psi,bar_psi] = eval_L_psi(controlArgs,input,N)


alpha = controlArgs.alpha;
beta_vec = controlArgs.betaVec;
lambda = controlArgs.lambda;
gamma = controlArgs.gamma;
P = controlArgs.P;

Q = controlArgs.Q;

beta = beta_vec(1);

sqrtQ = 1/sqrt(min(eig(P)));


gamma_k = [];

gamma_k = [gamma_k;0];


for i = 1:N
L(i) = eval_L(input.matrices(i).H,input.matrices(i).G);

psi1 = (N-i)*gamma+lambda*alpha-beta;
psi2 = lambda*alpha;
psi(i) = sqrt(max(psi1,psi2))*sqrtQ;
 
bar_psi(i) = sqrt((N-i)*gamma+lambda*alpha-beta_vec(i))*sqrtQ;


end

psi(N) = sqrt(alpha)/sqrtQ;
bar_psi(N) = sqrt(alpha)/sqrtQ;


end