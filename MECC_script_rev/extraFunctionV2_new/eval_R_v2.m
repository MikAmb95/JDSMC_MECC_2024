

function [R,R_vecU,R_vecX,delta,L,psi,psi_bar,bound_vec] = eval_R_v2(A,B,controlArg,input,dBar,N)




Q = controlArg.Q;
P = controlArg.P;
R = controlArg.R;

[L,psi,psi_bar,delta,bound_vec] = eval_gap_input(A,B,controlArg,input,dBar,N);


%%% R_vecU

for i = 1:N-1
  R_vecU(i) = delta(i)*(2*L(i)*psi(i)+delta(i));      
end

R_U = norm(R)*sum(R_vecU);

for i = 1:N
  R_vecX(i) = psi_bar(i)*psi_bar(i);      
end

R_X = max(norm(P),norm(Q))*sum(R_vecX);


R = R_U+R_X;