

function [L,psi,bar_psi,delta_out,bound_vec] = eval_gap_input(A,B,controlArg,input,dBar,N)


delta(1) = input.eBarVec(1);

bound_vec = zeros(3,1);


[L,psi,bar_psi] = eval_L_psi(controlArg,input,N);

gamma_k = [];

gamma_k = [gamma_k;0];

for k = 2:N
    
 gamma_k = [gamma_k;L(k)*(bar_psi(k)+psi(k))];

end
%gamma_k = eval_gamma_k2(controlArg,input,N);


for i = 2:N-1
l = input.ellBarVec(i);
H = input.matrices(i).H;
tau = input.omegaVec(i);
sigma = input.sigmaVec(i);
e = input.eBarVec(i-1);

kH = max(eig(H))/min(eig(H)); 
skH = sqrt(kH);
%First term \bar{tau}_k * \bar(e}_k-1 
% where \bar{tau}_k = (eta_k(lk) + tau_k)

eta1 = skH; 
eta2 = (1-(kH^-0.5))^((l-1)/2);
eta = eta1*eta2;
bound_vec(i,1) = eta*(1 + tau)*e;

%Second term sigma_k * \bar{d}

bound_vec(i,2) = eta*sigma*dBar;

%Third term gamma_g


bound_vec(i,3) = gamma_k(i);


%%% sum of the three bounds

delta(i) = sum(bound_vec(i,:));

end

delta_out = delta;
end