
function L = eval_L(H,G)

term1 = norm(H^(-0.5));
term2 = norm((H^(-0.5))*G);

L = term1*term2;