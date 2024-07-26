function [X] = reorderMPCVector(x,n,N)
% Reorders an MPC output that is N*n x 1 into a (n,N) matrix X
X = zeros(n,N);
if length(x) ~= n*N
    error('Mismatch in length of x and (n,N)')
end
for i = 1:N
    X(:,i) = x(1+n*(i-1):n+n*(i-1));
end
end

