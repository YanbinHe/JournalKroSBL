function [mu,Sigma] = posterior_compute_noiseless(gamma1,gamma2,gamma3,A1,A2,A3,y1,y2,y3)

Gamma1sqrt = diag(gamma1.^0.5);
Gamma2sqrt = diag(gamma2.^0.5);
Gamma3sqrt = diag(gamma3.^0.5);

temp1 = Gamma1sqrt*pinv(A1*Gamma1sqrt);
temp2 = Gamma2sqrt*pinv(A2*Gamma2sqrt);
temp3 = Gamma3sqrt*pinv(A3*Gamma3sqrt);

mu = kron(temp1*y1,kron(temp2*y2,temp3*y3));
Sigma = diag(kron(gamma1,kron(gamma2,gamma3))) - kron(temp1*A1*diag(gamma1),kron(temp2*A2*diag(gamma2),temp3*A3*diag(gamma3)));

end

