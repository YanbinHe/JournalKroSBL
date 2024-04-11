function [mean,var] = posterior_compute(sigmasq,gamma1,gamma2,gamma3,A,B,C,Ht,y)
% 

gamma1sqrt = diag(sqrt(gamma1));
gamma2sqrt = diag(sqrt(gamma2));
gamma3sqrt = diag(sqrt(gamma3));

Abar = gamma1sqrt*(A'*A)*gamma1sqrt;
Bbar = gamma2sqrt*(B'*B)*gamma2sqrt;
Cbar = gamma3sqrt*(C'*C)*gamma3sqrt;

[Q1,L1] = eig(Abar);
[Q2,L2] = eig(Bbar);
[Q3,L3] = eig(Cbar);

Par1 = kron(kron(gamma1sqrt*Q1,gamma2sqrt*Q2),gamma3sqrt*Q3);
Par2 = kron(kron(diag(L1),diag(L2)),diag(L3));

diaginv = (sigmasq^-1*Par2+ones(size(Par2,1),1)).^-1;

var = Par1*(diaginv.*Par1');
mean = sigmasq^-1*var*(Ht'*y);
end