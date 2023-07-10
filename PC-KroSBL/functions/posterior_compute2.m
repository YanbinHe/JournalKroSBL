function [mean,var] = posterior_compute2(sigmasq,gamma1,gamma2,gamma3,A,B,C,Ht,y)
% 

gA = (gamma1.*A');
gB = (gamma2.*B');
gC = (gamma3.*C');

G = diag(kron(gamma1,kron(gamma2,gamma3)));

Abar = A*gA;
Bbar = B*gB;
Cbar = C*gC;

[U1,P1] = eig(Abar);
[U2,P2] = eig(Bbar);
[U3,P3] = eig(Cbar);

ghhu = kron((gA*U1),kron(gB*U2,gC*U3));
P = kron(diag(P1),kron(diag(P2),diag(P3)));

diaginv = (sigmasq*ones(size(P,1),1)+P).^-1;
hy = Ht'*y;

ginv = (ghhu*(diaginv.*ghhu'));
var = G - ginv;
mean = sigmasq^-1*var*hy;
end