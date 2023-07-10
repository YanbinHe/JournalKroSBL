function [mean,var] = posterior_compute_inv(sigmasq,mu1,mu2,mu3,A,B,C,Ht,y)
% 
mu1inv = mu1.^-1;
mu2inv = mu2.^-1;
mu3inv = mu3.^-1;


gA = (mu1inv.*A');
gB = (mu2inv.*B');
gC = (mu3inv.*C');

Dinv = diag(kron(mu1inv,kron(mu2inv,mu3inv)));

Abar = A*gA;
Bbar = B*gB;
Cbar = C*gC;

[U1,S1] = eig(Abar);
[U2,S2] = eig(Bbar);
[U3,S3] = eig(Cbar);

diphihu = kron((gA*U1),kron(gB*U2,gC*U3));
P = kron(diag(S1),kron(diag(S2),diag(S3)));

diaginv = (sigmasq*ones(size(P,1),1)+P).^-1;
hy = Ht'*y;

ginv = (diphihu*(diaginv.*diphihu'));
var = Dinv - ginv;
mean = sigmasq^-1*var*hy;
end