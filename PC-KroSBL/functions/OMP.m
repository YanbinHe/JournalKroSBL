function [er,time,supp,H_re] = OMP(Res1,H,klevel,irs_pattern,A_irs_a,A_irs_d,A1,A2,H1,H2,y,K2)

time = 0;

tic;

residual = y;
set = [];

ite = 1;
while ite < klevel+1
    re(ite) = norm(residual);
    inner_pro = abs(H'*residual);
    [~,in] = max(inner_pro);
    set = [set in];
    
    Hset = H(:,set);
    P = Hset*(Hset'*Hset)^-1*Hset';
    
    residual = (eye(size(P))-P)*y;
    ite = ite + 1;
end

Hset = H(:,set);
coeff = pinv(Hset)*y;

mu_xx = zeros(size(H,2),1);
mu_xx(set) = coeff;

time = time + toc;

supp = zeros(size(H,2),1);
supp(set) = 1;

for i = 1:K2
    Htt = irs_pattern(:,i).'*kr(A_irs_a.',A_irs_d').';
    Ht(:,:,i) = Htt(:,1:Res1);
    H_re(:,:,i) = kron(kron(Ht(:,:,i),conj(A1)),A2)*mu_xx;
    Htrue(:,:,i) = vec(H2*diag(irs_pattern(:,i))*H1);
    er(i) = norm(H_re(:,:,i) - Htrue(:,:,i))/norm(Htrue(:,:,i));
end

er = mean(er);
end

