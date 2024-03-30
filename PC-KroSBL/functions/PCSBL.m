function [errorp,time,mu_new,alpha,H_re] = PCSBL(numItr,H,Res1,noise_var,y,K2,irs_pattern,A_irs_a,A_irs_d,A1,A2,H1,H2,SNR,thres) 
% modified based on the PC-SBL


time = 0;

a=0.5;
b=1e-4;

eta=0.9;

D=eye(Res1^3);
alpha_new=ones(Res1^3,1);

% var_new=inv(H'*H/noise_var+D);
diagD = diag(D);
Dinv = diag(diagD.^-1);
invPhic = H'* (noise_var * eye(size(H,1)) + H * Dinv * H')^-1 * H;
var_new = Dinv - Dinv * invPhic * Dinv;

mu_new=1/noise_var*var_new*(H'*y);

itrp=1;
normp = 1;

while(itrp < numItr + 1 && normp > thres)
    itrp;
    
    tic;
    mu_old=mu_new;
    
    mul=[mu_new(2:Res1^3);0];
    mur=[0;mu_new(1:Res1^3-1)];
    var=real(diag(var_new));
    varl=[var(2:Res1^3);0];
    varr=[0;var(1:Res1^3-1)];
    E=abs(mu_new).^2+eta*abs(mul).^2+eta*abs(mur).^2+var+eta*varl+eta*varr;
    alpha_new=(a+0.5)./(0.5*E+b);
    
    
    idx1=find(alpha_new>1e3*min(alpha_new));
    alpha_new(idx1)=1e10;
    
    
    alf=[alpha_new(2:Res1^3); 0];                                %   left-shifted version
    arf=[0; alpha_new(1:Res1^3-1)];                              %   right-shifted version
    diagD = alpha_new+eta*alf+eta*arf;
    

    
    % Modified: this is to reduce the dimensionality of matrix inversion
    Dinv = diag(diagD.^-1);
    invPhic = H'* (noise_var * eye(size(H,1)) + H * Dinv * H')^-1 * H;
    var_new = Dinv - Dinv * invPhic * Dinv;
    
%     var_new=inv(H'*H/noise_var+diag(diagD));
    mu_new=1/noise_var*var_new*(H'*y);
    
    time = time + toc;
    
    
    normp = norm(mu_new - mu_old)/norm(mu_old);
    
    itrp=itrp+1;
end

for i = 1:K2
    Htt = irs_pattern(:,i).'*kr(A_irs_a.',A_irs_d').';
    Ht(:,:,i) = Htt(:,1:Res1);
    H_re(:,:,i) = kron(kron(Ht(:,:,i),conj(A1)),A2)*mu_new;
    Htrue(:,:,i) = vec(H2*diag(irs_pattern(:,i))*H1);
    nmse(i) = (norm(H_re(:,:,i) - Htrue(:,:,i),'fro')/norm(Htrue(:,:,i),'fro'))^2;
end

errorp = sum(nmse(:))/K2
alpha = 1./alpha_new;
end

