%% Matlab Code for the N-Way OMP algorithm as described in the paper:
% "Computing Sparse Representations of Multidimensional Signals Using Kronecker Bases" by C. Caiafa & A. Cichocki,
% Neural Computation Journal, Vol. 25, No. 1 , pp. 186-220, 2013.
% Find a s-sparse representation of a tensor Y, i.e. 
% Y = X \times_1 D_1 \times_2 D_2...\times_N D_N
% where X(i_1,i_2,...,i_N)=0 for all i_n \notin I_n and |I_n|=s(n)

function [metrics] = KroOMP(A1,A2,A3,y,s,epsilon,x)

Y = reshape(y,[size(A3,1),size(A2,1),size(A1,1)]);
D = {};
D{1} = A3;
D{2} = A2;
D{3} = A1;

[I] = size(Y);
N = size(I,2);

for n = 1:N
    M(n) = size(D{n},2);
end

tic;

norma = norm(reshape(Y,[I(1),prod(I)/I(1)]),'fro');

R = Y; % initial residual
Ind = cell(1,N); % set of indices to select for each mode
% L = cell(1,N); % Triangular matrices used in Cholesky factorization of each mode
Dsub = cell(1,N); % subset of selected atoms per each mode

% auxiliar variables for searching the maximum of a tensor
v = cell(1,N); 
add = cell(1,N); 

ni = zeros(1,N); % number of selected indices in each mode

% Initialization
for n = 1:N
    Ind{n} = [];
%     L{n} = 1;
    ni(n) = 0;
end
X = zeros(M); % coefficients
cond = 0; % this condition is true when ni(n)=s(n) for all n (all nonzero coefficients where computed)
posi = zeros(1,N); % additional index to add

condchange = 1;
error = Inf;
% Main loop where selected indices in each mode are found
while ((condchange && (~cond) && (error > epsilon)))    
    niant = ni;
    proj = abs(double(ttensor(tensor(R),transp(D)))); % multiway correlation between dictionary and residual
    
    % detect the largest value in the last dimension
    for n = 1:N
        [proj,v{n}] = max(abs(proj));
        add{n} = 1;
    end
     
    posi(N) = v{N};
    add{N} = posi(N);
    % start from the last dimension, unfold the largest position till the
    % first dimension
    for n = N-1:-1:1        
        posi(n) = v{n}(add{:});
        add{n} = posi(n);
    end
    
    for n = 1:N
        if ((~ismember(posi(n),Ind{n})) && (ni(n) < s(n))) % to see if 
        Ind{n} = [Ind{n},posi(n)];
        Dsub{n} = [Dsub{n},D{n}(:,posi(n))];
        ni(n) = ni(n) + 1;     
        end
    end
    
    mat = Dsub{N};
    for n = N-1:-1:1
        mat = kron(mat,Dsub{n}); 
    end
    a = pinv(mat)*vec(Y);
      
    Contri = reshape(mat*a,size(Y));
    R = Y - Contri;
    
    error = norm(reshape(R,[I(1),prod(I)/I(1)]),'fro')/norma;
    disp([num2str(ni),'  ', num2str(error)])
    
    cond = 1;
    for n =1:N
        cond = cond && (ni(n) == s(n));
    end
    
    condchange = sum(ni-niant);
%     for n = 1:N
%         Dred{n}(:,Ind{n}(ni(n)))=0;
%     end
end

shape = zeros(1,N);

for i = 1:N
    shape(i) = length(Ind{i});
end

X(Ind{:}) = reshape(a,shape);
x_re = vec(double(X));

time = toc;
erroro = norm(x_re - x)/norm(x)

srr = recover_rate(x_re,x);
metrics={'error',erroro;
         'vector',x_re;
         'support recovery rate',srr;
         'time',time
         };
end

function [D] = transp(D)
N = size(D,2);
for n = 1:N
    D{n} = D{n}';
end
end