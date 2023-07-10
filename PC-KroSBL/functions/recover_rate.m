function r = recover_rate(x,xtilde)
% computing thr support recovery support

% using the equation defined in "Alternative to Extended Block Sparse 
% Bayesian Learning and Its Relation to Pattern-Coupled Sparse Bayesian 
% Learning"

delta = 1e-4;

suppx = find(abs(x) > delta*max(abs(x)));
suppxt= find(abs(xtilde) > delta*max(abs(xtilde)));

common = intersect(suppx,suppxt);
nomi = length(common);

suppx_d = suppx;
suppxt_d = suppxt;

[~,idx] = ismember(common,suppx);
suppx_d(idx) = [];
[~,idxt]= ismember(common,suppxt);
suppxt_d(idxt) = [];

diff = length(suppx_d) + length(suppxt_d);
deno = diff + length(suppx);

r = nomi/deno;
end

