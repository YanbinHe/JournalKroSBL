function r = recover_rate(x_re,xt)
% computing thr support recovery support

% using the equation defined in "Alternative to Extended Block Sparse 
% Bayesian Learning and Its Relation to Pattern-Coupled Sparse Bayesian 
% Learning"

delta = 1e-3;

suppxre = find(abs(x_re) > delta*max(abs(x_re))); % to eliminate the potential small terms that come from estimation inaccuracy
suppxt = find(abs(xt) > 0);

nomi = length(intersect(suppxre,suppxt));
deno = length(union(suppxre,suppxt));% the union of two sets

r = nomi/deno;
end

