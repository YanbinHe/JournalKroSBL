function theta = irs_reflection(dim)
theta = zeros(dim,1);

for i = 1:dim
    theta(i)=(rand >0.5)*2-1;
end

theta = diag(theta);
end

