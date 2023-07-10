function array_res = generate_steering(dim,angle)
% generate the array response of the antenna array
% dim: dimension of the vector
% angle: AoA or AoD, from [-1,1] on the discrete grid
array_res = zeros(dim,length(angle));
for i=1:dim
    array_res(i,:)=sqrt(1/dim)*exp(1j*2*pi*0.5*(i-1)*angle);
end

end

