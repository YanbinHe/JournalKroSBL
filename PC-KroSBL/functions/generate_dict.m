function dict = generate_dict(dim,NumResolu)

dict = zeros(dim,NumResolu);

for i=1:NumResolu
%     dict(:,i)=(sqrt(1/dim)*exp(1j*pi*cos(pi*(i-1)/NumResolu)*[0:dim-1])).';
    dict(:,i)=(sqrt(1/dim)*exp(1j*pi*(2*(i-1)/NumResolu-1)*[0:dim-1])).'; % from -1 to 1, uniformlly sampled
end 
end

