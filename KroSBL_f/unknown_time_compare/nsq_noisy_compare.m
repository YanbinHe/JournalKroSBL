% noisy compare
lenM = length(M1);

for avg = 1:AVG
    % this is for major revision purpose
    
    % generate measuring dictionaries
    A1_ori = randn(M1(end),N);
    A2_ori = randn(M2,N);
    A3_ori = randn(M3,N);
    
    % fig 1/2 (a) nmse and srr vs snr     
    k = 2;
    s = 5;

    x1 = zeros(N,1);
    x2 = zeros(N,1);
    x3 = zeros(N,1);    
    
    supp1 = randsample(1:N, K(k));
    x1(supp1) = randn(K(k),1);
    supp2 = randsample(1:N, K(k));
    x2(supp2) = randn(K(k),1);
    supp3 = randsample(1:N, K(k));
    x3(supp3) = randn(K(k),1);

    x = kron(x1,kron(x2,x3));
    suppTrue = find(abs(x)>0);

    y1 = A1_ori*x1;
    y2 = A2_ori*x2;
    y3 = A3_ori*x3;

    y_ori = kron(y1,kron(y2,y3));

    signal_power = norm(y_ori)^2/length(y_ori); %average signal power per asymbol
    
    for m = 1:lenM
        uNsq_simulationM;
    end

    %% save data for each trial
    filename = ['./results/unknown_nsq_nmse_m_', num2str(avg),'.mat'];
    save(filename, 'resultM')
end
