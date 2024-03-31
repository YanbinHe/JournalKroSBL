% noisy compare
lenK = length(K);
lenM = length(M1);
lenS = length(SNR);

for avg = 1:AVG

    
    % three different scenarios for (a) (b) (c) 
    % not doing all the combinations to reduce the computational time, may result
    % in mismatch value in final results

    % this is for major revision purpose

    % fig 1/2 (b) nmse and srr vs measurement      
    k = 2;
    s = 5;
    m = 3;

    x1 = zeros(N,1);
    x2 = zeros(N,1);
    x3 = zeros(N,1);    
    
    % generate measuring dictionaries
    A1_ori = randn(M1(m),N);
    A2_ori = randn(M2,N);
    A3_ori = randn(M3,N);

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

    uNsq_simulationM;

    %% save data for each trial
    filename = ['./results/unknown_nsq_time_try2_', num2str(avg),'.mat'];
    save(filename, 'resultM')
end
