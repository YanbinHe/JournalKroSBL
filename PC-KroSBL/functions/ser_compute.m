function [ser] = ser_compute(H_est,H_true,SNR,X,symbolmatrix,ms_ante,bs_ante,K2)

% H_est: estimated channel matrix (in the form of column vector)
% H_true: the true channel matrix (in the form of column vector)
% SNR: scenario snr
% X: transmitted symbols (row vector)
% inPut: the true bits stream for computing the Bit Error Rate (row vector)
% symbolmatrix: codebook

% We first generate a certain number of QPSK symbols x according to the
% variable sym. Then we transmit these symbols through the channel by
% multiplying them and obtain y. The detector T is designed by Linear MMSE 
% detector. The received symbols y_mmse are obtained by Ty. To reterive the
% symbols, we project the received symbols on the discrete constellation by
% computing the Euclidean distance.


ser = 0;
time_slots = length(X);

for i = 1:K2
    Htrue = reshape(H_true(:,:,i),bs_ante,ms_ante);
    Hest = reshape(H_est(:,:,i),bs_ante,ms_ante);
    
    % generate beamformer
    [s,~,d] = svd(Hest);
    precoder = d(:,1);
    decoder = s(:,1);
    
    SNRl = 10^(SNR/10);
 
    Yo = decoder'*Htrue*precoder*X;
    
    signal_power = norm(vec(Yo))^2 / numel(Yo);
    sigmasq = signal_power/SNRl;
    noise = sqrt(sigmasq / 2)*(randn(size(Yo))+1i*randn(size(Yo)));
    Y = Yo + noise;
    
    error = 0;
    for j = 1:time_slots % for all time slots
        y = Y(j); % for each received signal
        for k = 1:4
            norm_dis(k) = norm(y - decoder'*Hest*precoder*symbolmatrix(k));
        end
        [~,idx] = min(norm_dis);
        error = error + 1*(X(j) ~= symbolmatrix(idx));
    end
    
    ser = ser + error/numel(X);
end
ser = ser/K2;
end
