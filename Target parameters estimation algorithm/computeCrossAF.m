function [CAF, tau, fd] = computeCrossAF(sig1, sig2, fs, maxDoppler, maxDelay,tstart)
    % Calculate the cross-ambiguity function between two signals
    % Inputs:
    % sig1 - emitted signal, sig2 - echo signal
    % fs   - sampling frequency
    % maxDoppler - maximum Doppler frequency shift (Hz)
    % maxDelay   - maximum time delay (s)
    % Outputs:
    % CAF - value of Cross-Ambiguity Function
    % tau - value of delay axis(s)
    % fd -  value of Doppler axis(Hz)
    % tstart - starting time
    N = length(sig2);              % Define N using the length of sig2
    sig1 = [sig1, zeros(1, N - length(sig1))]; % Ensure that the length of sig1 is at least N
    sig1 = sig1(1:N);                          % Cut or fill to match the length

    fd = linspace(-maxDoppler, maxDoppler,  2*N);
    tau = tstart:1/fs:maxDelay-1/fs; 
    CAF = zeros(length(tau), length(fd));


    % Calculating Doppler frequency shift matrix
    dopplerMat = exp(-1j*2*pi.*(0:N-1)'/fs.*fd);
    % Calculate the cross-ambiguity function
    for i = 1:length(tau)
        delayIndex = round(tau(i) * fs);
        if delayIndex >= N
            delayedSig1 = zeros(1, N);  % Delay exceeds signal length, using all zero vector
        else
            delayedSig1 = [zeros(1, delayIndex), sig1(1:N - delayIndex)]; % Generate a delayed signal
        end

        % Calculate CAF under all Doppler frequency shifts
        for j = 1:length(fd)
            dopplerShiftedSig1 = conj(delayedSig1) .* dopplerMat(:,j).';
            CAF(i, j) = abs(sig2 * dopplerShiftedSig1.');
        end
    end
end
