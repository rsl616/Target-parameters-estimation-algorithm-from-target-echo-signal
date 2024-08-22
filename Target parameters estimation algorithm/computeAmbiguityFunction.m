function [ambiguity, delay, doppler] = computeAmbiguityFunction(signal, fs, maxDoppler, maxDelay)
    % Calculating auto-ambiguity function
    % Inputes:
    % signal - input signal
    % fs -     sampling frequency
    % maxDoppler - maximum Doppler frequency shift (Hz)
    % maxDelay   - maximum time delay (s)
    % Outputs:
    % ambiguity  - value of ambiguity function
    % delay      - value of delay axis
    % doppler    - value of Doppler axis

    N = length(signal);           % signal length
    dopplerShifts = linspace(-maxDoppler, maxDoppler, 2 * N);
    delays = 0:1/fs:maxDelay-1/fs;
    ambiguity = zeros(length(delays), length(dopplerShifts));

    % Calculating the influence of each Doppler shift in advance.
    dopplerMat = exp(1j * 2 * pi * (0:N-1).' / fs * dopplerShifts);

    % Calculating auto-ambiguity function
    for k = 1:length(dopplerShifts)
        % Doppler frequency shift is applied to the original signal
        dopplerSignal = signal .* dopplerMat(:, k).';

        for i = 1:length(delays)
            % Generating a delayed signal
            delayIndex = round(delays(i) * fs);
            if delayIndex > N
                delayedSignal = zeros(1, N);  % If the delay exceeds the signal length, 
                                                      % an all-zero vector is used
            else
                delayedSignal = [zeros(1, delayIndex), signal(1:end - delayIndex)];
            end

            % Truncate or stuff the delayed signal to ensure its length is consistent with N
            %delayedSignal = [delayedSignal, zeros(1, N - length(delayedSignal))];

            % Calculate the auto-ambiguity value under the current delay and Doppler frequency shift
            ambiguity(i, k) = abs(sum(delayedSignal(1:N) .* conj(dopplerSignal)));
        end
    end

    % Return delay and Doppler frequency shift
    delay = delays;
    doppler = dopplerShifts;
end

