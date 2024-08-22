%%%%
% Target parameters estimation based on time delay-Doppler deconvolution
% from the target echo signal
% Involving the following:
% Doppler frequency shift
% Auto-Ambiguity Function
% Cross-Ambiguity Function
% R-L(Richardson-Lucy) deconvolution
% illustration:
% The transmitting signal adopts PCW signal with hamming window
%%%%
%%
clc
clear
close all
fs = 20e3;      % sampling frequency 50kHz
T = 0.1;
B = 1/T;
p_wid = 0.02;   % pulse width
t = 0:1/fs:T-1/fs;
t1 = 0:1/fs:p_wid-1/fs;
A_pcw = 1;
f_pcw = 7e3;    % center frequency 7kHz
I = 4;          % pulse number
c = 1500;       % speed of sound
R = [125,140,153,160];  % range
v = [5,-8,8,-10];       % speed
rou = [0.8,1,0.7,0.6];
k_it = [10 20 40 60];   % iteration times of deconvolution algorithm
v_dopl = 2*f_pcw/c.*v; % Doppler frequency shift
%r_signal = zeros(I,length(t1));
r_signal_p = zeros(I,round(2*max(R)*fs/c)+length(t1));
y_pcw1 = A_pcw*sin(2*pi*f_pcw*t1);
y_pcw = [y_pcw1,zeros(1,length(t)-length(t1))];    % PCW signal
% temp = zeros(I,length(t1));
% real_echo = zeros(1,length(r_signal_p));
% imag_echo = zeros(1,length(r_signal_p));
w_han = hann(length(y_pcw1))';
s_s = w_han.*y_pcw1;   % emitted signal
% the echo signal is obtained as follows
for i = 1:I
    r_signal_p(i,(round(2*R(i)*fs/c)+1):(round(2*R(i)*fs/c)+length(t1))) = rou(i).*s_s.*exp(1j*2*pi*v_dopl(i).*t1);
end
r_s = sum(r_signal_p,1);
t_p = (0:1:length(r_s)-1)/fs;

N = 2^nextpow2(fs*p_wid);   % DFT points
Y = fft(s_s,N)/N*2;         % transmitting signal spectrum
f = fs/N*(0:1:N-1);         % frequency axis
P = abs(Y);                 % transmitting signal magnitude spectrum
%%
maxDelay1 = 0.1;            % maximum delay of AAF
maxDoppler1 = 190;          % maximum Doppler frequency of AAF
maxDelay2 = 0.25;           % maximum delay of CAF
maxDoppler2 = 190;          % maximum Doppler frequency of CAF
tstart = 0.12;              % starting delay
[a_fmag,delay_a,dopp_a] = computeAmbiguityFunction(s_s,fs,maxDoppler1,maxDelay1);
% pay attention to the order of two signals in the cross-ambiguity function
% the first signal is the transmission signal
% the second signal is the echo signal
[c_fmag,delay_c,dopp_c] = computeCrossAF(s_s, r_s, fs, maxDoppler2, maxDelay2,tstart);
% inversion filling for drawing
a_fmag_d = [fliplr(a_fmag'),a_fmag'];
delay_a_d = [-fliplr(delay_a),delay_a];
a_nom_max = 0;
for i = 1:length(a_fmag_d(:,1))
    a_nom_max = max(a_nom_max,max(a_fmag_d(i,:)));
end
a_fmag_d_nom = a_fmag_d./a_nom_max;
% the cross-ambiguity function does not need to be filled, but only normalized
c_nom_max = 0;
for i = 1:length(c_fmag(:,1))
    c_nom_max = max(c_nom_max,max(c_fmag(i,:)));
end
c_fmag_nom = c_fmag'./c_nom_max;
% dopp_a_n = -dopp_a;
% dopp_c_n = -dopp_c;
row_c = size(c_fmag_nom,1);
col_c = size(c_fmag_nom,2);
fsmd = zeros(row_c,col_c,4);
fsmd_nom = fsmd;
% no dB
% db_afmag = 20*log10(a_fmag_d_nom);
% db_cfmag = 20*log10(c_fmag_nom);
%%
% R-L deconvolution
new_psf = imresize(a_fmag_d_nom,[row_c col_c]);
for i = 1:4
    fsmd(:,:,i) = deconvlucy(c_fmag_nom,new_psf,k_it(i));
    fsmd_nom_max = 0;
    for j = 1:row_c
        fsmd_nom_max = max(fsmd_nom_max,max(fsmd(j,:,i)));
    end
    fsmd_nom(:,:,i) = fsmd(:,:,i)./fsmd_nom_max;
end
% fsmd = deconvlucy(c_fmag_nom,new_psf,10);
% fsmd_nom_max = 0;
% for j = 1:row_c
%     fsmd_nom_max = max(fsmd_nom_max,max(fsmd(j,:)));
% end
% fsmd_nom = fsmd./fsmd_nom_max;
% figure
% imagesc(delay_c*c/2, dopp_c*c/(2*f_pcw), abs(fsmd_nom));
% xlabel('距离 (m)');
% ylabel('速度 (m/s)');
% title('k = 10 反射性密度函数');
% colormap jet
% colorbar;
%%
% Auto-Cross-AF plot
figure;
imagesc(delay_a_d*c/2, dopp_a*c/(2*f_pcw), abs(a_fmag_d_nom));
xlabel('距离 (m)');
ylabel('速度 (m/s)');
title('自模糊度函数');
colormap jet
colorbar;
figure;
imagesc(delay_c*c/2, dopp_c*c/(2*f_pcw), abs(c_fmag_nom));
xlabel('距离 (m)');
ylabel('速度 (m/s)');
title('互模糊度函数');
colormap jet
colorbar;
%%
%R-L deconvolution
figure
imagesc(delay_c*c/2, dopp_c*c/(2*f_pcw), abs(fsmd_nom(:,:,1)));
xlabel('距离 (m)');
ylabel('速度 (m/s)');
title('k = 10 反射性密度函数');
colormap jet
colorbar;
figure
imagesc(delay_c*c/2, dopp_c*c/(2*f_pcw), abs(fsmd_nom(:,:,2)));
xlabel('距离 (m)');
ylabel('速度 (m/s)');
title('k = 20 反射性密度函数');
colormap jet
colorbar;
figure
imagesc(delay_c*c/2, dopp_c*c/(2*f_pcw), abs(fsmd_nom(:,:,3)));
xlabel('距离 (m)');
ylabel('速度 (m/s)');
title('k = 40 反射性密度函数');
colormap jet
colorbar;
figure
imagesc(delay_c*c/2, dopp_c*c/(2*f_pcw), abs(fsmd_nom(:,:,4)));
xlabel('距离 (m)');
ylabel('速度 (m/s)');
title('k = 60 反射性密度函数');
colormap jet
colorbar;
%%
figure
plot(t_p,real(r_s));
xlabel('time/s'); ylabel('amplitude/v'); title('receive signal')
figure
subplot(211)
plot(t1,s_s);
xlabel('time/s'); ylabel('amplitude/v'); title('emit signal');
subplot(212)
plot(f(1:N/2+1),P(1:N/2+1));
xlabel('频率(Hz)');  ylabel('幅度');  title('频谱');
%%
% figure
% imagesc(delay*c/2, dopp*c/(2*f_pcw), abs(afmag_nom));
% xlabel('距离 (m)');
% ylabel('速度 (m/s)');
% title('normalize 自模糊度函数');
% colorbar;
% figure;
% imagesc(delay1*c/2, dopp1*c/(2*f_pcw), abs(cfmag_nom));
% xlabel('距离 (m)');
% ylabel('速度 (m/s)');
% title('normalize 互模糊度函数');
% colorbar;
% %%
% figure
% imagesc(delay*c/2, dopp*c/(2*f_pcw), abs(db_afmag));
% xlabel('距离 (m)');
% ylabel('速度 (m/s)');
% title('db 自模糊度函数');
% colorbar;
% figure;
% imagesc(delay1*c/2, dopp1*c/(2*f_pcw), abs(db_cfmag));
% xlabel('距离 (m)');
% ylabel('速度 (m/s)');
% title('db 互模糊度函数');
% colorbar;