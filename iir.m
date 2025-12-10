% ========================================================================
% IIR BANDPASS BUTTERWORTH 
%   Fs      = 48 kHz
%   Fpass   = [1  4] kHz
%   Fstop   = < 0.8 kHz và > 4.2 kHz
%   Astop   = 60 dB
%   Apass   = 1  dB
% ========================================================================

clc; clear; close all;

%% 1. THÔNG SỐ THIẾT KẾ
Fs     = 48000;      % Hz
Fpass1 = 1000;       % Hz
Fpass2 = 4000;       % Hz
Fstop1 = 800;        % Hz
Fstop2 = 4200;       % Hz

Astop = 60;         % dB - suy hao dải chặn
Apass = 1;          % dB - gợn sóng dải thông

% Chuẩn hóa tần số (0–1)
Wp = [Fpass1 Fpass2]/(Fs/2);   % dải thông
Ws = [Fstop1 Fstop2]/(Fs/2);   % dải chặn

%% 2. THIẾT KẾ BỘ LỌC IIR BUTTERWORTH (DẠNG SOS CHO ỔN ĐỊNH)

% Tính bậc và tần số cắt
[N, Wc] = buttord(Wp, Ws, Apass, Astop);

% Thiết kế Butterworth bandpass - dùng ZPK rồi chuyển sang SOS
[z, p, k] = butter(N, Wc, 'bandpass');
[sos, g]  = zp2sos(z, p, k);      % dạng biquad + gain

fprintf('Bậc bộ lọc IIR Butterworth bandpass: N = %d\n', N);

%% 3. TẠO TÍN HIỆU MULTITONE KIỂM THỬ

T_duration = 0.05;                         % 50 ms
t = 0:1/Fs:T_duration - 1/Fs;

% Trong dải thông 1–4 kHz
x_pass1 = 1.0 * sin(2*pi*2500*t);          % 2.5 kHz
x_pass2 = 0.5 * sin(2*pi*3500*t);          % 3.5 kHz

% Ngoài dải thông
x_stop1 = 0.3 * sin(2*pi*500*t);           % 500 Hz
x_stop2 = 0.2 * sin(2*pi*6000*t);          % 6 kHz

pre_x = x_pass1 + x_pass2 + x_stop1 + x_stop2;

%% 4. LỌC TÍN HIỆU

after_x = sosfilt(sos, pre_x);            % lọc với SOS
after_x = g * after_x;                    % nhân thêm gain

%% 5. PHÂN TÍCH ĐÁP ỨNG & TÍN HIỆU

figure('Name','Phân tích Bộ lọc IIR và Tín hiệu');

% --- A. Đáp ứng biên độ bộ lọc ---
subplot(3,1,1);
[H, f] = freqz(sos, 4096, Fs);            % dùng SOS
H = g * H;
mag_dB = 20*log10(abs(H)+eps);
plot(f/1000, mag_dB);
title('Đáp ứng Biên độ Bộ lọc IIR Thông Dải (Butterworth)');
xlabel('Tần số (kHz)');
ylabel('Biên độ (dB)');
grid on;
xlim([0 Fs/2000]);   % 0–24 kHz

hold on;
plot([Fstop1 Fstop1]/1000, [-Astop 5], 'r--');
plot([Fpass1 Fpass1]/1000, [-Astop 5], 'g--');
plot([Fpass2 Fpass2]/1000, [-Astop 5], 'g--');
plot([Fstop2 Fstop2]/1000, [-Astop 5], 'r--');
legend('Đáp ứng Bộ lọc', 'Fstop1', 'Fpass1', 'Fpass2', 'Fstop2', ...
       'Location','SouthWest');

% --- B. Tín hiệu trước & sau lọc (thời gian) ---
subplot(3,2,3);
plot(t, pre_x);
title('Tín hiệu Đầu vào (trước lọc)');
xlabel('Thời gian (s)');
ylabel('Biên độ'); grid on;

subplot(3,2,4);
plot(t, after_x);
title('Tín hiệu Sau lọc (IIR Bandpass)');
xlabel('Thời gian (s)');
ylabel('Biên độ'); grid on;

% --- C. Phổ trước lọc ---
Y_pre = fft(pre_x);
L = length(pre_x);
P2_pre = abs(Y_pre/L);
P1_pre = P2_pre(1:floor(L/2)+1);
P1_pre(2:end-1) = 2*P1_pre(2:end-1);
f_fft = Fs*(0:floor(L/2))/L;

subplot(3,2,5);
plot(f_fft/1000, 20*log10(P1_pre+eps));
title('Phổ Tín hiệu Đầu vào (dB)');
xlabel('Tần số (kHz)');
ylabel('Biên độ (dB)');
grid on; xlim([0 Fs/2000]);

% --- D. Phổ sau lọc ---
Y_after = fft(after_x);
P2_after = abs(Y_after/L);
P1_after = P2_after(1:floor(L/2)+1);
P1_after(2:end-1) = 2*P1_after(2:end-1);

subplot(3,2,6);
plot(f_fft/1000, 20*log10(P1_after+eps));
title('Phổ Tín hiệu Sau lọc (dB)');
xlabel('Tần số (kHz)');
ylabel('Biên độ (dB)');
grid on; xlim([0 Fs/2000]);

%% 6. CHỈ SỐ ĐÁNH GIÁ

fprintf('\nCÁC CHỈ SỐ ĐÁNH GIÁ BỘ LỌC IIR BANDPASS\n');

% 1. Độ trễ nhóm (trung bình trong dải thông)
[Gd, f_gd] = grpdelay(sos, 1024, Fs);     % Gd: theo mẫu
idx_pass = (f_gd >= Fpass1) & (f_gd <= Fpass2);
Gd_mean_samples = mean(Gd(idx_pass));
Gd_mean_sec     = Gd_mean_samples / Fs;

fprintf('* Độ trễ nhóm trung bình: %.2f mẫu ~ %.6f s\n', ...
        Gd_mean_samples, Gd_mean_sec);

% 2. Suy hao dải chặn & gợn sóng dải thông
pass_idx = (f >= Fpass1) & (f <= Fpass2);
stop_idx = (f <= Fstop1) | (f >= Fstop2);

pass_ripple = max(mag_dB(pass_idx)) - min(mag_dB(pass_idx));  % dB
stop_atten  = -max(mag_dB(stop_idx));                         % dB

fprintf('* Suy hao dải chặn đo được  ~ %.2f dB (yêu cầu >= %.1f dB)\n', ...
        stop_atten, Astop);
fprintf('* Gợn sóng dải thông đo được ~ %.2f dB (yêu cầu <= %.1f dB)\n', ...
        pass_ripple, Apass);

max_pre   = max(abs(pre_x));
max_after = max(abs(after_x));
fprintf('* Max|x(n)| = %.3f, Max|y(n)| = %.3f\n', max_pre, max_after);
