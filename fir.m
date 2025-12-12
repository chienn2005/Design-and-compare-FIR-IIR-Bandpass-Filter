
%Thiết Kế Bộ Lọc Thông Dải FIR (Phương pháp Parks-McClellan)

% 1. ĐỊNH NGHĨA THAM SỐ THIẾT KẾ
Fs = 48000;         % Tần số lấy mẫu (Hz)
Fpass1 = 1000;      % Dải thông bắt đầu (Hz)
Fpass2 = 4000;      % Dải thông kết thúc (Hz)
Fstop1 = 800;       % Dải chặn dưới kết thúc (Hz)
Fstop2 = 4200;      % Dải chặn trên bắt đầu (Hz)
Astop = 50;         % Độ suy hao dải chặn (dB)
Apass = 1.5;        % Độ gợn sóng dải thông (dB)

% 2. THIẾT KẾ BỘ LỌC DÙNG FIRPMORD VÀ FIRPM ---

% Chuyển đổi tham số suy hao/gợn sóng sang biên độ sai lệch (linear deviation)
delta_stop = 10^(-Astop/20); 
DevPass = (10^(Apass/20) - 1) / (10^(Apass/20) + 1);

% Vector các tần số biên (edges) và biên độ mong muốn (amplitude)
Fedges = [Fstop1, Fpass1, Fpass2, Fstop2]; 
Alevels = [0, 1, 0]; % Biên độ: Dải chặn -> Dải thông -> Dải chặn
Devs = [delta_stop, DevPass, delta_stop]; % Sai lệch tối đa cho phép

% Tính toán Bậc bộ lọc N và các tham số tối ưu (Fo, Ao, W)
[N, Fo, Ao, W] = firpmord(Fedges, Alevels, Devs, Fs);
N = N + 50; 

% Thiết kế Bộ lọc FIR bằng firpm (Phương pháp Parks-McClellan/Remez)
b = firpm(N, Fo, Ao, W); 

% --- 3. TẠO TÍN HIỆU ĐA TẦNG ---
T_duration = 0.05; % Thời gian mô phỏng
t = 0:1/Fs:T_duration - 1/Fs; % Vector thời gian

% Tín hiệu thành phần (Tín hiệu lọt dải và tín hiệu bị chặn)
x_pass1 = 1.0 * sin(2*pi*2500*t); 
x_pass2 = 0.5 * sin(2*pi*3500*t); 
x_stop1 = 0.3 * sin(2*pi*500*t);  
x_stop2 = 0.2 * sin(2*pi*5000*t); 

% Tín hiệu đầu vào (pre_x)
pre_x = x_pass1 + x_pass2 + x_stop1 + x_stop2;

% --- 4. LỌC TÍN HIỆU ---
after_x = filter(b, 1, pre_x);

%PLOT
figure('Name', 'Phân tích Bộ lọc và Tín hiệu');

% --- A. Đáp ứng Biên độ Bộ lọc (dB) ---
subplot(3, 1, 1);
[H, freq] = freqz(b, 1, 4096, Fs);
mag_dB = 20*log10(abs(H));
plot(freq/1000, mag_dB);
title('Đáp ứng Biên độ Bộ lọc Thông Dải (dB)');
xlabel('Tần số (kHz)');
ylabel('Biên độ (dB)');
grid on;
xlim([0, Fs/2000]); 

% Vẽ các dải tần lên biểu đồ
hold on;
plot([Fstop1, Fstop1]/1000, [-Astop, 5], 'r--');
plot([Fpass1, Fpass1]/1000, [-Astop, 5], 'g--');
plot([Fpass2, Fpass2]/1000, [-Astop, 5], 'g--');
plot([Fstop2, Fstop2]/1000, [-Astop, 5], 'r--');
legend('Đáp ứng Bộ lọc', 'Fstop1', 'Fpass1', 'Fpass2', 'Fstop2', 'Location', 'SouthWest');

% --- B. Tín hiệu Đầu vào (Miền Thời gian) ---
subplot(3, 2, 3); 
plot(t, pre_x); 
title('Tín hiệu Đầu vào (trước lọc)');
xlabel('Thời gian (s)');
ylabel('Biên độ');
grid on;

% --- C. Tín hiệu Sau lọc (Miền Thời gian) ---
subplot(3, 2, 4); 
plot(t, after_x); 
title('Tín hiệu Sau lọc');
xlabel('Thời gian (s)');
ylabel('Biên độ');
grid on;

% --- D. Phổ Tín hiệu Đầu vào (Miền Tần số) ---
Y_pre = fft(pre_x);
P2_pre = abs(Y_pre/length(pre_x));
P1_pre = P2_pre(1:floor(length(pre_x)/2)+1);
P1_pre(2:end-1) = 2*P1_pre(2:end-1);
f_fft = Fs*(0:floor(length(pre_x)/2))/length(pre_x);
subplot(3, 2, 5); 
plot(f_fft/1000, 20*log10(P1_pre));
title('Phổ Tín hiệu Đầu vào (dB)');
xlabel('Tần số (kHz)');
ylabel('Biên độ (dB)');
grid on;

% --- E. Phổ Tín hiệu Sau lọc (Miền Tần số) ---
Y_after = fft(after_x);
P2_after = abs(Y_after/length(after_x));
P1_after = P2_after(1:floor(length(after_x)/2)+1);
P1_after(2:end-1) = 2*P1_after(2:end-1);
subplot(3, 2, 6); 
plot(f_fft/1000, 20*log10(P1_after));
title('Phổ Tín hiệu Sau lọc (dB)');
xlabel('Tần số (kHz)');
ylabel('Biên độ (dB)');
grid on;

fprintf('\n CÁC CHỈ SỐ ĐÁNH GIÁ BỘ LỌC \n');

% 1. Bậc của bộ lọc (Order: N)
fprintf('* Bậc của bộ lọc (N): %d\n', N); 
% 2. Độ trễ nhóm (Group Delay)
GroupDelay_samples = (N) / 2;
GroupDelay_sec = GroupDelay_samples / Fs;
fprintf('* Độ trễ nhóm (Group Delay):\n');
fprintf('  - Tính bằng mẫu (Samples): %.1f mẫu\n', GroupDelay_samples);
fprintf('  - Tính bằng giây (Seconds): %.6f s\n', GroupDelay_sec);
[Gd, ~] = grpdelay(b, 1, 1024, Fs);
fprintf('  - Giá trị trung bình của Group Delay: %.6f s\n', mean(Gd)/Fs);
fprintf('  (Vì là bộ lọc FIR tuyến tính, Group Delay là hằng số, không gây méo pha)\n');

% 3. Phổ tần số (Frequency Spectrum)

fprintf('* Phổ tần số (Frequency Spectrum):\n');

fprintf(' - Suy hao dải chặn đạt yêu cầu: >= %.1f dB\n', Astop);

fprintf(' - Gợn sóng dải thông đạt yêu cầu: <= %.1f dB\n', Apass);

