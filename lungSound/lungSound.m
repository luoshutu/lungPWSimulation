%% 目的：实现肺部信号的解调
% By Luoshutu
clear; clc;
%% 参数设置
f0          = 3.5e6;    % 换能器中心频率 [Hz]
fs          = 40e6;     % 采样频率
prf         = 1000;     % 脉冲重复频率Hz/s

%% 加载数据
load data_m

%% IQ解调
[mm,nn]  = size(data_m);
Data_I = zeros(mm,nn);
Data_Q = zeros(mm,nn);
t = 0:nn-1;
for i = 1:mm
    Data_I(i,:) = cos(2*pi*f0/fs*t).*(data_m(i,:));
    Data_Q(i,:) = sin(2*pi*f0/fs*t).*(data_m(i,:));
end

%% 显示M模式
hd   = design(fdesign.lowpass('N,F3dB',16,5e6,100e6),'butter');
Data_I_fil = filter(hd,Data_I.');
Data_Q_fil = filter(hd,Data_Q.');
Data_Amp = sqrt((Data_I_fil).*(Data_I_fil) + (Data_Q_fil).*(Data_Q_fil));

figure;
imagesc(20*log(1 + abs(Data_Amp)));
title('M模式');
colormap(gray);

%% 相位解调
xx = 0:(1/prf):(1/prf)*(mm-1);
z  = hilbert(sum(real(data_m(:,400:420)')));                 %希尔伯特变换对的利用---通过实部来求虚部
% z = sum(data_m(:,400:420)');
x1 = z.*exp(-j*2*pi*prf*xx);
% x1    = Data_e.*exp(-j*2*pi*fc*t);
pha   = angle(x1);          % 求取复数信号的相位
demod = unwrap(pha);        % 相位校正  
hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,1e4,20e6,100e6),'butter');
demod_fil = filter(hd,demod);
figure;
plot((demod_fil));
%% 求arctan2
for xy = 1:mm
    De(xy) = atan2(sum(real(Data_I(xy,400:420))),sum(real(Data_Q(xy,400:420))));
end
Dt = unwrap(De);
hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,1e4,10e6,100e6),'butter');
Dt_fil = filter(hd,Dt);
figure
plot(Dt,'b');
hold on;
plot(Dt_fil,'r');
legend('滤波前','滤波后');

