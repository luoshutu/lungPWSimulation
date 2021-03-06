%% 调相信号的解调
%% 清理
clear
clc
close all

%% 参数
f0 = 2;     % 原始信号频率，单位Hz
fc = 1e3;   % 载波信号频率，单位Hz
fs = 1e5;   % 信号采样频率，单位Hz

%% 设置信号
t           = 0:(1/fs):8;                   % 时间
sigSrc      = 0.5*sin(2*pi*f0.*t) + 0.5*sin(2*pi*2*f0.*t) + 0.3*cos(2*pi*0.5*f0.*t+pi/4);       % 原始信号
sigCar      = sin(2*pi*fc.*t);              % 载波信号
sigPhaMod	= sin(2*pi*fc.*t + 100*sigSrc);     % 调相信号

figure;
plot(sigSrc,'r');

%% IQ解调
% sigDemod = asin(sigPhaMod) - asin(sigCar);
sigDemod = sigPhaMod;

Data_I = cos(2*pi*fc/fs.*t).*(sigDemod);
Data_Q = sin(2*pi*fc/fs.*t).*(sigDemod);

Data_e = Data_I + i.*Data_Q;
% %% 滤波
% hd   = design(fdesign.lowpass('N,F3dB',16,100,fs),'butter');
% Data_I_fil = filter(hd, Data_I);
% Data_Q_fil = filter(hd, Data_Q);

% Data_a = sqrt((Data_I_fil.*Data_I_fil)+(Data_Q_fil.*Data_Q_fil));

%% 解调
z  = hilbert(sigPhaMod);                 %希尔伯特变换对的利用---通过实部来求虚部
x1 = z.*exp(-j*2*pi*fc*t);
% x1    = Data_e.*exp(-j*2*pi*fc*t);
pha   = angle(x1);    % 求取复数信号的相位
demod = (1/100)*unwrap(pha);      % 相位校正  
figure;
plot(demod);
hold on;
plot(sigSrc);

% hd   = design(fdesign.lowpass('N,F3dB',16,100,fs),'butter');
% demod_fil = filter(hd, demod);
% figure;
% plot(demod_fil);
