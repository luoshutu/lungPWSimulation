%% 目的：仿真实现肺部信号的解调
% By Luoshutu
clear; clc;
%% 设置起始信息
Start = 3;
Num   = 1;
% figure;
% hold on;
for cir_map = Start:Start + Num -1
%% 参数设置
f0          = cir_map*(1e6);
% f0          = 7e6;      % 换能器中心频率 [Hz]
fs          = 100e6;	% 仿真使用的采样频率
c           = 1540;     % 声速 [m/s]
lambda      = c/f0;     % 波长 [m]
prf         = 500;      %脉冲重复频率Hz/s

%% 仿体设置
N            = 10000;                        %仿体点数
z_size       = 60/1000;                      %仿体深度[m]
MobileRange  = ceil(4e-3 / (z_size / N));    %点运动范围，4mm
phantom_half = zeros(N,1);    
amp          = randn(N/2,1);

% skin
z_skin = randperm(round(N*(30/1000)),round(N*(5/1000)));
for i = 1:z_skin
    phantom_half(i,1) = amp(i,1);
end

% fat
z_fat = randperm(round(N*(170/1000)),round(N*(10/1000))) + round(N*(30/1000));
for j = z_fat
	phantom_half(j,1) = amp(j,1)*2; 
end

%muscle
z_muscle = randperm(round(N*(300/1000)),round(N*(200/1000))) + round(N*(200/1000));
for k = z_muscle
	phantom_half(k,1) = amp(k,1); 
end

%% 心跳信号设置
heartRate = 80/60;  %心跳频率，每分钟80次
heartNum  = ceil(prf / heartRate); %一次心跳所占采样线数
heartSignalAmp = round([0.00 1.00 0.95 0.73 0.20 -0.20 0.06 0.08 -0.03 0.00 0.00 0.00].*(1e-4 / (z_size / N)));
heartSignal = roundn(interp(heartSignalAmp,round(heartNum/length(heartSignalAmp))),-2);
figure;
plot(heartSignal);
grid on;

%% humming信号(哼哼)
% f_hum   = 80;                               %哼哼频率，单位Hz
% v_hum   = 3e-5;                             %哼哼最大幅值，单位m
% % humming = (v_hum/(z_size/N)) * sin(2*pi*f_hum*(0:1/prf:3/f_hum));
% % humming = (v_hum/(z_size/N)) * sawtooth(2*pi*f_hum*(0:1/prf:3/f_hum),0.4);
% humming = (v_hum/(z_size/N)) * square(2*pi*f_hum*(0:1/prf:3/f_hum),50);        %占空比为50%的方波
% % humming = (v_hum/(z_size/N)) * (2 * rand(1,round(prf/f_hum)) - 1);
% figure;
% plot(humming);
% title('Humming');
% grid on;
%% 激励函数
excitation_pulse = sin(2*pi*f0*(0:1/fs:8/f0));
excitation_pulse = excitation_pulse.*hanning(max(size(excitation_pulse)))';

%% 卷积
repeatNumber   = 1024;                                                  % 发射次数
echo           = zeros(N + length(excitation_pulse) - 1,repeatNumber);  % 卷积结果缓存
PHANTOM        = zeros(N, repeatNumber);                                % 仿体缓存
dam            = 5e3;                                                   % 仿体衰减,值越小衰减越大
p_phantom_pha  = rand(N/2, 1) * 2 * pi;                                 % 每个点的相位

for i = 1:repeatNumber                              %卷积次数
    phantom = phantom_half;
    POINT_mobile = round(heartSignal(mod(i,length(heartSignal))+1));% + round(humming(mod(i,length(humming)) + 1));
    POINT_position = round(N/2 - POINT_mobile); %THE POINT位置
    %THE POINT后面的仿体散射源为THE POINT前的镜像，且会随着THE POINT运动而运动
    for p = 1:(POINT_position)
        phantom(POINT_position + p) = phantom(POINT_position - p + 1) / (((N-p)/dam)^3); %随距离3次方衰减
    end
    
    %THE POINT幅值，随扫描次数变化，为20到30倍之前仿体最大值之间的随机值
    POINT_Am  = (randperm(10,1)+20)*max(phantom_half);
    phantom(POINT_position,1) = POINT_Am;
    phantom = phantom(1:N);
    
%     clear j;
%     if POINT_position+4999 > N
%     	phantom(POINT_position:end) = phantom(POINT_position:end).*exp(j*p_phantom_pha(1:(N - POINT_position + 1)));
%     else
%         phantom(POINT_position:POINT_position+4999) = phantom(POINT_position:POINT_position+4999).*exp(j*p_phantom_pha);
%     end
    
    PHANTOM(:,i) = phantom;         %存储仿体
    
    %卷积
    echo(:,i) = conv(excitation_pulse,phantom.').';
end

% figure;
% imagesc(20*log(0.1 + abs(PHANTOM)));
% title('PHANTON');
% colormap(gray);
% figure;
% imagesc(20*log(0.1 + abs(echo)));
% title('Echo');
% colormap(gray);
%% IQ解调
[mm,nn]  = size(echo);
Data_I = zeros(nn,mm);
Data_Q = zeros(nn,mm);
t = 0:mm-1;

for i = 1:repeatNumber
    Data_I(i,:) = cos(2*pi*f0/fs*t).*(echo(:,i).');
    Data_Q(i,:) = sin(2*pi*f0/fs*t).*(echo(:,i).');
end

%% 显示M模式
hd   = design(fdesign.lowpass('N,F3dB',16,5e6,100e6),'butter');
Data_I_fil = filter(hd,Data_I.');
Data_Q_fil = filter(hd,Data_Q.');
Data_Amp = sqrt((Data_I_fil).*(Data_I_fil) + (Data_Q_fil).*(Data_Q_fil));
D=20;   %  下采样抽取率
figure;
% subplot(Num,2,2*(cir_map-Start)+1);
imagesc(20*log(1 + abs(Data_Amp(1:D:end,:))));
title(['M模式',num2str(cir_map),'M']);
colormap(gray);

%% 相位解调
% xx = 0:(1/prf):(1/prf)*(repeatNumber-1);
% % z  = hilbert(sum(real(echo(6000:6020,:))));                 %希尔伯特变换对的利用---通过实部来求虚部
% z = sum(echo(6000:6020,:));
% x1 = z.*exp(-j*2*pi*prf*xx);
% % x1    = Data_e.*exp(-j*2*pi*fc*t);
% pha   = angle(x1);          % 求取复数信号的相位
% demod = unwrap(pha);        % 相位校正  
% hd   = design(fdesign.lowpass('N,F3dB',16,1.5e6,100e6),'butter');
% demod_fil = filter(hd,demod);
% figure;
% plot(demod_fil);
%% 求arctan2
for xy = 1:repeatNumber
    De(xy) = atan2(sum((Data_I(xy,6000:6020))),sum((Data_Q(xy,6000:6020))));
end
Dt = unwrap(De);
hd   = design(fdesign.lowpass('N,F3dB',16,1.5e6,100e6),'butter');
Dt_fil = filter(hd,Dt);
figure
plot(Dt_fil);
% axis([0 400 -3 4]);
end

