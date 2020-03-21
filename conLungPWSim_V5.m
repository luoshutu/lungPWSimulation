clear
clc

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
heartSignalAmp = round([0.00 1.00 0.95 0.73 0.20 -0.20 0.06 0.08 -0.03 0.00].*(1e-4 / (z_size / N)));
heartSignal = interp(heartSignalAmp,round(heartNum/length(heartSignalAmp)));
figure;
plot(heartSignal);
grid on;

%% 激励函数
excitation_pulse = sin(2*pi*f0*(0:1/fs:8/f0));
excitation_pulse = excitation_pulse.*hanning(max(size(excitation_pulse)))';

%% 卷积
repeatNumber   = 2048;                                                  % 发射次数
echo           = zeros(N + length(excitation_pulse) - 1,repeatNumber);  % 卷积结果缓存
PHANTOM        = zeros(N, repeatNumber);                                % 仿体缓存
dam            = 5e3;                                                   % 仿体衰减,值越小衰减越大
p_phantom_pha  = rand(N/2, 1) * 2 * pi;                                 % 每个点的相位

for i = 1:repeatNumber                              %卷积次数
    phantom = phantom_half;
    POINT_mobile = round(heartSignal(mod(i,length(heartSignal))+1));
    POINT_position = round(N/2);% - POINT_mobile); %THE POINT位置
    %THE POINT后面的仿体散射源为THE POINT前的镜像，且会随着THE POINT运动而运动
    for p = 1:(POINT_position)
        phantom(POINT_position + p) = phantom(POINT_position - p + 1) / (((N-p)/dam)^3); %随距离3次方衰减
    end
    
    %THE POINT幅值，随扫描次数变化，为20到30倍之前仿体最大值之间的随机值
    POINT_Am  = (randperm(10,1)+20)*max(phantom_half);
    phantom(POINT_position,1) = POINT_Am;
    phantom = phantom(1:N);
    
    clear j;
    if POINT_position+4999 > N
    	phantom(POINT_position:end) = phantom(POINT_position:end).*exp(j*p_phantom_pha(1:(N - POINT_position + 1)));
    else
        phantom(POINT_position:POINT_position+4999) = phantom(POINT_position:POINT_position+4999).*exp(j*p_phantom_pha);
    end
    
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

% 滤波
hd   = design(fdesign.lowpass('N,F3dB',16,5e6,100e6),'butter');
Data_I_fil = filter(hd,Data_I.');
Data_Q_fil = filter(hd,Data_Q.');

%% 显示M模式
Data_Amp = sqrt((Data_I_fil).*(Data_I_fil) + (Data_Q_fil).*(Data_Q_fil));
D=20;   %  下采样抽取率
figure;
% subplot(Num,2,2*(cir_map-Start)+1);
imagesc(20*log(1 + abs(Data_Amp(1:D:end,:))));
title(['M模式',num2str(cir_map),'M']);
colormap(gray);

%% 累加
slowTimeSignal_I = zeros(repeatNumber,1);
slowTimeSignal_Q = zeros(repeatNumber,1);

for j = 1:repeatNumber
    slowTimeSignal_I(j,1) = sum(Data_I_fil(floor(mm/2 + MobileRange):floor(mm/2 + MobileRange*2),j));
    slowTimeSignal_Q(j,1) = sum(Data_Q_fil(floor(mm/2 + MobileRange):floor(mm/2 + MobileRange*2),j));
end

% hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,1000,6e6,100e6),'butter');
% slowTimeSignal_I = filter(hd,slowTimeSignal_I);
% slowTimeSignal_Q = filter(hd,slowTimeSignal_Q);

%% 慢时间短时傅里叶变换
clear i;
slowTimeSignal = slowTimeSignal_I + i*slowTimeSignal_Q;
% hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,10000,5e5,1e6),'butter');
% slowTimeSignal = filter(hd,slowTimeSignal);

[S,F,T,~] = spectrogram(slowTimeSignal,256,250,256);
figure;       
S = fftshift(S,1);
P = 20*log(100 + abs(S));
[x,y]=size(P); 
% subplot(Num,2,2*(cir_map-Start)+2);
imagesc(P);
colorbar;
set(gca,'YDir','normal');
% colormap(gray);
xlabel('时间 t/s');ylabel('频率 f/Hz');
axis off;
% axis([1 y x/2 - 20 x/2 + 20]);
title(['短时傅里叶时频图',num2str(cir_map),'M']);

% lateralRes  = 1 / (f0/0.5e6); %unit : mm 假设3兆时，横向分辨率为1毫米
% step = 20;
% ind  = 0;
% for i = 1: step: repeatNumber - 128
%     ind = ind+1;
%     v = heartSignal(mod(i,length(heartSignal))+1);
% 
%     distancePerPulse = v/prf;  % v unit mm/s
%     n = lateralRes/distancePerPulse;
%     n = floor (n); 
%     nn(i) = n;
%     lateralFilter = sinc(-abs(n)/2 :0.1: abs(n)/2)/n; %构建一个和频率及速度还有PRF都相关的慢时间方向上的滤波器
%  
%     slowsig = slowTimeSignal (i:i+128);
% %     slowsig = echo(5100,i:i+128);
%     slowsig = conv (slowsig, lateralFilter, 'same');
%     im(:,ind)= fftshift(abs(fft(slowsig)));
% end
% im(65,:)=0;
% figure;
% imagesc(20*log(0.1+im))
% title(['PW图',num2str(cir_map),'M']);
end

