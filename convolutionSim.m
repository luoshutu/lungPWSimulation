clear
clc

%% 参数设置
f0             = 4e6;                 % 换能器中心频率 [Hz]
fs             = 100e6;               % 仿真使用的采样频率
c              = 1540;                % 声速 [m/s]
lambda         = c/f0;                % 波长 [m]

%% 仿体设置
N       = 10000;
amp     = randn(N/2,1);
z_size  = 60/1000;   %  Height of phantom [m]

phantom_half = zeros(N,1);

% skin
z_skin = randperm(round(N*(30/1000)),round(N*(5/1000)));
for i = 1:z_skin
    phantom_half(i,1) = amp(i,1);
end
% phantom_half(1:round(N*(30/1000)),1) = phantom_half(1:round(N*(30/1000)),1)*hamming(round(N*(30/1000))).';

% fat
z_fat = randperm(round(N*(170/1000)),round(N*(10/1000))) + round(N*(30/1000));
for j = z_fat
	phantom_half(j,1) = amp(j,1)*2; 
end
% phantom_half(round(N*(30/1000)):round(N*(200/1000)),1) = phantom_half(round(N*(30/1000)):round(N*(200/1000)),1)*hamming(round(N*(200/1000)) - round(N*(30/1000)) + 1).';

%muscle
z_muscle = randperm(round(N*(300/1000)),round(N*(200/1000))) + round(N*(200/1000));
for k = z_muscle
	phantom_half(k,1) = amp(k,1); 
end
% phantom_half(round(N*(200/1000)):round(N*(495/1000)),1) = phantom_half(round(N*(200/1000)):round(N*(495/1000)),1)*hamming(round(N*(495/1000)) - round(N*(200/1000)) + 1).';

% figure(1);
% plot(phantom_half);
% title('Phantom');

% THE POINT
MobileRange = ceil(lambda / (z_size / N));  %点运动范围，一个波长
POINT_Am    = 100*max(phantom_half);

%% 发射脉冲
impulse_response = sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response = impulse_response.*hanning(max(size(impulse_response)))';

% figure(2);
% plot(impulse_response);
% title('ImpulseResponse');
%% 卷积
repeatNumber = 1024;
echo         = zeros(N + length(impulse_response) - 1,repeatNumber);
PHANTOM      = zeros(N,repeatNumber);
dam          = 4000;%衰减系数
mobile_temp  = ceil(MobileRange/2);
mobileRan    = -floor(MobileRange/2):1:floor(MobileRange/2);
mobileRan_temp = zeros(1,MobileRange);
mixMN = 5;
for i = 1:repeatNumber
    phantom = phantom_half;

    POINT_mobile = mobileRan_temp(mobile_temp);
%     moile(i) = POINT_mobile;
%     Line_mobile = 0;%randperm(8,1) - 4;
    for p = 1:(N/2 + POINT_mobile)
        if ((p - POINT_mobile) > 0) %&& ((p - POINT_mobile) <= (N/2 - 1))
            phantom(N - p + 1,1) = phantom(p - POINT_mobile,1)/(((N-p)/dam)^3);
        end
    end
    
    POINT_position = N/2 + POINT_mobile;
    if mod(MobileRange,2) == 0
        phantom(POINT_position-floor(MobileRange/2):POINT_position+floor(MobileRange/2),1) = POINT_Am*hamming(MobileRange + 1);
    else
        phantom(POINT_position-floor(MobileRange/2):POINT_position+floor(MobileRange/2),1) = POINT_Am*hamming(MobileRange);
    end
%     phantom(POINT_position,1) = POINT_Am;
    if mobile_temp < mixMN + 1
        mobileRan_temp = mobileRan((mobile_temp+mixMN):end);
    elseif mobile_temp > MobileRange - mixMN
        mobileRan_temp = mobileRan(1:(mobile_temp-mixMN));
    else
        mobileRan_temp = [mobileRan(1:mobile_temp-mixMN) mobileRan(mobile_temp+mixMN:end)];
    end
    mobile_temp = randperm(length(mobileRan_temp),1);
    
    phantom(:,1) = (phantom(:,1)+0.1);
    PHANTOM(:,i) = phantom(:,1);
%     hold on;
%     plot(phantom);
    echo(:,i) = conv(impulse_response,phantom.').';
end

figure;
imagesc(20*log(0.1 + abs(PHANTOM)));
title('PHANTON');
colormap(gray);
figure;
imagesc(20*log(0.1 + abs(echo)));
title('Echo');
colormap(gray);
%% IQ解调
[mm,nn]  = size(echo);
Data_I = zeros(nn,mm);
Data_Q = zeros(nn,mm);
t = 0:mm-1;

for i = 1:repeatNumber
    Data_I(i,:) = cos(2*pi*f0/fs*t).*(echo(:,i).');
    Data_Q(i,:) = sin(2*pi*f0/fs*t).*(echo(:,i).');
end

Data_I_fil = filter(Filter_IIR_L_2M,Data_I.');
Data_Q_fil = filter(Filter_IIR_L_2M,Data_Q.');

Data_Amp = sqrt((Data_I_fil).*(Data_I_fil) + (Data_Q_fil).*(Data_Q_fil));
D=20;   %  下采样抽取率
figure;
imagesc(20*log(1 + abs(Data_Amp(1:D:end,:))));
colormap(gray);
%% 累加
slowTimeSignal_I = zeros(repeatNumber,1);
slowTimeSignal_Q = zeros(repeatNumber,1);
for j = 1:repeatNumber
%     for k = 1:MobileRange
%         slowTimeSignal_I(j,1) = slowTimeSignal_I(j,1) + Data_I_fil(round(length(echo)/2 - ceil(MobileRange/2) - 1) + k,j);
%         
%         slowTimeSignal_Q(j,1) = slowTimeSignal_Q(j,1) + Data_Q_fil(round(length(echo)/2 - ceil(MobileRange/2) - 1) + k,j);
%     end
    
    slowTimeSignal_I(j,1) = sum(Data_I_fil(mm/2 - floor(MobileRange/2):mm/2 + floor(MobileRange/2),j));
    slowTimeSignal_Q(j,1) = sum(Data_Q_fil(mm/2 - floor(MobileRange/2):mm/2 + floor(MobileRange/2),j));
end

% slowTimeSignal_Amp = sqrt((slowTimeSignal_I).*(slowTimeSignal_I) + (slowTimeSignal_Q).*(slowTimeSignal_Q));
% slowTimeSignal_Amp = filter(Filter_IIR_H_100Hz,slowTimeSignal_Amp);
%% 慢时间傅里叶变换

% figure;
% plot(slowTimeSignal);
% title('slowTimeSignal');
% figure;
% spectrogram(slowTimeSignal,256,250,256,'yaxis');
clear i j;
slowTimeSignal = slowTimeSignal_I + i*slowTimeSignal_Q;
[S,F,T,P] = spectrogram(slowTimeSignal,256,250,256);
figure;       
S = fftshift(S,1);
P=20*log10(abs(S));
[x,y]=size(S); 
imagesc(P)
colorbar;
set(gca,'YDir','normal');
% colormap(gray);
xlabel('时间 t/s');
ylabel('频率 f/Hz');
title('短时傅里叶时频图');
