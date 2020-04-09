%% 目的：实现肺部信号的解调
% By Luoshutu
clear; clc;
%% 参数设置
f0          = 3.5e6;    % 换能器中心频率 [Hz]
fs          = 40e6;     % 采样频率
prf         = 2000;     % 脉冲重复频率Hz/s

%% 加载数据
% load data_m           % 颈动脉数据

%% Read File
framenum    = 0:20;
linenum     = 2560;
channelnum  = 8;
packetnum   = 320;
pointnum    = 2048;

path     = 'F:\BEE_testSoftware\BEE8CH_Normal\left_breast_448us_8pulse_3.5MHz';
filename = '\history_data_';

Data     = zeros(linenum*length(framenum),pointnum);

for i=1:length(framenum)
    fp    = fopen([path,filename,num2str(framenum(i)),'.bin']);
    
    for j=1:packetnum
        
        tempdata  = fread(fp,pointnum*channelnum,'int16');
        data      = reshape(tempdata,channelnum,pointnum);
        
        start     = (j-1)*channelnum+1+(i-1)*linenum;
        stop      = j*channelnum+(i-1)*linenum;
        
        Data( start:stop , :) = data;
        %Data((i-1)*linenum*length(framenum)/channelnum+1:i*linenum*length(framenum)/channelnum,:)=data((),:);     
    end
    
    fclose(fp);
end

%% IQ解调
[mm,nn]  = size(Data);
Data_I = zeros(mm,nn);
Data_Q = zeros(mm,nn);
t = 0:nn-1;
for i = 1:mm
    Data_I(i,:) = cos(2*pi*f0/fs*t).*(Data(i,:));
    Data_Q(i,:) = sin(2*pi*f0/fs*t).*(Data(i,:));
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

%% 求arctan2,进行解调
sig_Arctan2 = zeros(1,mm);
for xy = 1:mm
    sig_Arctan2(xy) = atan2(sum(Data_I(xy,1000:1120)), sum(Data_Q(xy,1000:1120)));
end
% figure
% plot(sig_Arctan2);
sig_demod = unwrap(sig_Arctan2);               % 相位校正

hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,1,1.5e6,100e6),'butter'); 
sig_demod_fil = filter(hd, sig_demod);
figure;
plot(sig_demod,'b');
hold on;
plot(sig_demod_fil,'r');
legend('滤波前','滤波后');
% axis([1 mm -100 100]);

%% 画频谱
Fs       = fs;
sig      = sig_demod;
L        = length(sig);
NFFT     = 2^nextpow2(L); % Next power of 2 from length of y
Y        = fft(sig,NFFT)/L;
f        = Fs/2*linspace(0,1,NFFT/2+1);

figure;
plot(f,2*abs(Y(1:NFFT/2+1)))
title('IFilter','fontsize',14)
xlabel('Frequency (MHz)','fontsize',10)
ylabel('|Y(f)|','fontsize',10)


