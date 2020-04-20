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
framenum    = 0:5;
linenum     = 2560;
channelnum  = 8;
packetnum   = 320;
pointnum    = 2048;

path     = 'D:\_MyProject\MATLAB\right_breast_448us_8pulse_3.5MHz_2';
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
hd   = design(fdesign.lowpass('N,F3dB',16,5e6,fs),'butter');
Data_I_fil = filter(hd,Data_I.');
Data_Q_fil = filter(hd,Data_Q.');
Data_Amp = sqrt((Data_I_fil).*(Data_I_fil) + (Data_Q_fil).*(Data_Q_fil));

figure;
imagesc(20*log(1 + abs(Data_Amp)));
title('M模式');
colormap(gray);
clear Data_Amp;
%% 求arctan2,进行解调
% sound_buff = zeros(1,5000-300+1);
for qq = 9
    hd   = design(fdesign.lowpass('N,F3dB',16,qq * 1e6,fs),'butter');
    Data_I_fil = filter(hd,Data_I.');
    Data_Q_fil = filter(hd,Data_Q.');

    [mm,nn] = size(Data_I);
    sig_Arctan2 = zeros(1,mm);
    pha = zeros(1,mm);
    pha_ = zeros(1,mm);
    I_sig(1,:) = sum(Data_I_fil(300:400,:));
    Q_sig(1,:) = sum(Data_Q_fil(300:400,:));
    
    for xy = 2:mm
%         sig_Arctan2(xy) = atan2(sum(Data_I(xy,300:400)), sum(Data_Q(xy,300:400)));
        sig_Arctan2(xy) = atan2(I_sig(xy), Q_sig(xy));
        pha(xy) = atan((I_sig(xy)*Q_sig(xy-1) - I_sig(xy-1)*Q_sig(xy)) / (I_sig(xy)*I_sig(xy-1) + Q_sig(xy-1)*Q_sig(xy)));
        pha_(xy) = atan((I_sig(xy)*Q_sig(xy-1) + I_sig(xy-1)*Q_sig(xy)) / (I_sig(xy)*I_sig(xy-1) - Q_sig(xy-1)*Q_sig(xy)));
    end
    % figure
    % plot(sig_Arctan2);
    sig_demod = unwrap(sig_Arctan2);               % 相位校正
    
    hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,60,900,prf),'butter'); 
    sig_demod_fil = filter(hd, sig_demod);
    figure;
    plot(sig_demod,'b');
    hold on;
    plot(sig_demod_fil(300:end),'r');
    legend('滤波前','滤波后');
    % axis([1 mm -100 100]);
%     sound(200*sig_demod_fil,prf);
%     sound_buff(1,:) = 200*sig_demod_fil(300:5000);
end
% sound(sound_buff(6,:),prf)
% 9
pha = unwrap(pha);
pha_ = unwrap(pha_);
hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,60,900,prf),'butter'); 
phafil = filter(hd, pha);
pha_fil = filter(hd, pha_);
% sig_fL = zeros(1,mm);
% for ii = 3:mm
%     sig_fL(ii) = atan(sqrt((phafil(ii) - phafil(ii-1))/(pha_fil(ii) + pha_fil(ii-1))))/(pi*prf);
% end
% sig_fL_fil = filter(hd,sig_fL(3:end));
% sound(200*(phafil(300:end)),prf);
figure;    
plot(sig_demod_fil(300:end));
hold on;
plot(phafil(300:end));
% plot(pha_fil(300:end));
legend('sig_demod_fil','phafil');
%% 声音信号处理
% %1000:5000有效 6000:8000噪声
% hd   = design(fdesign.highpass('N,F3dB',16,80,prf),'butter'); 
% lungSound_src = filter(hd, sig_demod);
% bg_noise = lungSound_src(6000:8000-1);
% len_noi = length(bg_noise);
% lungSound_sig = lungSound_src(1000:5000-1);
% % sound(200*lungSound_sig,prf);
% %% 频谱相减
% bg_amp = abs(fft(bg_noise));
% lung_sound = [];
% a=2;b=2;%设定alpha和belta的值
% clear i;
% for pp = 0:floor(length(lungSound_sig)/len_noi)-1
%     sound_sig_seg = lungSound_sig(1+pp*len_noi:len_noi+pp*len_noi);
%     sss_fft = fft(sound_sig_seg);
%     sss_amp = abs(sss_fft);
%     sss_pha = unwrap(angle(sss_fft));
%     SSS = sss_amp.^a - bg_amp.^a;
%     SSS = abs(SSS).^(1/b);
% %     SSS = SSS.*exp(-i*sss_pha);
%     sss = ifft(SSS);
% %     thr = mean(sss);
% %     for qq = 1:len_noi
% %         if abs(sss(qq)) > thr
% %             sss(qq) = sss(qq)/4;
% %         end
% %     end
%     lung_sound = [lung_sound sss];
% end
% 
% hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,80,800,prf),'butter'); 
% lung_sound_fil = filter(hd, lung_sound);
% % figure;
% % plot(lung_sound);
% % sound(200*real(lung_sound_fil),prf)
% % audiowrite(lung_sound,prf,'NewWorld.wav');%将重新生成的增强后的语音信号存储为音频信号
% 
% %% 画频谱
% Fs       = prf;
% sig      = sig_demod_fil(300:5000);
% L        = length(sig);
% NFFT     = 2^nextpow2(L); % Next power of 2 from length of y
% Y        = fft(sig,NFFT)/L;
% f        = Fs/2*linspace(0,1,NFFT/2+1);
% 
% figure;
% plot(f,2*abs(Y(1:NFFT/2+1)))
% title('IFilter','fontsize',14)
% xlabel('Frequency (Hz)','fontsize',10)
% ylabel('|Y(f)|','fontsize',10)
