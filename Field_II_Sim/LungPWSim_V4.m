%% 目的：使用Field II实现肺部的PW信号仿真 
% 二维仿真，calc_h二维点扩散方程
% By Luoshutu
clear; clc;

%% 设置仿真起始
Start = 7;      % 从中心频率为Start*e6开始仿真
Num   = 1;      % 仿真
% figure;
% hold on;
for cir_map = Start:Start + Num -1
%% 参数设置
f0             = cir_map*(1e6);
% f0           = 7e6;                 % 换能器中心频率 [Hz]
fs             = 40e6;                % 仿真使用的采样频率
c              = 1540;                % 声速 [m/s]
lambda         = c/f0;                % 波长 [m]
prf            = 30;                  % 帧频Hz/s
lateralRes     = 1 / (f0/0.5e6);      % unit : mm 假设3兆时，横向分辨率为1毫米

%% 仿体设置
N            = 1000;                 % 仿体点数
z_size       = 60/1000;               % 仿体深度[m]
MobileRange  = ceil(4e-3 / (z_size / N));    %点运动范围，4mm

%% 心跳信号设置
% 肺部会有搏动症，使用M超看应与之前做的单条线的仿真一致
heartRate = 65/60;  %心跳频率，每分钟65次
heartNum  = ceil(prf / heartRate); %一次心跳所占采样线数
heartSignalAmp = round([0.00 1.00 0.95 0.73 0.20 -0.20 0.06 0.08 -0.03 0.00].*ceil(4e-4 / (z_size / N)));
heartSignal = interp(heartSignalAmp,round(heartNum/length(heartSignalAmp)));
figure;
plot(heartSignal);
grid on;

%% 仿体随时间变化为二维
lineNumber   = 1024;                              % 仿体线数
imgFrame     = 60;                                 % 一共产生几帧图像
PHANTOM      = zeros(imgFrame, N, lineNumber);    % 仿体缓存
% 胸模线
sternumLine = round(hamming(10).*ceil(5e-3 / (z_size / N)));
sternumLine = interp(sternumLine,round(2*lineNumber/length(sternumLine)));
% 呼吸波形
breathRate  = 1;  %呼吸频率,Hz
breathNum   = ceil(prf / breathRate); %一次呼吸所占帧数
breathWave  = round(hamming(10).*ceil(1e-3 / (z_size / N)));
breathWave  = interp(breathWave,round(breathNum/length(breathWave)));

% figure;
% plot(sternumLine);
% grid on;
for frame_cnt = 1:imgFrame
    % 每一帧图片，胸模线会随心跳产生的移动
    frame_mobile = round(heartSignal(mod(frame_cnt, length(heartSignal))+1)); 
    breathMove = round(breathWave(mod(frame_cnt, length(breathWave))+1));
    for i = 1:lineNumber                          % 获取二维仿体
        % 仿体线之间的相对移动加上随心跳的移动
        ster_line_position = round(sternumLine(mod(i, length(sternumLine))+lineNumber/2));
%         ster_line_position = round(sternumLine(mod(i+15*breathMove, length(sternumLine))+lineNumber/2));
        ster_line_mobile = frame_mobile + ster_line_position + breathMove;
        phantom = phantom_line_2(N, ster_line_mobile, i, 20*breathMove);
        PHANTOM(frame_cnt, :, i) = phantom;       % 存储仿体
    end
end
%% 仿体显示
figure;
filename='D:\\_MyProject\\MATLAB\\lungPWSimulation\\Field_II_Sim\\Picture\\phantom9.gif';          %输出路径+保存的文件名.gif
% while(1)
for frame_cnt = 1:imgFrame
%     figure(frame_cnt);
    imgPhantom(:,:) = PHANTOM(frame_cnt,:,:);
    image(20*(1+abs(imgPhantom)));
    colormap(gray);
    title('PHANTOM');
    axis([1 lineNumber 1 N*(650/1000)]);
    pause(0.03);
    set(gcf,'color','w');  %设置背景为白色
    set(gca,'units','pixels','Visible','off');
    frame = getframe(gcf); 
    im = frame2im(frame);     %将影片动画转换为编址图像,因为图像必须是index索引图像
    imshow(im);
    [I,map] = rgb2ind(im,20); %将真彩色图像转化为索引图像
    if frame_cnt==1
        imwrite(I,map,filename,'gif','Loopcount',inf,'DelayTime',0.1);     %Loopcount只是在i==1的时候才有用
    else
        imwrite(I,map,filename,'gif','WriteMode','append','DelayTime',0.1);%DelayTime:帧与帧之间的时间间隔
    end
end
% end
close all;
%% 画声场图
% fp = cal_field(N,z_size,lineNumber,fs,f0);

figure;
imagesc(20*log(1 + abs(fp)));
title('声场');
% colormap(gray);

%% 卷积得回波
echoBuff = zeros(imgFrame, N, lineNumber);    % 缓存
timeGain = ones(1,N);                         % 时间增益曲线

for idx_p = 1:N
    timeGain(1,idx_p) = (1e10)*((idx_p/(N*(650/1000)/2) + 1)^2);
end

for frame_cnt = 1:imgFrame
    singleFrame(:,:) = PHANTOM(frame_cnt,:,:);
    echoTemp = abs(conv2(fp, singleFrame));
    echoBuff(frame_cnt,:,:) = echoTemp(1:N,1:lineNumber).*timeGain';
end
%% Echo显示
figure;
filename_E='D:\\_MyProject\\MATLAB\\lungPWSimulation\\Field_II_Sim\\Picture\\echo9_c.gif';          %输出路径+保存的文件名.gif
% while(1)
for frame_cnt = 1:imgFrame
%     figure(frame_cnt);
    imgEcho(:,:) = echoBuff(frame_cnt,:,:);
    imagesc(20*log(10 + real(imgEcho)));
%     colormap(gray);
    title('卷积回波');
    axis([1 lineNumber 1 N*(650/1000)]);
    pause(0.03);
    set(gcf,'color','w');  %设置背景为白色
    set(gca,'units','pixels','Visible','off');
    frame = getframe(gcf); 
    im = frame2im(frame);     %将影片动画转换为编址图像,因为图像必须是index索引图像
    imshow(im);
    [I,map] = rgb2ind(im,20); %将真彩色图像转化为索引图像
    if frame_cnt==1
        imwrite(I,map,filename_E,'gif','Loopcount',inf,'DelayTime',0.1);     %Loopcount只是在i==1的时候才有用
    else
        imwrite(I,map,filename_E,'gif','WriteMode','append','DelayTime',0.1);%DelayTime:帧与帧之间的时间间隔
    end
end
% end
close all;

%% IQ解调
[fra,mm,nn]  = size(echoBuff);
Data_I = zeros(nn,mm);
Data_Q = zeros(nn,mm);
dataIBuff = zeros(fra,mm,nn);
dataQBuff = zeros(fra,mm,nn);

hd   = design(fdesign.lowpass('N,F3dB',16,10e6,fs),'butter');
for idx = 1:fra
    t = 0:mm-1;
    for i = 1:nn
        IQTemp = echoBuff(idx,:,i);
        Data_I(i,:) = cos(2*pi*f0/fs*t).*(IQTemp);
        Data_Q(i,:) = sin(2*pi*f0/fs*t).*(IQTemp);
    end

    % 滤波
    Data_I_fil = filter(hd,Data_I.');
    Data_Q_fil = filter(hd,Data_Q.');
    
    dataIBuff(idx,:,:) = Data_I_fil;
    dataQBuff(idx,:,:) = Data_Q_fil;
end
%% 求模
% Data_Amp = sqrt((Data_I).*(Data_I) + (Data_Q).*(Data_Q));
Data_Amp = sqrt((dataIBuff).*(dataIBuff) + (dataQBuff).*(dataQBuff));

%% 显示M模式
D=2;   %  下采样抽取率
clear imgBModel
filename_B='D:\\_MyProject\\MATLAB\\lungPWSimulation\\Field_II_Sim\\Picture\\BModel9_c.gif';
figure;
% while(1)
for frame_cnt = 1:imgFrame
%     figure;
    imgBModel(:,:) = Data_Amp(frame_cnt,1:end,:);
    imagesc(20*log(1 + abs(imgBModel)));
%     colormap(gray);
    title(['B模式',num2str(cir_map),'M']);
    axis([1 lineNumber 1 N*(650/1000)]);
    pause(0.03);
    set(gcf,'color','w');  %设置背景为白色
    set(gca,'units','pixels','Visible','off');
    frame = getframe(gcf); 
    im = frame2im(frame);     %将影片动画转换为编址图像,因为图像必须是index索引图像
    imshow(im);
    [I,map] = rgb2ind(im,20); %将真彩色图像转化为索引图像
    if frame_cnt==1
        imwrite(I,map,filename_B,'gif','Loopcount',inf,'DelayTime',0.1);     %Loopcount只是在i==1的时候才有用
    else
        imwrite(I,map,filename_B,'gif','WriteMode','append','DelayTime',0.1);%DelayTime:帧与帧之间的时间间隔
    end
end
% end
close all;
%% 累加
% slowTimeSignal_I = zeros(lineNumber,1);
% slowTimeSignal_Q = zeros(lineNumber,1);
% 
% for j = 1:lineNumber
%     slowTimeSignal_I(j,1) = sum(Data_I_fil(floor(mm/2 + MobileRange):floor(mm/2 + MobileRange*4),j));
%     slowTimeSignal_Q(j,1) = sum(Data_Q_fil(floor(mm/2 + MobileRange):floor(mm/2 + MobileRange*4),j));
% end
% 
% % hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,1000,6e6,100e6),'butter');
% % slowTimeSignal_I_fil = filter(hd,slowTimeSignal_I);
% % slowTimeSignal_Q_fil = filter(hd,slowTimeSignal_Q);
% 
% %% 慢时间短时傅里叶变换
% clear i;
% slowTimeSignal = slowTimeSignal_I + i*slowTimeSignal_Q;
% % hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,10,50e6,100e6),'butter');
% % slowTimeSignal = filter(hd,slowTimeSignal);
% 
% [S,F,T,~] = spectrogram(slowTimeSignal,256,250,256);
% figure;       
% S = fftshift(S,1);
% P = 20*log(1 + abs(S));
% [x,y] = size(P);
% % subplot(Num,2,2*(cir_map-Start)+2);
% imagesc(P); colorbar;
% set(gca,'YDir','normal');
% % colormap(gray);
% xlabel('时间 t/s');ylabel('频率 f/Hz');
% % axis off;
% % axis([1 y x/2 - 20 x/2 + 20]);
% title(['短时傅里叶时频图',num2str(cir_map),'M']);
end


