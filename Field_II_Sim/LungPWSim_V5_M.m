%% Ŀ�ģ�ʹ��Field IIʵ�ַβ���PW�źŷ��� 
% ��ά���棬calc_h��ά����ɢ����
% By Luoshutu
clear; clc;

%% ���÷�����ʼ
Start = 7;      % ������Ƶ��ΪStart*e6��ʼ����
Num   = 1;      % ����
% figure;
% hold on;
for cir_map = Start:Start + Num -1
%% ��������
f0             = cir_map*(1e6);
% f0           = 7e6;                 % ����������Ƶ�� [Hz]
fs             = 40e6;                % ����ʹ�õĲ���Ƶ��
c              = 1540;                % ���� [m/s]
lambda         = c/f0;                % ���� [m]
prf            = 500;                 % ֡ƵHz/s
lateralRes     = 1 / (f0/0.5e6);      % unit : mm ����3��ʱ������ֱ���Ϊ1����

%% ��������
N            = 1000;                 % �������
z_size       = 60/1000;               % �������[m]
MobileRange  = ceil(4e-3 / (z_size / N));    %���˶���Χ��4mm

%% �����ź�����
% �β����в���֢��ʹ��M����Ӧ��֮ǰ���ĵ����ߵķ���һ��
heartRate = 65/60;  %����Ƶ�ʣ�ÿ����65��
heartNum  = ceil(prf / heartRate); %һ��������ռ��������
heartSignalAmp = round([0.00 1.00 0.95 0.73 0.20 -0.20 0.06 0.08 -0.03 0.00].*ceil(4e-4 / (z_size / N)));
heartSignal = interp(heartSignalAmp,round(heartNum/length(heartSignalAmp)));
figure;
plot(heartSignal);
grid on;

%% ������ʱ��仯Ϊ��ά
lineNumber   = 15;                                 % ��������
imgFrame     = 2048;                              % һ��������֡ͼ��
PHANTOM      = zeros(imgFrame, N, lineNumber);    % ���建��
% ��ģ��
sternumLine = round(hamming(10).*ceil(5e-3 / (z_size / N)));
sternumLine = interp(sternumLine,round(2*lineNumber/length(sternumLine)));
% ��������
breathRate  = 1;  %����Ƶ��,Hz
breathNum   = ceil(prf / breathRate); %һ�κ�����ռ֡��
breathWave  = round(hamming(10).*ceil(1e-3 / (z_size / N)));
breathWave  = interp(breathWave,round(breathNum/length(breathWave)));

% figure;
% plot(sternumLine);
% grid on;
for frame_cnt = 1:imgFrame
    % ÿһ֡ͼƬ����ģ�߻��������������ƶ�
    frame_mobile = round(heartSignal(mod(frame_cnt, length(heartSignal))+1)); 
    breathMove = round(breathWave(mod(frame_cnt, length(breathWave))+1));
    for i = 1:lineNumber                          % ��ȡ��ά����
        % ������֮�������ƶ��������������ƶ�
        ster_line_position = round(sternumLine(mod(i, length(sternumLine))+round(lineNumber/2)));
%         ster_line_position = round(sternumLine(mod(i+15*breathMove, length(sternumLine))+lineNumber/2));
        ster_line_mobile = frame_mobile + ster_line_position + breathMove;
        phantom = phantom_line_2(N, ster_line_mobile, i, abs(20*breathMove));
        PHANTOM(frame_cnt, :, i) = phantom;       % �洢����
    end
end
%% ������ʾ
figure;
clear imgPhantom;
imgPhantom(:,:) = PHANTOM(:,:,1)';
image(20*(1+abs(imgPhantom)));
colormap(gray);
title('PHANTOM');
set(gca,'units','pixels','Visible','off');
axis([1 imgFrame 1 N*(650/1000)]);

%% ������ͼ
% fp = cal_field(N,z_size,1024,fs,f0);

figure;
imagesc(20*log(1 + abs(fp)));
title('����');
% colormap(gray);

%% ����ûز�
echoBuff = zeros(imgFrame, N, lineNumber);    % ����
timeGain = ones(1,N);                         % ʱ����������

for idx_p = 1:N
    timeGain(1,idx_p) = (1e10)*((idx_p/(N*(650/1000)/2) + 1)^2);
end

for frame_cnt = 1:imgFrame
    singleFrame(:,:) = PHANTOM(frame_cnt,:,:);
    echoTemp = abs(conv2(fp, singleFrame));
    echoBuff(frame_cnt,:,:) = echoTemp(1:N,1:lineNumber).*timeGain';
end
%% Echo��ʾ
figure;
clear imgEcho;
imgEcho(:,:) = echoBuff(:,:,1)';
imagesc(20*log(10 + abs(imgEcho)));
colormap(gray);
title('����ز�');
set(gca,'units','pixels','Visible','off');
axis([1 imgFrame 1 N*(650/1000)]);

%% IQ���
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

    % �˲�
    Data_I_fil = filter(hd,Data_I.');
    Data_Q_fil = filter(hd,Data_Q.');
    
    dataIBuff(idx,:,:) = Data_I_fil;
    dataQBuff(idx,:,:) = Data_Q_fil;
end
%% ��ģ
% Data_Amp = sqrt((Data_I).*(Data_I) + (Data_Q).*(Data_Q));
Data_Amp = sqrt((dataIBuff).*(dataIBuff) + (dataQBuff).*(dataQBuff));

%% ��ʾMģʽ
D=2;   %  �²�����ȡ��
clear imgBModel
figure;
imgBModel(:,:) = Data_Amp(:,1:end,1)';
imagesc(20*log(1 + abs(imgBModel)));
colormap(gray);
title(['Bģʽ',num2str(cir_map),'M']);
axis([1 imgFrame 1 N*(650/1000)]);
set(gca,'units','pixels','Visible','off');
%% �ۼ�
slowTimeSignal_I = zeros(imgFrame,1);
slowTimeSignal_Q = zeros(imgFrame,1);

for idx_j = 1:imgFrame
    slowTimeSignal_I(idx_j,1) = sum(dataIBuff(idx_j,floor(N/2 + 50):floor(N/2 + 50*2),1));
    slowTimeSignal_Q(idx_j,1) = sum(dataQBuff(idx_j,floor(N/2 + 50):floor(N/2 + 50*2),1));
end

% hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,1000,6e6,100e6),'butter');
% slowTimeSignal_I_fil = filter(hd,slowTimeSignal_I);
% slowTimeSignal_Q_fil = filter(hd,slowTimeSignal_Q);

%% ��ʱ���ʱ����Ҷ�任
clear i;
slowTimeSignal = slowTimeSignal_I + i*slowTimeSignal_Q;
% hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,10,50e6,100e6),'butter');
% slowTimeSignal = filter(hd,slowTimeSignal);

[S,F,T,~] = spectrogram(slowTimeSignal,256,250,256);
figure;       
S = fftshift(S,1);
P = 20*log(100 + abs(S));
[x,y] = size(P);
% subplot(Num,2,2*(cir_map-Start)+2);
imagesc(P); colorbar;
set(gca,'YDir','normal');
colormap(gray);
xlabel('ʱ�� t/s');ylabel('Ƶ�� f/Hz');
% axis off;
% axis([1 y x/2 - 20 x/2 + 20]);
title(['��ʱ����ҶʱƵͼ',num2str(cir_map),'M']);
end


