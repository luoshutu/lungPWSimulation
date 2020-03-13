clear
clc

Start = 8;
Num   = 1;
% figure;
% hold on;
for cir_map = Start:Start + Num -1
%% ��������
f0          = cir_map*(1e6);
% f0          = 7e6;      % ����������Ƶ�� [Hz]
fs          = 100e6;	% ����ʹ�õĲ���Ƶ��
c           = 1540;     % ���� [m/s]
lambda      = c/f0;     % ���� [m]

%% ��������
N            = 10000;                        %�������
z_size       = 60/1000;                      %�������[m]
MobileRange  = ceil(5e-3 / (z_size / N));    %���˶���Χ��3mm
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

% figure(1);
% plot(phantom_half);
% title('phantom_half');

%% ��������
excitation_pulse = sin(2*pi*f0*(0:1/fs:8/f0));
excitation_pulse = excitation_pulse.*hanning(max(size(excitation_pulse)))';

% figure(2);
% plot(excitation_pulse);
% title('Excitation Pulse');
%% ���
repeatNumber   = 1024;                                                  %�������
echo           = zeros(N + length(excitation_pulse) - 1,repeatNumber);  %����������
PHANTOM        = zeros(N,repeatNumber);                                 %���建��
dam            = 5000;                                                  %����˥��,ֵԽС˥��Խ��
mobileRan      = -floor(MobileRange/2):1:floor(MobileRange/2);          %������λ���ƶ�ֵ����
mobile_temp    = ceil(MobileRange/2);                                   %���ȡֵ
mixMN          = ceil(1e-3 / (z_size / N));                             %��ÿ���˶�����һ�ε�λ�����ټ������
mobileRan_temp = zeros(1,MobileRange);                                  %ȥ����һ�����ĵ�λ������mixMN�������ƶ�ֵ����

for i = 1:repeatNumber                              %�������
    phantom = phantom_half;
    POINT_mobile = mobileRan_temp(mobile_temp);     %����THE POINT�˶�����
    POINT_position = N/2 + POINT_mobile;            %THE POINTλ��
%     moile(i) = POINT_mobile;
    Line_mobile = randperm(64,1) - 32;
    %THE POINT����ķ���ɢ��ԴΪTHE POINTǰ�ľ����һ�����THE POINT�˶����˶�
    for p = 1:(POINT_position - floor(MobileRange/2))
        phantom(POINT_position + floor(MobileRange/2) + p) = phantom(POINT_position - floor(MobileRange/2) - p + 1) / (((N-p)/dam)^3); %�����3�η�˥��
    end
    
    %THE POINT��ֵ����ɨ������仯��Ϊ20��30��֮ǰ�������ֵ֮������ֵ
    POINT_Am  = (randperm(10,1)+20)*max(phantom_half);  
    %
    if mod(MobileRange,2) == 0
        phantom(N/2 + Line_mobile - 32 : N/2 + Line_mobile + 32,1) = POINT_Am;%.*hamming(MobileRange + 1);
    else
        phantom(N/2 + Line_mobile - 32 : N/2 + Line_mobile + 32,1) = POINT_Am;%.*hamming(MobileRange);
    end
%     phantom(POINT_position,1) = POINT_Am;
    %ÿ���ƶ����ƶ���Χ���ƣ����ƶ�ֵ������ȥ����һ������λ������mixMN����
    if mobile_temp < mixMN + 1
        mobileRan_temp = mobileRan((mobile_temp+mixMN):end);
    elseif mobile_temp > MobileRange - mixMN
        mobileRan_temp = mobileRan(1:(mobile_temp-mixMN));
    else
        mobileRan_temp = [mobileRan(1:mobile_temp-mixMN) mobileRan(mobile_temp+mixMN:end)];
    end
    mobile_temp = randperm(length(mobileRan_temp),1);
    
    phantom = (phantom(1:N,1)+0.1);  
    PHANTOM(:,i) = phantom;         %�洢����
    
    %���
    echo(:,i) = conv(excitation_pulse,phantom.').';
end

figure;
imagesc(20*log(0.1 + abs(PHANTOM)));
title('PHANTON');
colormap(gray);
figure;
imagesc(20*log(0.1 + abs(echo)));
title('Echo');
colormap(gray);
%% IQ���
[mm,nn]  = size(echo);
Data_I = zeros(nn,mm);
Data_Q = zeros(nn,mm);
t = 0:mm-1;

for i = 1:repeatNumber
    Data_I(i,:) = cos(2*pi*f0/fs*t).*(echo(:,i).');
    Data_Q(i,:) = sin(2*pi*f0/fs*t).*(echo(:,i).');
end

% �˲�
Data_I_fil = filter(Filter_IIR_L_2M,Data_I.');
Data_Q_fil = filter(Filter_IIR_L_2M,Data_Q.');

%% ��ʾMģʽ
Data_Amp = sqrt((Data_I_fil).*(Data_I_fil) + (Data_Q_fil).*(Data_Q_fil));
D=20;   %  �²�����ȡ��
figure;
% subplot(Num,2,2*(cir_map-Start)+1);
imagesc(20*log(1 + abs(Data_Amp(1:D:end,:))));
title(['Mģʽ',num2str(cir_map),'M']);
colormap(gray);

%% �ۼ�
slowTimeSignal_I = zeros(repeatNumber,1);
slowTimeSignal_Q = zeros(repeatNumber,1);

for j = 1:repeatNumber
    slowTimeSignal_I(j,1) = sum(Data_I_fil(floor(mm/2 + MobileRange):floor(mm/2 + MobileRange*3),j));
    slowTimeSignal_Q(j,1) = sum(Data_Q_fil(floor(mm/2 + MobileRange):floor(mm/2 + MobileRange*3),j));
end
%% ��ʱ���ʱ����Ҷ�任
clear i;
slowTimeSignal = slowTimeSignal_I + i*slowTimeSignal_Q;
[S,F,T,~] = spectrogram(slowTimeSignal,256,250,256);

figure;       
S = fftshift(S,1);
P = 20*log10(1 + abs(S));
[x,y]=size(S); 
% subplot(Num,2,2*(cir_map-Start)+2);
imagesc(P);
colorbar;
set(gca,'YDir','normal');
% colormap(gray);
xlabel('ʱ�� t/s');ylabel('Ƶ�� f/Hz');
title(['��ʱ����ҶʱƵͼ',num2str(cir_map),'M']);
end

