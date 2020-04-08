clear
clc

Start = 3;
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
prf         = 500;      %�����ظ�Ƶ��Hz/s

%% ��������
N            = 10000;                        %�������
z_size       = 60/1000;                      %�������[m]
MobileRange  = ceil(4e-3 / (z_size / N));    %���˶���Χ��4mm
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

%% �����ź�����
heartRate = 80/60;  %����Ƶ�ʣ�ÿ����80��
heartNum  = ceil(prf / heartRate); %һ��������ռ��������
heartSignalAmp = round([0.00 1.00 0.95 0.73 0.20 -0.20 0.06 0.08 -0.03 0.00].*(1e-4 / (z_size / N)));
heartSignal = interp(heartSignalAmp,round(heartNum/length(heartSignalAmp)));
figure;
plot(heartSignal);
grid on;

%% ��������
excitation_pulse = sin(2*pi*f0*(0:1/fs:8/f0));
excitation_pulse = excitation_pulse.*hanning(max(size(excitation_pulse)))';

%% ���
repeatNumber   = 2048;                                                  % �������
echo           = zeros(N + length(excitation_pulse) - 1,repeatNumber);  % ����������
PHANTOM        = zeros(N, repeatNumber);                                % ���建��
dam            = 5e3;                                                   % ����˥��,ֵԽС˥��Խ��
p_phantom_pha  = rand(N/2, 1) * 2 * pi;                                 % ÿ�������λ

for i = 1:repeatNumber                              %�������
    phantom = phantom_half;
    POINT_mobile = round(heartSignal(mod(i,length(heartSignal))+1));
    POINT_position = round(N/2);% - POINT_mobile); %THE POINTλ��
    %THE POINT����ķ���ɢ��ԴΪTHE POINTǰ�ľ����һ�����THE POINT�˶����˶�
    for p = 1:(POINT_position)
        phantom(POINT_position + p) = phantom(POINT_position - p + 1) / (((N-p)/dam)^3); %�����3�η�˥��
    end
    
    %THE POINT��ֵ����ɨ������仯��Ϊ20��30��֮ǰ�������ֵ֮������ֵ
    POINT_Am  = (randperm(10,1)+20)*max(phantom_half);
    phantom(POINT_position,1) = POINT_Am;
    phantom = phantom(1:N);
    
    clear j;
    if POINT_position+4999 > N
    	phantom(POINT_position:end) = phantom(POINT_position:end).*exp(j*p_phantom_pha(1:(N - POINT_position + 1)));
    else
        phantom(POINT_position:POINT_position+4999) = phantom(POINT_position:POINT_position+4999).*exp(j*p_phantom_pha);
    end
    
    PHANTOM(:,i) = phantom;         %�洢����
    
    %���
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
hd   = design(fdesign.lowpass('N,F3dB',16,5e6,100e6),'butter');
Data_I_fil = filter(hd,Data_I.');
Data_Q_fil = filter(hd,Data_Q.');

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
    slowTimeSignal_I(j,1) = sum(Data_I_fil(floor(mm/2 + MobileRange):floor(mm/2 + MobileRange*2),j));
    slowTimeSignal_Q(j,1) = sum(Data_Q_fil(floor(mm/2 + MobileRange):floor(mm/2 + MobileRange*2),j));
end

% hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,1000,6e6,100e6),'butter');
% slowTimeSignal_I = filter(hd,slowTimeSignal_I);
% slowTimeSignal_Q = filter(hd,slowTimeSignal_Q);

%% ��ʱ���ʱ����Ҷ�任
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
xlabel('ʱ�� t/s');ylabel('Ƶ�� f/Hz');
axis off;
% axis([1 y x/2 - 20 x/2 + 20]);
title(['��ʱ����ҶʱƵͼ',num2str(cir_map),'M']);

% lateralRes  = 1 / (f0/0.5e6); %unit : mm ����3��ʱ������ֱ���Ϊ1����
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
%     lateralFilter = sinc(-abs(n)/2 :0.1: abs(n)/2)/n; %����һ����Ƶ�ʼ��ٶȻ���PRF����ص���ʱ�䷽���ϵ��˲���
%  
%     slowsig = slowTimeSignal (i:i+128);
% %     slowsig = echo(5100,i:i+128);
%     slowsig = conv (slowsig, lateralFilter, 'same');
%     im(:,ind)= fftshift(abs(fft(slowsig)));
% end
% im(65,:)=0;
% figure;
% imagesc(20*log(0.1+im))
% title(['PWͼ',num2str(cir_map),'M']);
end

