%% Ŀ�ģ�ʹ��Field IIʵ�ַβ���PW�źŷ��� 
% ��ά���棬calc_h��ά����ɢ����
% By Luoshutu
clear; clc;
%% ��ʼ��
path(path,'D:/MATLAB/Field_II_ver_3_24_windows_gcc');
field_init;

%% ���÷�����ʼ
Start = 3;      % ������Ƶ��ΪStart*e6��ʼ����
Num   = 1;      % ����
% figure;
% hold on;
for cir_map = Start:Start + Num -1
%% ��������
f0             = cir_map*(1e6);
% f0           = 7e6;                 % ����������Ƶ�� [Hz]
fs             = 100e6;               % ����ʹ�õĲ���Ƶ��
c              = 1540;                % ���� [m/s]
lambda         = c/f0;                % ���� [m]
prf            = 500;                 % �����ظ�Ƶ��Hz/s
lateralRes     =  1 / (f0/0.5e6);     % unit : mm ����3��ʱ������ֱ���Ϊ1����

width          = 0.3/1000;            % ��Ԫ��� [m]
element_height = 5/1000;              % ��Ԫ�߶� [m]
kerf           = 0.1/1000;            % ��Ԫ��϶��� [m]
Ne             = 128;                 % ��Ԫ����
focus          = [0 0 70]/1000;       % �̶�����λ�� [m] 
focal_depth    = focus(3);            % �������
array_size     = (kerf+width)*Ne;     % ��Ԫ�ܿ��

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
heartSignalAmp = round([0.00 1.00 0.95 0.73 0.20 -0.20 0.06 0.08 -0.03 0.00].*ceil(1e-3 / (z_size / N)));
heartSignal = interp(heartSignalAmp,round(heartNum/length(heartSignalAmp)));
figure;
plot(heartSignal);
grid on;
%% ������ʱ��仯Ϊ��ά
repeatNumber   = 1024;                                                  %�������
PHANTOM        = zeros(N,repeatNumber);                                 %���建��
dam            = 5000;                                                  %����˥��,ֵԽС˥��Խ��
mobileRan      = -floor(MobileRange/2):1:floor(MobileRange/2);          %������λ���ƶ�ֵ����
mobile_temp    = ceil(MobileRange/2);                                   %���ȡֵ
mixMN          = ceil(1.5e-3 / (z_size / N));                           %��ÿ���˶�����һ�ε�λ�����ټ������
mobileRan_temp = zeros(1,MobileRange);                                  %ȥ����һ�����ĵ�λ������mixMN�������ƶ�ֵ����

for i = 1:repeatNumber                              %�������
    phantom = phantom_half;
    POINT_mobile = mobileRan_temp(mobile_temp);     %����THE POINT�˶�����
    POINT_position = round(N/2 + POINT_mobile + heartSignal(mod(i,length(heartSignal))+1));            %THE POINTλ��
    Line_mobile = randperm(64,1) - 32;
    %THE POINT����ķ���ɢ��ԴΪTHE POINTǰ�ľ����һ�����THE POINT�˶����˶�
    for p = 1:(POINT_position - floor(MobileRange/2))
        phantom(POINT_position + floor(MobileRange/2) + p) = phantom(POINT_position - floor(MobileRange/2) - p + 1) / (((N-p)/dam)^3); %�����3�η�˥��
    end
    
    %THE POINT��ֵ����ɨ������仯��Ϊ20��30��֮ǰ�������ֵ֮������ֵ
    POINT_Am  = (randperm(10,1)+20)*max(phantom_half);  
    phantom(N/2 + Line_mobile - 3 : N/2 + Line_mobile + 3,1) = POINT_Am;

    %ÿ���ƶ����ƶ���Χ���ƣ����ƶ�ֵ������ȥ����һ������λ������mixMN����
    if mobile_temp < mixMN + 1
        mobileRan_temp = mobileRan((mobile_temp+mixMN):end);
    elseif mobile_temp > MobileRange - mixMN
        mobileRan_temp = mobileRan(1:(mobile_temp-mixMN));
    else
        mobileRan_temp = [mobileRan(1:mobile_temp-mixMN) mobileRan(mobile_temp+mixMN:end)];
    end
    mobile_temp = randperm(length(mobileRan_temp),1);
      
    PHANTOM(:,i) = phantom(1:N);         %�洢����
end

figure;
imagesc(20*log(0.1 + abs(PHANTOM)));
title('PHANTON');
colormap(gray);

%% ����������
% Generate aperture for emission
% Generate aperture for emission
set_sampling(fs);
emit_aperture = xdc_linear_array (Ne, width, element_height, kerf, 1, 1,focus);
focusPlane = zeros(Ne,1);
xdc_focus_times(emit_aperture, 0, focusPlane.');
%% ����������Ӧ�Լ���������
% Set the impulse response and excitation of the emit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:8/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (emit_aperture, impulse_response);

excitation=sin(2*pi*f0*(0:1/fs:8/f0));
xdc_excitation (emit_aperture, excitation);

%% ��������
pha_Pos    = zeros(1, 3); 
fp         = zeros(N, repeatNumber);
for p = 1:N
    for q = 1:repeatNumber
        if PHANTOM(p,q) ~= 0 
            pha_Pos(1, 1) = (q - repeatNumber/2) / repeatNumber * 127 * (kerf+width) + width;
            pha_Pos(1, 2) = 0;
            pha_Pos(1, 3) = p * (z_size / N);
            [hp, start_time] = calc_hp(emit_aperture,pha_Pos);
            fp(p, q) = max(hp(:,1));
        end
    end
end
%% 
fp = abs(fp);
ffp = fp(2000:end,:);
env = ffp/max(max(ffp));
env = 100*log(env + 10);
env = env - min(min(env));
env = 64*env/max(max(env));
figure;
image(env); %axis image;
colormap(jet(64)); colorbar;

%% IQ���
[mm,nn]  = size(fp);
Data_I = zeros(nn,mm);
Data_Q = zeros(nn,mm);
t = 0:mm-1;

for i = 1:repeatNumber
    Data_I(i,:) = cos(2*pi*f0/fs*t).*(fp(:,i).');
    Data_Q(i,:) = sin(2*pi*f0/fs*t).*(fp(:,i).');
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
    slowTimeSignal_I(j,1) = sum(Data_I_fil(floor(mm/2 + MobileRange):floor(mm/2 + MobileRange*4),j));
    slowTimeSignal_Q(j,1) = sum(Data_Q_fil(floor(mm/2 + MobileRange):floor(mm/2 + MobileRange*4),j));
end

hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,1000,6e6,100e6),'butter');
slowTimeSignal_I_fil = filter(hd,slowTimeSignal_I);
slowTimeSignal_Q_fil = filter(hd,slowTimeSignal_Q);

%% ��ʱ���ʱ����Ҷ�任
clear i;
slowTimeSignal = slowTimeSignal_I + i*slowTimeSignal_Q;
hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,10,50e6,100e6),'butter');
slowTimeSignal = filter(hd,slowTimeSignal);

[S,F,T,~] = spectrogram(slowTimeSignal,256,250,256);
figure;       
S = fftshift(S,1);
P = 2000000000*log(1 + abs(S));
[x,y]=size(S); 
% subplot(Num,2,2*(cir_map-Start)+2);
image(P);
colorbar;
set(gca,'YDir','normal');
% colormap(gray);
xlabel('ʱ�� t/s');ylabel('Ƶ�� f/Hz');
axis off;
% axis([1 299 100 156]);
title(['��ʱ����ҶʱƵͼ',num2str(cir_map),'M']);
% step = 20;
% ind  = 0;
% for i = 1: step: repeatNumber - 128
%     ind = ind+1;
%     v = heartSignal(mod(i,length(heartSignal))+1) + 1;
% 
%     distancePerPulse = v/prf;  % v unit mm/s
%     n =lateralRes/distancePerPulse;
%     n = floor (n); 
%     lateralFilter = ones (1,n)/n; %����һ����Ƶ�ʼ��ٶȻ���PRF����ص���ʱ�䷽���ϵ��˲���
%  
%     slowsig = slowTimeSignal (i:i+128);
%     slowsig = conv (slowsig, lateralFilter, 'same');
%     im(:,ind)= fftshift(abs(fft(slowsig)));
% end
% im(65,:)=0;
% figure;
% imagesc(20*log(0.1+im))
% title(['PWͼ',num2str(cir_map),'M']);
end


