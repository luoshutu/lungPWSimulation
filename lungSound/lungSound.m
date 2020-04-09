%% Ŀ�ģ�ʵ�ַβ��źŵĽ��
% By Luoshutu
clear; clc;
%% ��������
f0          = 3.5e6;    % ����������Ƶ�� [Hz]
fs          = 40e6;     % ����Ƶ��
prf         = 1000;     % �����ظ�Ƶ��Hz/s

%% ��������
load data_m

%% IQ���
[mm,nn]  = size(data_m);
Data_I = zeros(mm,nn);
Data_Q = zeros(mm,nn);
t = 0:nn-1;
for i = 1:mm
    Data_I(i,:) = cos(2*pi*f0/fs*t).*(data_m(i,:));
    Data_Q(i,:) = sin(2*pi*f0/fs*t).*(data_m(i,:));
end

%% ��ʾMģʽ
hd   = design(fdesign.lowpass('N,F3dB',16,5e6,100e6),'butter');
Data_I_fil = filter(hd,Data_I.');
Data_Q_fil = filter(hd,Data_Q.');
Data_Amp = sqrt((Data_I_fil).*(Data_I_fil) + (Data_Q_fil).*(Data_Q_fil));

figure;
imagesc(20*log(1 + abs(Data_Amp)));
title('Mģʽ');
colormap(gray);

%% ��λ���
xx = 0:(1/prf):(1/prf)*(mm-1);
z  = hilbert(sum(real(data_m(:,400:420)')));                 %ϣ�����ر任�Ե�����---ͨ��ʵ�������鲿
% z = sum(data_m(:,400:420)');
x1 = z.*exp(-j*2*pi*prf*xx);
% x1    = Data_e.*exp(-j*2*pi*fc*t);
pha   = angle(x1);          % ��ȡ�����źŵ���λ
demod = unwrap(pha);        % ��λУ��  
hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,1e4,20e6,100e6),'butter');
demod_fil = filter(hd,demod);
figure;
plot((demod_fil));
%% ��arctan2
for xy = 1:mm
    De(xy) = atan2(sum(real(Data_I(xy,400:420))),sum(real(Data_Q(xy,400:420))));
end
Dt = unwrap(De);
hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,1e4,10e6,100e6),'butter');
Dt_fil = filter(hd,Dt);
figure
plot(Dt,'b');
hold on;
plot(Dt_fil,'r');
legend('�˲�ǰ','�˲���');

