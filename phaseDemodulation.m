%% �����źŵĽ��
%% ����
clear
clc
close all

%% ����
f0 = 2;     % ԭʼ�ź�Ƶ�ʣ���λHz
fc = 1e3;   % �ز��ź�Ƶ�ʣ���λHz
fs = 1e5;   % �źŲ���Ƶ�ʣ���λHz

%% �����ź�
t           = 0:(1/fs):8;                   % ʱ��
sigSrc      = 0.5*sin(2*pi*f0.*t) + 0.5*sin(2*pi*2*f0.*t) + 0.3*cos(2*pi*0.5*f0.*t+pi/4);       % ԭʼ�ź�
sigCar      = sin(2*pi*fc.*t);              % �ز��ź�
sigPhaMod	= sin(2*pi*fc.*t + 100*sigSrc);     % �����ź�

figure;
plot(sigSrc,'r');

%% IQ���
% sigDemod = asin(sigPhaMod) - asin(sigCar);
sigDemod = sigPhaMod;

Data_I = cos(2*pi*fc/fs.*t).*(sigDemod);
Data_Q = sin(2*pi*fc/fs.*t).*(sigDemod);

Data_e = Data_I + i.*Data_Q;
% %% �˲�
% hd   = design(fdesign.lowpass('N,F3dB',16,100,fs),'butter');
% Data_I_fil = filter(hd, Data_I);
% Data_Q_fil = filter(hd, Data_Q);

% Data_a = sqrt((Data_I_fil.*Data_I_fil)+(Data_Q_fil.*Data_Q_fil));

%% ���
z  = hilbert(sigPhaMod);                 %ϣ�����ر任�Ե�����---ͨ��ʵ�������鲿
x1 = z.*exp(-j*2*pi*fc*t);
% x1    = Data_e.*exp(-j*2*pi*fc*t);
pha   = angle(x1);    % ��ȡ�����źŵ���λ
demod = (1/100)*unwrap(pha);      % ��λУ��  
figure;
plot(demod);
hold on;
plot(sigSrc);

% hd   = design(fdesign.lowpass('N,F3dB',16,100,fs),'butter');
% demod_fil = filter(hd, demod);
% figure;
% plot(demod_fil);