%% Ŀ�ģ�ʹ��Field IIʵ�ַβ���PW�źŷ��� 
% ��ά����
% By Luoshutu
clear; clc;

%% ��ʼ��
path(path,'D:/MATLAB/Field_II_ver_3_24_windows_gcc');
field_init;

%% ���÷�����ʼ
Start = 3;
Num   = 1;
% figure;
% hold on;
for cir_map = Start:Start + Num -1
%% ��������
f0             = cir_map*(1e6);
% f0          = 7e6;      % ����������Ƶ�� [Hz]
fs             = 100e6;	% ����ʹ�õĲ���Ƶ��
c              = 1540;     % ���� [m/s]
lambda         = c/f0;     % ���� [m]
prf            = 500;      % �����ظ�Ƶ��Hz/s
z_size         = 60/1000;  % �������[m]

width          = 0.3/1000;            % ��Ԫ��� [m]
element_height = 5/1000;              % ��Ԫ�߶� [m]
kerf           = 0.1/1000;            % ��Ԫ��϶��� [m]
Ne             = 128;                 % ��Ԫ����
focus          = [0 0 50]/1000;       % �̶�����λ�� [m] 
focal_depth    = focus(3);            % �������
array_size     = (kerf+width)*Ne;     % ��Ԫ�ܿ��

%% ���ط���
load phantom;
phantomSrc = phantom;
phantom = phantomSrc(1:3:8000,1:50:end);

[m_p, n_p]   = size(phantom);

%% �����ź�����
heartRate = 80/60;  % ����Ƶ�ʣ�ÿ����80��
heartNum  = ceil(prf / heartRate); % һ��������ռ��������
heartSignalAmp = round([0.00 1.00 0.95 0.73 0.20 -0.20 0.06 0.08 -0.03 0.00].*(1e-4 / (z_size / m_p)));
heartSignal = interp(heartSignalAmp,round(heartNum/length(heartSignalAmp)));
% figure;
% plot(heartSignal);
% grid on;

%% humming�ź�(�ߺ�)
f_hum   = 80;                               % �ߺ�Ƶ�ʣ���λHz
v_hum   = 3e-5;                             % �ߺ�����ֵ����λm
% humming = (v_hum/(z_size/m_p)) * sin(2*pi*f_hum*(0:1/prf:3/f_hum));
% humming = (v_hum/(z_size/m_p)) * sawtooth(2*pi*f_hum*(0:1/prf:3/f_hum),0.4);
humming = (v_hum/(z_size/m_p)) * square(2*pi*f_hum*(0:1/prf:3/f_hum),50);        % ռ�ձ�Ϊ50%�ķ���
% humming = (v_hum/(z_size/m_p)) * (2 * rand(1,round(prf/f_hum)) - 1);
% figure;
% plot(humming);
% title('Humming');
% grid on;
%% ��������
excitation_pulse = sin(2*pi*f0*(0:1/fs:8/f0));
excitation_pulse = excitation_pulse.*hanning(max(size(excitation_pulse)))';

%% ���
echo         = zeros(m_p + length(excitation_pulse) - 1,n_p);  % ����������

for i = 1:n_p                              % �������
    % ���
    echo(:,i) = conv(excitation_pulse,phantom(:,i).').';
end

% figure;
% imagesc(20*log(0.1 + abs(phantom)));
% title('PHANTON');
% colormap(gray);
% figure;
% imagesc(20*log(0.1 + abs(echo)));
% title('Echo');
% colormap(gray);

%% ����������
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

%% ����תɢ��Դ����
phantom_positions  = zeros(m_p * n_p,3);

for i = 1:m_p
    for j = 1:n_p
        phantom_positions((i-1) * n_p + j, 1) = (j - n_p/2) * (array_size / n_p);
        phantom_positions((i-1) * n_p + j, 3) = i * (kerf+width);
    end
end

phantom_amplitudes = phantom';
phantom_amplitudes = phantom_amplitudes(:);
% phantom_amplitudes = double(255 - phantom_amplitudes);
imgFFF = reshape(phantom_amplitudes, n_p, m_p)';
figure;
imagesc(abs(imgFFF));
% axis image;
axis off;
colormap(gray);
%% ��������

[v,t]= calc_scat_multi(emit_aperture,emit_aperture,phantom_positions,phantom_amplitudes);%[0 0 0; 0 0 15]/1000,[0 ;3]);
[N,M]= size(v);
v=v/max(v(:));

%% ��ͼ
% fp = abs(v);
% env = fp/max(max(fp));
% env = log(env + 0.1);
% env = env - min(min(env));
% env = 64*env/max(max(env));
% figure;
% image(env); %axis image;
% colormap(jet(64)); colorbar;

%% ��������
speed               = c;                            % ��������1540
txFrequency         = f0;                           % ��������Ƶ��
channel_number      = 32;                           % ͨ����Ŀ32
channel_number_half = channel_number/2;             % ͨ����Ŀһ��
beam_number         = 4;                            % ��������4
point_number        = n_p;                            % �����������2944
point_space         = c/(2*fs);
% point_space         = (m * (kerf+width))/N;         % �����������2.4640e-05
ele_number          = Ne;                           % ��Ԫ��Ŀ128
ele_space           = kerf + width;                 % ��Ԫ���3.0*e-04
rx_fn               = 2;                          % ����F��2.8
%% ��ʱ
Delay = zeros(channel_number, point_number, beam_number);

for i = 1:beam_number
    for j = 1:point_number
        for k = 1:channel_number
            focus_depth  = j * point_space;
            ele_distance = abs(k - (channel_number_half + (i - 1) / beam_number)) * ele_space;
            Delay(k, j, i) = int16((sqrt(focus_depth * focus_depth + ele_distance * ele_distance - 2 * focus_depth * ele_distance *cos(pi/2 - 0)) - focus_depth) / speed * fs);
        end
    end
end

%% �׾��ͼ�Ȩ
Apodization_base    = hamming(channel_number);
Apodization         = zeros(channel_number,point_number);   % ��Ȩ - ��С�԰�
Aperture            = zeros(channel_number,point_number);   % �׾� - ��ͬ��ȶ�Ӧ��ͬ�Ŀ׾���С

ape_ele_number = zeros(1,point_number);
for x = 1:point_number
    ape_ele_number(1,x) = int16((x * point_space / rx_fn) / ele_space) + 1;
    if ape_ele_number(1,x) <= 32
        for y = 1:ape_ele_number(1,x)
            if mod(y, 2) == 0
                Aperture(channel_number_half + (y / 2), x) = 1;
            else
                Aperture(channel_number_half - floor(y / 2), x) = 1;
            end
        end
    else
        Aperture(1:channel_number, x) = 1;
    end
end

for z = 1:point_number
    Apodization(:, z) = Apodization_base .* Aperture(:, z);
end
%% Beamforming
Data                = v;                            % ��������

RFData              = zeros(point_number,ele_number * beam_number);  % ��ʼ��RF����
tempZero            = zeros(point_number,channel_number_half);       % ��������
data                = [tempZero Data tempZero];     % ����������

for m = 1:ele_number
    for n = 1:beam_number
        for p = 1:point_number
            preSynthesisData = zeros(1,channel_number);
            for q = 1:channel_number
                delay_number = Delay(q, p, n) + p;
                if(delay_number > point_number)
                    delay_number = point_number;
                end
                preSynthesisData(1, q) = data(delay_number, m + q);
            end
            RFData(p, (m - 1) * beam_number + n) = preSynthesisData * Apodization(:, p);
        end
    end
end

figure;
img_rf = abs(hilbert(RFData));
imagesc(img_rf);
% set(gcf,'Units','centimeter','Position',[5 2 11 11*((N*point_space)/(M*(kerf + width)))]);%0.076 0.0512
set(gcf,'Units','centimeter','Position',[5 2 11 11*(200/190)]);
% colormap(jet(64));
axis off;
colormap(gray);
%% �ͷ�Ϊ�������������Դ
xdc_free (emit_aperture);

%% IQ���
[mm,nn]  = size(echo);
Data_I = zeros(nn,mm);
Data_Q = zeros(nn,mm);
t = 0:mm-1;

for i = 1:n_p
    Data_I(i,:) = cos(2*pi*f0/fs*t).*(echo(:,i).');
    Data_Q(i,:) = sin(2*pi*f0/fs*t).*(echo(:,i).');
end

% �˲�
hd   = design(fdesign.lowpass('N,F3dB',16,5e6,100e6),'butter');
Data_I_fil = filter(hd,Data_I.');
Data_Q_fil = filter(hd,Data_Q.');

%% ��ʾMģʽ
Data_Amp = sqrt((Data_I_fil).*(Data_I_fil) + (Data_Q_fil).*(Data_Q_fil));
D=20;   % �²�����ȡ��
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
P = 20*log(10 + abs(S));
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

lateralRes  = 1 / (f0/0.5e6); % unit : mm ����3��ʱ������ֱ���Ϊ1����
step = 5;
ind  = 0;
for i = 1: step: repeatNumber - 256
    ind = ind+1;
    v = heartSignal(mod(i,length(heartSignal))+1);

    distancePerPulse = v/prf;  % v unit mm/s
    n = lateralRes/distancePerPulse;
    n = floor (n); 
    nn(i) = n;
    lateralFilter = sinc(-abs(n)/2 :0.1: abs(n)/2)/n; % ����һ����Ƶ�ʼ��ٶȻ���PRF����ص���ʱ�䷽���ϵ��˲���
 
    slowsig = slowTimeSignal (i:i+256);
%     slowsig = echo(5100,i:i+128);
%     slowsig = conv (slowsig, lateralFilter, 'same');
    im(:,ind)= fftshift(abs(fft(slowsig)));
end
% im(65,:)=0;
figure;
imagesc(20*log(0.1+im))
title(['PWͼ',num2str(cir_map),'M']);
end

