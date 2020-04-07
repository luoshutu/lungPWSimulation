%% 目的：使用Field II实现肺部的PW信号仿真 
% 二维仿真
% By Luoshutu
clear; clc;

%% 初始化
path(path,'D:/MATLAB/Field_II_ver_3_24_windows_gcc');
field_init;

%% 设置仿真起始
Start = 3;
Num   = 1;
% figure;
% hold on;
for cir_map = Start:Start + Num -1
%% 参数设置
f0             = cir_map*(1e6);
% f0          = 7e6;      % 换能器中心频率 [Hz]
fs             = 100e6;	% 仿真使用的采样频率
c              = 1540;     % 声速 [m/s]
lambda         = c/f0;     % 波长 [m]
prf            = 500;      % 脉冲重复频率Hz/s
z_size         = 60/1000;  % 仿体深度[m]

width          = 0.3/1000;            % 阵元宽度 [m]
element_height = 5/1000;              % 阵元高度 [m]
kerf           = 0.1/1000;            % 阵元间隙宽度 [m]
Ne             = 128;                 % 阵元数量
focus          = [0 0 50]/1000;       % 固定焦点位置 [m] 
focal_depth    = focus(3);            % 焦点深度
array_size     = (kerf+width)*Ne;     % 阵元总宽度

%% 加载仿体
load phantom;
phantomSrc = phantom;
phantom = phantomSrc(1:3:8000,1:50:end);

[m_p, n_p]   = size(phantom);

%% 心跳信号设置
heartRate = 80/60;  % 心跳频率，每分钟80次
heartNum  = ceil(prf / heartRate); % 一次心跳所占采样线数
heartSignalAmp = round([0.00 1.00 0.95 0.73 0.20 -0.20 0.06 0.08 -0.03 0.00].*(1e-4 / (z_size / m_p)));
heartSignal = interp(heartSignalAmp,round(heartNum/length(heartSignalAmp)));
% figure;
% plot(heartSignal);
% grid on;

%% humming信号(哼哼)
f_hum   = 80;                               % 哼哼频率，单位Hz
v_hum   = 3e-5;                             % 哼哼最大幅值，单位m
% humming = (v_hum/(z_size/m_p)) * sin(2*pi*f_hum*(0:1/prf:3/f_hum));
% humming = (v_hum/(z_size/m_p)) * sawtooth(2*pi*f_hum*(0:1/prf:3/f_hum),0.4);
humming = (v_hum/(z_size/m_p)) * square(2*pi*f_hum*(0:1/prf:3/f_hum),50);        % 占空比为50%的方波
% humming = (v_hum/(z_size/m_p)) * (2 * rand(1,round(prf/f_hum)) - 1);
% figure;
% plot(humming);
% title('Humming');
% grid on;
%% 激励函数
excitation_pulse = sin(2*pi*f0*(0:1/fs:8/f0));
excitation_pulse = excitation_pulse.*hanning(max(size(excitation_pulse)))';

%% 卷积
echo         = zeros(m_p + length(excitation_pulse) - 1,n_p);  % 卷积结果缓存

for i = 1:n_p                              % 卷积次数
    % 卷积
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

%% 换能器设置
% Generate aperture for emission
set_sampling(fs);
emit_aperture = xdc_linear_array (Ne, width, element_height, kerf, 1, 1,focus);
focusPlane = zeros(Ne,1);
xdc_focus_times(emit_aperture, 0, focusPlane.');

%% 设置脉冲响应以及激励脉冲
% Set the impulse response and excitation of the emit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:8/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (emit_aperture, impulse_response);

excitation=sin(2*pi*f0*(0:1/fs:8/f0));
xdc_excitation (emit_aperture, excitation);

%% 仿体转散射源坐标
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
%% 计算声场

[v,t]= calc_scat_multi(emit_aperture,emit_aperture,phantom_positions,phantom_amplitudes);%[0 0 0; 0 0 15]/1000,[0 ;3]);
[N,M]= size(v);
v=v/max(v(:));

%% 画图
% fp = abs(v);
% env = fp/max(max(fp));
% env = log(env + 0.1);
% env = env - min(min(env));
% env = 64*env/max(max(env));
% figure;
% image(env); %axis image;
% colormap(jet(64)); colorbar;

%% 参数定义
speed               = c;                            % 传播声速1540
txFrequency         = f0;                           % 超声发射频率
channel_number      = 32;                           % 通道数目32
channel_number_half = channel_number/2;             % 通道数目一半
beam_number         = 4;                            % 波束个数4
point_number        = n_p;                            % 纵向采样点数2944
point_space         = c/(2*fs);
% point_space         = (m * (kerf+width))/N;         % 纵向采样点间距2.4640e-05
ele_number          = Ne;                           % 阵元数目128
ele_space           = kerf + width;                 % 阵元间距3.0*e-04
rx_fn               = 2;                          % 接收F数2.8
%% 延时
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

%% 孔径和加权
Apodization_base    = hamming(channel_number);
Apodization         = zeros(channel_number,point_number);   % 加权 - 减小旁瓣
Aperture            = zeros(channel_number,point_number);   % 孔径 - 不同深度对应不同的孔径大小

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
Data                = v;                            % 接收数据

RFData              = zeros(point_number,ele_number * beam_number);  % 初始化RF数组
tempZero            = zeros(point_number,channel_number_half);       % 补零数据
data                = [tempZero Data tempZero];     % 补零后的数据

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
%% 释放为换能器分配的资源
xdc_free (emit_aperture);

%% IQ解调
[mm,nn]  = size(echo);
Data_I = zeros(nn,mm);
Data_Q = zeros(nn,mm);
t = 0:mm-1;

for i = 1:n_p
    Data_I(i,:) = cos(2*pi*f0/fs*t).*(echo(:,i).');
    Data_Q(i,:) = sin(2*pi*f0/fs*t).*(echo(:,i).');
end

% 滤波
hd   = design(fdesign.lowpass('N,F3dB',16,5e6,100e6),'butter');
Data_I_fil = filter(hd,Data_I.');
Data_Q_fil = filter(hd,Data_Q.');

%% 显示M模式
Data_Amp = sqrt((Data_I_fil).*(Data_I_fil) + (Data_Q_fil).*(Data_Q_fil));
D=20;   % 下采样抽取率
figure;
% subplot(Num,2,2*(cir_map-Start)+1);
imagesc(20*log(1 + abs(Data_Amp(1:D:end,:))));
title(['M模式',num2str(cir_map),'M']);
colormap(gray);

%% 累加
slowTimeSignal_I = zeros(repeatNumber,1);
slowTimeSignal_Q = zeros(repeatNumber,1);

for j = 1:repeatNumber
    slowTimeSignal_I(j,1) = sum(Data_I_fil(floor(mm/2 + MobileRange):floor(mm/2 + MobileRange*2),j));
    slowTimeSignal_Q(j,1) = sum(Data_Q_fil(floor(mm/2 + MobileRange):floor(mm/2 + MobileRange*2),j));
end

% hd   = design(fdesign.bandpass('N,F3dB1,F3dB2',16,1000,6e6,100e6),'butter');
% slowTimeSignal_I = filter(hd,slowTimeSignal_I);
% slowTimeSignal_Q = filter(hd,slowTimeSignal_Q);

%% 慢时间短时傅里叶变换
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
xlabel('时间 t/s');ylabel('频率 f/Hz');
axis off;
% axis([1 y x/2 - 20 x/2 + 20]);
title(['短时傅里叶时频图',num2str(cir_map),'M']);

lateralRes  = 1 / (f0/0.5e6); % unit : mm 假设3兆时，横向分辨率为1毫米
step = 5;
ind  = 0;
for i = 1: step: repeatNumber - 256
    ind = ind+1;
    v = heartSignal(mod(i,length(heartSignal))+1);

    distancePerPulse = v/prf;  % v unit mm/s
    n = lateralRes/distancePerPulse;
    n = floor (n); 
    nn(i) = n;
    lateralFilter = sinc(-abs(n)/2 :0.1: abs(n)/2)/n; % 构建一个和频率及速度还有PRF都相关的慢时间方向上的滤波器
 
    slowsig = slowTimeSignal (i:i+256);
%     slowsig = echo(5100,i:i+128);
%     slowsig = conv (slowsig, lateralFilter, 'same');
    im(:,ind)= fftshift(abs(fft(slowsig)));
end
% im(65,:)=0;
figure;
imagesc(20*log(0.1+im))
title(['PW图',num2str(cir_map),'M']);
end

