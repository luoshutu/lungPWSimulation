%% 功能：实现BEE8的单通道M模式成像
% By Luoshutu
%% 清理空间
clc;
close all;
clear;

%% Read File
framenum    = 0:9;
linenum     = 2560;
channelnum  = 8;
packetnum   = 320;
pointnum    = 2048;

path     = 'D:\_MyProject\MATLAB\history_blood';
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

%% 解调
data_m = zeros(linenum*length(framenum)/8,pointnum);
for i = 1:linenum*length(framenum)
    if(rem(i,8) == 4)
        data_m((i + 4)/8,:) = Data(i,:);
    end
end
    
Data_I = zeros(linenum*length(framenum)/8,pointnum);
Data_Q = zeros(linenum*length(framenum)/8,pointnum);

n1 = 0:pointnum-1;

for i = 1:linenum*length(framenum)/8
    
    Data_I(i,:) = cos(2*pi*4e6/40e6*n1).*(data_m(i,:));
    Data_Q(i,:) = sin(2*pi*4e6/40e6*n1).*(data_m(i,:));
    
end

%% 求幅值
Data_Amplitude = sqrt((Data_I).*(Data_I) + (Data_Q).*(Data_Q));

%% 对数变换
Image_Initial = 40*log(50+abs(Data_Amplitude));
% figure(5);
% imagesc(Image_Initial);                    
% colormap('gray');
%% 归一化
Am_max = max(max(Image_Initial));
Am_min = min(min(Image_Initial));
Image_Final1 = uint16(round((255-0).*(Image_Initial-Am_min)./(Am_max-Am_min)+0));
figure(1);
image(Image_Final1(:,:).');
colormap(gray(256));


