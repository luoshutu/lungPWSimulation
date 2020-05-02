clear all; close all; 
num_lines = 128; %每幅图像有多少条线
windowlength = 100; % 每个焦点里面有多少个散射源
phantom_pts = 128* 10000; %一个长长的一维仿体
step = 77; %每一帧，仿体移动多少个点
p_phantom_amp = randn (1,phantom_pts); %每个点的散射强度 
p_phantom_pha = rand (1,phantom_pts) * 2* pi; %每个点的相位
PSF = fspecial('gaussian', [1,windowlength], 15);  %PSF的样子。 不写SINC了， 高斯也差不多
for time = 1: 20
for i =1 : num_lines
    window = (i-1)*windowlength+1+time*step:(i-1)*windowlength+windowlength+time*step; 
    liver_phantom = p_phantom_amp(window); %肝脏的散射源多，每个细胞都是
    lung_phantom = liver_phantom; 
    lung_phantom (1:7:end) = 0;
    lung_phantom (4:7:end) = 0; %肺的散射源少。 每个肺泡才是。 
    %散射源少， 下面的求和就比较容易随机 （所以既有可追踪的滑动， 也有随机出现的亮点--蚂蚁）
    %散射源多， 求和就比较容易稳定 （所以更加的可以追踪）
    %偷懒解释： 上面代码， 其实不应该把肺仿体的一部分散射源置零。 更好的做法应该是把那些散射源的相位取他们邻居的相位
    echostrength_liver (time,i) = sum(abs(PSF.* liver_phantom.*exp(j*p_phantom_pha(window))));
    echostrength_lung (time,i) = sum(abs(PSF.* lung_phantom.*exp(j*p_phantom_pha(window))));
end
imagesc(echostrength_lung(time,:)); pause(0.1)
end 
% figure;
%   subplot(121); imagesc(echostrength_liver); title('liver sliding')
%   subplot(122); imagesc(echostrength_lung); title('lung sliding')
