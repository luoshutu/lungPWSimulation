clear all; close all; 
num_lines = 128; %ÿ��ͼ���ж�������
windowlength = 100; % ÿ�����������ж��ٸ�ɢ��Դ
phantom_pts = 128* 10000; %һ��������һά����
step = 77; %ÿһ֡�������ƶ����ٸ���
p_phantom_amp = randn (1,phantom_pts); %ÿ�����ɢ��ǿ�� 
p_phantom_pha = rand (1,phantom_pts) * 2* pi; %ÿ�������λ
PSF = fspecial('gaussian', [1,windowlength], 15);  %PSF�����ӡ� ��дSINC�ˣ� ��˹Ҳ���
for time = 1: 20
for i =1 : num_lines
    window = (i-1)*windowlength+1+time*step:(i-1)*windowlength+windowlength+time*step; 
    liver_phantom = p_phantom_amp(window); %�����ɢ��Դ�࣬ÿ��ϸ������
    lung_phantom = liver_phantom; 
    lung_phantom (1:7:end) = 0;
    lung_phantom (4:7:end) = 0; %�ε�ɢ��Դ�١� ÿ�����ݲ��ǡ� 
    %ɢ��Դ�٣� �������;ͱȽ�������� �����Լ��п�׷�ٵĻ����� Ҳ��������ֵ�����--���ϣ�
    %ɢ��Դ�࣬ ��;ͱȽ������ȶ� �����Ը��ӵĿ���׷�٣�
    %͵�����ͣ� ������룬 ��ʵ��Ӧ�ðѷη����һ����ɢ��Դ���㡣 ���õ�����Ӧ���ǰ���Щɢ��Դ����λȡ�����ھӵ���λ
    echostrength_liver (time,i) = sum(abs(PSF.* liver_phantom.*exp(j*p_phantom_pha(window))));
    echostrength_lung (time,i) = sum(abs(PSF.* lung_phantom.*exp(j*p_phantom_pha(window))));
end
imagesc(echostrength_lung(time,:)); pause(0.1)
end 
% figure;
%   subplot(121); imagesc(echostrength_liver); title('liver sliding')
%   subplot(122); imagesc(echostrength_lung); title('lung sliding')
