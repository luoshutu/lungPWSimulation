function [outputArg1] = phantom_line_2(N, ster_line_mobile, seed, lungSlip)
% PHANTOM_LINE ���ڲ������壬�·β�����һ���ߵķ��塣
% NΪ���������ster_line_mobileΪ��ģ�������ƶ�λ�ã�
% seedΪ������ӣ���֤��ģ���ϵķ���ÿһ�����в�𣬵�ÿһ֮֡����������
% lungSlipΪ�λ�����
%   �˴���ʾ��ϸ˵��
phantom        = zeros(N,1);
rng(seed);
amp            = rand(N, 1);                      % ��ģ��֮�Ϸ�����������

rng(seed+lungSlip);                               % ��֤ÿ�������λ������ͬ
p_phantom_pha  = rand(N, 1) * 2 * pi;             % ����β�ÿ�������λ

% skin
rng(seed);
z_skin = randperm(round(N*(30/1000)),round(N*(26/1000)));
for idx_skin = z_skin
    phantom(idx_skin,1) = amp(idx_skin,1) + 1;
%     phantom(idx_skin,1) = 2;
end

% fat
z_fat = randperm(round(N*(170/1000)),round(N*(70/1000))) + round(N*(30/1000));
for idx_fat = z_fat
    phantom(idx_fat,1) = amp(idx_fat,1)*2; 
end

%muscle
z_muscle = randperm(round(N*(600/1000)),round(N*(400/1000))) + round(N*(200/1000));
for idx_muscle = z_muscle
    phantom(idx_muscle,1) = amp(idx_muscle,1); 
end

POINT_position = round(N/2 - ster_line_mobile - 5); %THE POINT��λ��

% ����ģ���˶���֬���Լ������һЩ��΢�˶�
sub_line = N/2 - POINT_position;
cnt = floor(POINT_position/2);
for idx_m = floor(POINT_position/2):POINT_position-1
    sub_interval = round(((POINT_position-1)-floor(POINT_position/2))/sub_line);
    if mod(idx_m-floor(POINT_position/2)+1,sub_interval) == 0
        cnt = cnt + 1;
    end  
    phantom(idx_m) = phantom(cnt);
    cnt = cnt + 1;
end

% ��ģ���·�����ɢ��ԴΪ��֮�Ϸ���ľ����һ�������ģ�ߵ��˶����˶�
for p = 1:(POINT_position)
    phantom(POINT_position + p) = phantom(POINT_position - p + 1) ;
end

% ��ģ�߷�ֵ�ı仯��Ϊ1��5��֮ǰ�������ֵ֮������ֵ
% POINT_Am  = (randperm(3,1))*max(phantom);
POINT_Am = 4;
phantom(POINT_position - 3:POINT_position,1) = POINT_Am;
phantom = phantom(1:N);

phantom(POINT_position-2:N) = phantom(POINT_position-2:N).*exp(1i*p_phantom_pha(POINT_position-2:N));

outputArg1 = phantom;
end

