function [outputArg1] = phantom_line(N, ster_line_mobile, seed)
% PHANTOM_LINE ���ڲ������壬�·β�����һ���ߵķ��塣
% NΪ���������ster_line_mobileΪ��ģ�������ƶ�λ�ã�
% seedΪ������ӣ���֤��ģ���ϵķ���ÿһ�����в�𣬵�ÿһ֮֡����������
%   �˴���ʾ��ϸ˵��
phantom        = zeros(N,1);
rng(seed);
amp            = rand(N, 1);                      % ��ģ��֮�Ϸ�����������

rng('shuffle');                                     % ��֤ÿ�������λ������ͬ
p_phantom_pha  = rand(N/2, 1) * 2 * pi;             % ����β�ÿ�������λ

% skin
rng(seed);
% z_skin = randperm(round(N*(30/1000)),round(N*(30/1000)));
for idx_skin = 1:round(N*(30/1000))
%     phantom(idx_skin,1) = amp(idx_skin,1) + 0.5;
    phantom(idx_skin,1) = 2;
end

% fat
z_fat = randperm(round(N*(170/1000)),round(N*(70/1000))) + round(N*(30/1000));
for idx_fat = z_fat
	phantom(idx_fat,1) = amp(idx_fat,1)*2; 
end

%muscle
z_muscle = randperm(round(N*(400/1000)),round(N*(160/1000))) + round(N*(200/1000));
for idx_muscle = z_muscle
	phantom(idx_muscle,1) = amp(idx_muscle,1); 
end

POINT_position = round(N/2 - ster_line_mobile); %THE POINT��λ��

% ��ģ���·�����ɢ��ԴΪ��֮�Ϸ���ľ����һ�������ģ�ߵ��˶����˶�
for p = 1:(POINT_position)
    phantom(POINT_position + p) = phantom(POINT_position - p + 1) ;
end

% ��ģ�߷�ֵ�ı仯��Ϊ1��5��֮ǰ�������ֵ֮������ֵ
% POINT_Am  = (randperm(3,1))*max(phantom);
POINT_Am = 10;
phantom(POINT_position,1) = POINT_Am;
phantom = phantom(1:N);

if POINT_position+(N/2 - 1) > N
    phantom(POINT_position:end) = phantom(POINT_position:end).*exp(1i*p_phantom_pha(1:(N - POINT_position + 1)));
else
    phantom(POINT_position:POINT_position+N/2-1) = phantom(POINT_position:POINT_position+N/2-1).*exp(1i*p_phantom_pha);
end

outputArg1 = phantom;
end

