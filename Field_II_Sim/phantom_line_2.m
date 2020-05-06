function [outputArg1] = phantom_line_2(N, ster_line_mobile, seed, lungSlip)
% PHANTOM_LINE 用于产生仿体，仿肺部纵向一条线的仿体。
% N为仿体点数，ster_line_mobile为胸模线上下移动位置，
% seed为随机种子，保证胸模线上的仿体每一条线有差别，但每一帧之间无随机差别
% lungSlip为肺滑距离
%   此处显示详细说明
phantom        = zeros(N,1);
rng(seed);
amp            = rand(N, 1);                      % 胸模线之上仿体的随机幅度

rng(seed+lungSlip);                               % 保证每次随机相位都不相同
p_phantom_pha  = rand(N, 1) * 2 * pi;             % 仿体肺部每个点的相位

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

POINT_position = round(N/2 - ster_line_mobile - 5); %THE POINT的位置

% 随胸模线运动，脂肪以及肌肉的一些轻微运动
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

% 胸模线下方仿体散射源为其之上仿体的镜像，且会随着胸模线的运动而运动
for p = 1:(POINT_position)
    phantom(POINT_position + p) = phantom(POINT_position - p + 1) ;
end

% 胸模线幅值的变化，为1到5倍之前仿体最大值之间的随机值
% POINT_Am  = (randperm(3,1))*max(phantom);
POINT_Am = 4;
phantom(POINT_position - 3:POINT_position,1) = POINT_Am;
phantom = phantom(1:N);

phantom(POINT_position-2:N) = phantom(POINT_position-2:N).*exp(1i*p_phantom_pha(POINT_position-2:N));

outputArg1 = phantom;
end

