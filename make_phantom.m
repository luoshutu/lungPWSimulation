%% 清理
clear;
clc;

%% 读入框架图
phantom_im  = 255 - rgb2gray(imread('phant.jpg'));
[m_p, n_p] = size(phantom_im);
for i = 1:m_p
    for j = 1:n_p
        if phantom_im(i, j) < 200
            phantom_im(i, j) = 0;
        end
        % 将所有像素值限定在235以下，方便之后添加分层线
        if phantom_im(i, j) > 0
            phantom_im(i, j) = phantom_im(i, j) - 20;
        end
    end
end
% figure;
% image(phantom_im);
% axis('image');

phantom = double(phantom_im);
%% 添加分层线
seed = [1 310 1190 3870];
phantom(seed(1),:) = 255;    % 探头与皮肤

% 皮肤与脂肪
phantom(seed(2),1) = 255;
for i = 2:n_p
    [~, y] = max(phantom(seed(2) - 2:seed(2) + 2, i));
    seed(2) = y + seed(2) - 3;
    phantom(seed(2), i) = 255;
end

% 脂肪与肌肉
phantom(seed(3),1) = 255;
for i = 2:n_p
    [~, y] = max(phantom(seed(3) - 2:seed(3) + 2, i));
    seed(3) = y + seed(3) - 3;
    phantom(seed(3), i) = 255;
end

% 胸膜线
phantom(seed(4),1) = 255;
for i = 2:n_p
    [~, y] = max(phantom(seed(4) - 2:seed(4) + 2, i));
    seed(4) = y + seed(4) - 3;
    phantom(seed(4), i) = 255;
end

%% 添加散射源
for i = 1:n_p
    y = find(phantom(:, i) == 255);  % 找到每一条深度线上对应分层线的点的坐标
    
    % skin
    len_skin = y(2) - y(1) +1;
    z_skin = randperm(len_skin, round(len_skin * (95/100)));
    rand('seed',i);
    phantSkin(z_skin) = uint8(255 * rand(length(z_skin),1));
    for j = z_skin
        if phantom(j, i) == 0
            phantom(j, i) = phantSkin(j);
        end
    end

    % fat
    len_fat = y(3) - y(2) + 1;
    z_fat = randperm(len_fat, round(len_fat * (90/100))) + y(2);
    rand('seed',i+0.1);
    phantFat(z_fat) = uint8(100 * rand(length(z_fat),1));
    for p = z_fat
        if phantom(p, i) == 0
            phantom(p, i) = phantFat(p);
        end
    end

    % muscle
    len_muscle = y(4) - y(3) + 1;
    z_muscle = randperm(len_muscle, round(len_muscle * (95/100))) + y(3);
    rand('seed',i+0.2);
    phantMuscle(z_muscle) = uint8(150 * rand(length(z_muscle),1));
    for q = z_muscle
        if phantom(q, i) == 0
            phantom(q, i) = phantMuscle(q);
        end
    end

    % lung
    rand('seed',i+0.3);
    p_phantom_pha  = rand(ceil(m_p / 2), 1) * 2 * pi;                                % 每个点的相位
    clear j;
    for r = 1:y(4) - 1
        phantom(y(4) + r, i) = phantom(y(4) - r, i);
    end
    phantom(y(4):2*y(4), i) = phantom(y(4):2*y(4), i).*exp(j*p_phantom_pha(1:y(4)+1));
end

figure;
image(abs(phantom));
axis('image');

save phantom;
