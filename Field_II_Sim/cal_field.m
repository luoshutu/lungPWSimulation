function [outputArg1] = cal_field(N,z_size,lineNumber,fs,f0)
% 用于计算所需要的声场
%% 初始化
% path(path,'D:/MATLAB/Field_II_ver_3_24_windows_gcc');
% field_init;
%% 参数设置
width          = 0.3/1000;            % 阵元宽度 [m]
element_height = 5/1000;              % 阵元高度 [m]
kerf           = 0.1/1000;            % 阵元间隙宽度 [m]
Ne             = 128;                 % 阵元数量
focus          = [0 0 50]/1000;       % 固定焦点位置 [m] 
focal_depth    = focus(3);            % 焦点深度
array_size     = (kerf+width)*Ne;     % 阵元总宽度

%% 换能器设置
% Generate aperture for emission
set_sampling(fs);
emit_aperture = xdc_linear_array (Ne, width, element_height, kerf, 1, 1,focus);
% focusPlane = zeros(Ne,1);
% xdc_focus_times(emit_aperture, 0, focusPlane.');
%% 设置脉冲响应以及激励脉冲
% Set the impulse response and excitation of the emit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:8/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (emit_aperture, impulse_response);

excitation=sin(2*pi*f0*(0:1/fs:8/f0));
xdc_excitation (emit_aperture, excitation);

%% 计算声场
fieldWidth    = 15;                                   % 所计算声场的水平方向点数,奇数
fieldLength   = 15;%round((fieldWidth * (array_size / repeatNumber)) / (z_size / N));
fieldRadius   = floor(fieldWidth / 2);                % 所计算声场的水平方向半径
pha_Pos       = zeros(1, 3);                          % 计算声场的点位置
depthCenter   = round(focal_depth / (z_size / N));    % 整个发射声场的深度方向中心点
horCenter     = round(lineNumber / 2);                % 整个发射声场的水平方向中心点
fp            = zeros(fieldLength, fieldWidth);       % 声场缓存 

if mod(fieldLength, 2) == 0
    fieldLength = fieldLength + 1;
end

ind = 0;
for p = depthCenter - floor(fieldLength / 2):depthCenter + floor(fieldLength / 2)
    ind = ind + 1;
    for q = horCenter - fieldRadius:horCenter + fieldRadius
        pha_Pos(1, 1) = (q - lineNumber/2) * (array_size / lineNumber);
        pha_Pos(1, 2) = 0;
        pha_Pos(1, 3) = p * (z_size / N);
        [hp, ~] = calc_hp(emit_aperture,pha_Pos);
        fp(ind, q - (horCenter - fieldRadius) + 1) = max(hp(:,1));
    end
end

outputArg1 = fp;
end

