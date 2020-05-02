function [outputArg1] = cal_field(N,z_size,lineNumber,fs,f0)
% ���ڼ�������Ҫ������
%% ��ʼ��
% path(path,'D:/MATLAB/Field_II_ver_3_24_windows_gcc');
% field_init;
%% ��������
width          = 0.3/1000;            % ��Ԫ��� [m]
element_height = 5/1000;              % ��Ԫ�߶� [m]
kerf           = 0.1/1000;            % ��Ԫ��϶��� [m]
Ne             = 128;                 % ��Ԫ����
focus          = [0 0 50]/1000;       % �̶�����λ�� [m] 
focal_depth    = focus(3);            % �������
array_size     = (kerf+width)*Ne;     % ��Ԫ�ܿ��

%% ����������
% Generate aperture for emission
set_sampling(fs);
emit_aperture = xdc_linear_array (Ne, width, element_height, kerf, 1, 1,focus);
% focusPlane = zeros(Ne,1);
% xdc_focus_times(emit_aperture, 0, focusPlane.');
%% ����������Ӧ�Լ���������
% Set the impulse response and excitation of the emit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:8/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (emit_aperture, impulse_response);

excitation=sin(2*pi*f0*(0:1/fs:8/f0));
xdc_excitation (emit_aperture, excitation);

%% ��������
fieldWidth    = 15;                                   % ������������ˮƽ�������,����
fieldLength   = 15;%round((fieldWidth * (array_size / repeatNumber)) / (z_size / N));
fieldRadius   = floor(fieldWidth / 2);                % ������������ˮƽ����뾶
pha_Pos       = zeros(1, 3);                          % ���������ĵ�λ��
depthCenter   = round(focal_depth / (z_size / N));    % ����������������ȷ������ĵ�
horCenter     = round(lineNumber / 2);                % ��������������ˮƽ�������ĵ�
fp            = zeros(fieldLength, fieldWidth);       % �������� 

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

