%% CAPON algorithm for DOA estimation
clc,clear,close all;

%% 均匀阵列，阵元数=8，半波长，角度分别为10°，20°，30°
wavelength = 1;       % 单位波长
d = wavelength / 2;   % 半波长
elements = 8;         % 阵列数量
element_position = 0:d:(elements-1)*d; % 阵元位置
doas = 3;             % 波达方向数量
theta = [10 20 30];   % 角度
% snr = 15;            
snr = 30;             % 输入信噪比
snapshot_num = 500;   % 快拍数
% 导向矩阵
A = exp(-1i*2*pi*element_position.'*sin(theta*pi/180));

%% CAPON算法步骤
% Step1：根据N个接收信号矢量得到协方差矩阵的估计值；
S = randn(doas,snapshot_num);
X0 = A * S;
X  = awgn(X0,snr,'measured');
Rxx = X*X'/snapshot_num;

% Step2：通过对P_Capon1/a^H(θ)*Rxx^(-1)*a(θ)计算谱函数，...
%        通过寻找峰值来得到波达方向的估计值。
for sch_ang = 1:361  %搜索范围-90°至90°
    angle(sch_ang) = (sch_ang-181)/2;
    phi_angle = angle(sch_ang)*pi/180;
    a = exp(-1i*2*pi*element_position*sin(phi_angle)).';
    Spec(sch_ang) = 1./(a'*pinv(Rxx)*a);
end

%% 绘图
Spec = abs(Spec);
Spec_max = max(Spec);
Spec = 10*log10(Spec/Spec_max);
plot(angle,Spec,'Linewidth',2);
grid on
xlabel('角度 (°)'),ylabel('空间谱 (dB)');
title('CAPON for DOA Estimation');
axis([-90 90 -25 0]),set(gca, 'XTick',-90:30:90);

