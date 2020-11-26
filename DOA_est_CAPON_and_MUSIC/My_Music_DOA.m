%% MUSIC algorithm for DOA estimation
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

%% MUSIC算法步骤
% Step1：根据N个接收信号矢量得到协方差矩阵的估计值；
S = randn(doas,snapshot_num);
X0 = A * S;
X  = awgn(X0,snr,'measured');
Rxx = X*X'/snapshot_num;

% Step2：对Rxx进行特征值分解；按特征值从大到小排序，...
%        把与信号个数 K 相等的最大特征值对应的特征向量看作信号子空间，...
%        剩下 (M-K) 个特征值对应的特征向量看作噪声空间，则Rxx写作两部分和。
InvS = inv(Rxx);
[EVector,EValue] = eig(Rxx); % Rxx*EVector = EVector*EValue (8*8)
EVA = diag(EValue);
[EVA,I] = sort(EVA);
EVA=fliplr(EVA);
EVector=fliplr(EVector(:,I));

% Step3：使θ变化，按照P_music(θ)=1/a^H(θ)*U_N*U_N^H*a(θ)计算谱函数，...
%        通过寻找峰值来得到波达方向的估计值。
for sch_ang = 1:361  %搜索范围-90°至90°
    angle(sch_ang) = (sch_ang-181)/2;
    phi_angle = angle(sch_ang)*pi/180;
    a = exp(-1i*2*pi*element_position*sin(phi_angle)).';
    L = doas;
    EN = EVector(:,(L+1):elements);
    Spec(sch_ang) = 1/(a'*EN*EN'*a);
end

%% 绘图
Spec = abs(Spec);
Spec_max = max(Spec);
Spec = 10*log10(Spec/Spec_max);
plot(angle,Spec,'Linewidth',2);
grid on
xlabel('角度 (°)'),ylabel('空间谱 (dB)');
title('Normal MUSIC for DOA Estiamtion');
axis([-90 90 -45 0]),set(gca, 'XTick',-90:30:90);


