%% CAPON algorithm for DOA estimation
clc,clear,close all;

%% �������У���Ԫ��=8���벨�����Ƕȷֱ�Ϊ10�㣬20�㣬30��
wavelength = 1;       % ��λ����
d = wavelength / 2;   % �벨��
elements = 8;         % ��������
element_position = 0:d:(elements-1)*d; % ��Ԫλ��
doas = 3;             % ���﷽������
theta = [10 20 30];   % �Ƕ�
% snr = 15;            
snr = 30;             % ���������
snapshot_num = 500;   % ������
% �������
A = exp(-1i*2*pi*element_position.'*sin(theta*pi/180));

%% CAPON�㷨����
% Step1������N�������ź�ʸ���õ�Э�������Ĺ���ֵ��
S = randn(doas,snapshot_num);
X0 = A * S;
X  = awgn(X0,snr,'measured');
Rxx = X*X'/snapshot_num;

% Step2��ͨ����P_Capon1/a^H(��)*Rxx^(-1)*a(��)�����׺�����...
%        ͨ��Ѱ�ҷ�ֵ���õ����﷽��Ĺ���ֵ��
for sch_ang = 1:361  %������Χ-90����90��
    angle(sch_ang) = (sch_ang-181)/2;
    phi_angle = angle(sch_ang)*pi/180;
    a = exp(-1i*2*pi*element_position*sin(phi_angle)).';
    Spec(sch_ang) = 1./(a'*pinv(Rxx)*a);
end

%% ��ͼ
Spec = abs(Spec);
Spec_max = max(Spec);
Spec = 10*log10(Spec/Spec_max);
plot(angle,Spec,'Linewidth',2);
grid on
xlabel('�Ƕ� (��)'),ylabel('�ռ��� (dB)');
title('CAPON for DOA Estimation');
axis([-90 90 -25 0]),set(gca, 'XTick',-90:30:90);

