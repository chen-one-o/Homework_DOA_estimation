% Capon algorithm for DOA estimation
clc
clear all
close all
sensor_number=10; %阵元数
source_number=2; %信源数
T=200;  %阵列采样快拍数
sig1=exp(j*(0.1*pi*(0:T-1))).';
sig2=exp(j*(0.2*pi*(0:T-1))).';  %窄带复正弦信号
S0=[sig1.';sig2.'];  %信号矩阵形式
snr=10; %信噪比dB
source_doa1=-45*pi/180;
source_doa2=60*pi/180;
source_doa=[source_doa1 source_doa2]; %信号入射角度
A=zeros(sensor_number,source_number);
A=[(exp(-j*pi*(0:sensor_number-1)*sin(source_doa(1)))).' (exp(-j*pi*(0:sensor_number-1)*sin(source_doa(2)))).'];  %方向矩阵
%%%高斯白噪声*********************************************************
n=zeros(sensor_number,T); 
real_noise0=randn(sensor_number,T);
imag_noise0=randn(sensor_number,T);
mean_real_noise=mean(mean(real_noise0));
mean_imag_noise=mean(mean(imag_noise0));
real_noise=real_noise0-mean_real_noise;
imag_noise=imag_noise0-mean_imag_noise;
noise =(real_noise+j*imag_noise)/(2^0.5);
%*******************************************************************
PWn1=noise*noise'/T;
PWn2=diag(PWn1);
PWn=mean(PWn2); % 噪声功率
PWs1=S0*S0'/T;
PWs2=diag(PWs1);
PWs=mean(PWs2);   %信号功率
aa=(PWn*(10^(snr/10))/PWs)^0.5;    %信号相对幅值
S=aa*S0;
x=A*S+noise;
R=x*x'/T;   %阵列协方差矩阵
%**********谱峰搜索过程******************************************
search_doa=-90:90;
for i=1:length(search_doa)
    a=exp(-j*pi*(0:sensor_number-1)*sin(search_doa(i)*pi/180)).';
    Pcapon(i)=1/abs(a'*pinv(R)*a);
end
%**************************************************************

plot(search_doa,10*log(Pcapon),'r');
xlabel('入射角度/dB');
ylabel('空间谱');
legend('MP Spectrum');
title('MP算法');
hold on
grid on


