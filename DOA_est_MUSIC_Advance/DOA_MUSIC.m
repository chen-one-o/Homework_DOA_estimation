%% Music算法估计阵列信号波达方向DOA
clc,clear,close all;

%% 均匀阵列，阵元数=8，半波长，角度分别为10°，20°，30°
wavelength = 1; 
d = wavelength / 2; 
design_ula = design_array_1d('ula', 8, d);
doas = linspace(pi/18, pi/6, 3);
power_source = 1;
power_noise = 1;
snapshot_count = 1000;
source_count = length(doas);

%% 
% 随机模型
[~, R] = snapshot_gen_sto(design_ula, doas, wavelength, snapshot_count, power_noise, power_source);

% 阵源数目探测
[~, l] = eig(0.5*(R+R'), 'vector');
l = flipud(l);
n_mdl = sn_mdl(l, design_ula.element_count, snapshot_count);
n_aic = sn_aic(l, design_ula.element_count, snapshot_count);
fprintf('There are %d sources.\n', source_count);
fprintf('# of sources estimated by MDL = %d\n', n_mdl);
fprintf('# of sources estimated by AIC = %d\n', n_aic);

%% normal MUSIC
tic;
sp_normal2 = music_1d(R, source_count, design_ula, wavelength, 180, 'RefineEstimates', true);
toc;
sp_normal2.true_positions = doas;
fprintf('[MUSIC] Estimated DOAs:\n');
disp(sp_normal2.x_est);
plot_sp(sp_normal2, 'title', 'Normal MUSIC','PlotType','cartesian');


