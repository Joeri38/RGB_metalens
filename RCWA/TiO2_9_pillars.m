clc;clear all;

%addpath('./reticolo_allege');
addpath('reticolo_allege');

input = struct;
input.angle_delta = 0;
input.angle_theta = 0;
input.substrate_material = 'SiO2';
input.superstrate_material = 'air';
input.num_fourier = 12; % OG: 12
input.material_inside_circle = 'TiO2_ALD';
input.material_outside_circle = 'air';
input.cent_x = 0;
input.cent_y = 0;
input.patterned_layer_t_range = 0.8e-6; % [m] membrane thickness
input.period_range = 0.75e-6; % [m] OG: 0.64e-6
input.staircase_N = 20;

res = 2;
input.diameter_1_range = linspace(0.05e-6, 0.25e-6, res);
input.diameter_2_range = linspace(0.05e-6, 0.25e-6, res);
input.diameter_3_range = linspace(0.05e-6, 0.25e-6, res);
input.diameter_4_range = linspace(0.05e-6, 0.25e-6, res);
input.diameter_5_range = linspace(0.05e-6, 0.25e-6, res);
input.diameter_6_range = linspace(0.05e-6, 0.25e-6, res);
input.diameter_7_range = linspace(0.05e-6, 0.25e-6, res);
input.diameter_8_range = linspace(0.05e-6, 0.25e-6, res);
input.diameter_9_range = linspace(0.05e-6, 0.25e-6, res);
%input.diameter_1_range = 0.05e-6;
%input.diameter_2_range = 0.10e-6;
%input.diameter_3_range = 0.10e-6;
%input.diameter_4_range = 0.15e-6;
%input.diameter_5_range = 0.17e-6;
%input.diameter_6_range = 0.19e-6;
%input.diameter_7_range = 0.20e-6;
%input.diameter_8_range = 0.21e-6;
%input.diameter_9_range = 0.23e-6;
input.lambda_range = 0.658e-6; % [m] 3.5 um
input.lambda0 = 0.658e-6; % [m] 3.5 um
input.debug = true;
input.multiplex = true;

total_tic = tic;
[all_meta_atoms,completed_jobs] = get_9_square_meta_atom_GD(input);

% to plot supercell: 
%get_supercell_meta_atom_GD.m %line 466: parfor
%get_supercell_meta_atom_GD.m %line 523: settings.debug = true;

%%
diameter_1 = transpose(input.diameter_1_range);
phase = mod(angle(completed_jobs.amplitude_TE),2*pi); % [rad]
transmission = abs(completed_jobs.amplitude_TE).^2; % [1]
phase_TM = mod(angle(completed_jobs.amplitude_TM),2*pi); % [rad]
transmission_TM = abs(completed_jobs.amplitude_TM).^2; % [1]
z_real = real(completed_jobs.amplitude_TE);
z_imag = imag(completed_jobs.amplitude_TE);
n_eff = completed_jobs.max_real_n_eff;

save_data = true;
if save_data
    save("library/TiO2 pillars/z_real_" + num2str(input.lambda0*1e6) + "um_p=" + num2str(input.period_range*1e9) + "nm_t=" + num2str(input.patterned_layer_t_range*1e9) + "nm_res=" + num2str(res) + ".mat", 'z_real');
    save("library/TiO2 pillars/z_imag_" + num2str(input.lambda0*1e6) + "um_p=" + num2str(input.period_range*1e9) + "nm_t=" + num2str(input.patterned_layer_t_range*1e9) + "nm_res=" + num2str(res) + ".mat", 'z_imag');
end

total_time = toc(total_tic);
fprintf('Total time elapsed: %g s, time per job: %g s\n', total_time);

figure;

subplot(1,2,1);
plot(input.diameter_1_range*1e6,phase(1:res), 'o');
%plot(input.diameter_range*1e6,z_real(1:res), 'o');
xlabel('Diameter 1 [um]');
ylabel('Phase');
%ylabel('Real(z)');
title("Lambda " + num2str(input.lambda0*1e6) + "um, period " + num2str(input.period_range*1e9) + "nm, thickness " + num2str(input.patterned_layer_t_range*1e9) + "nm");
%ylim([0, 2.05]);

subplot(1,2,2);
plot(input.diameter_1_range*1e6,transmission(1:res).^2, 'o');
%plot(input.diameter_range*1e6,z_imag(1:res), 'o');
xlabel('Diameter 1 [um]');
ylabel('Transmission');
%ylabel('Imag(z)');
title("Lambda " + num2str(input.lambda0*1e6) + "um, period " + num2str(input.period_range*1e9) + "nm, thickness " + num2str(input.patterned_layer_t_range*1e9) + "nm");
%ylim([0, 105]);

saveas(gcf, "figs/TiO2 9 pillars/lamb=" + num2str(input.lambda0*1e6) + "um_p=" + num2str(input.period_range*1e9) + "nm_t=" + num2str(input.patterned_layer_t_range*1e9) + "nm.fig");