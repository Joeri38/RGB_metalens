clc;clear all;

addpath('./reticolo_allege');
addpath('reticolo_allege');

input = struct;
input.angle_delta = 0;
input.angle_theta = 0;
input.substrate_material = 'air';
input.superstrate_material = 'air';
input.num_fourier = 12;
input.material_inside_circle = 'air';
input.material_outside_circle = 'Si';
input.cent_x = 0;
input.cent_y = 0;
input.patterned_layer_t_range = 1.5e-6; % [m] membrane thickness
input.period_range = 0.320e-6; % [m]
input.staircase_N = 20; 
input.diameter_range = linspace(0.1e-6,0.3e-6,100); % [m] range of diameters from 0 to 2.5e-6 in 100 steps
input.lambda_range = 0.76e-6; % [m]
input.lambda0 = 0.76e-6; % [m]
input.debug = false;
input.multiplex = true;

[all_meta_atoms,completed_jobs] = get_circular_meta_atom_GD(input);


%%

phases = angle(completed_jobs.amplitude_TE); % [rad]
transmission = abs(completed_jobs.amplitude_TE).^2; % [1]

figure;

subplot(1,2,1);
plot(input.diameter_range*1e6,phases/pi);
xlabel('Diameter [um]');
ylabel('Phase [rad]');
title('Small library phase');


subplot(1,2,2);
plot(input.diameter_range*1e6,100*transmission);
xlabel('Diameter [um]');
ylabel('Transmission [%]');
title('Small library transmission');
