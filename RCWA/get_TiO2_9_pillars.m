addpath('reticolo_allege');

disp('Hello from MATLAB, size of designs is:');

% Load designs to simulate
load('../designs.mat');
display(size(designs));

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

input.diameters = designs;
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
phase = mod(angle(completed_jobs.amplitude_TE),2*pi); % [rad]
transmission = abs(completed_jobs.amplitude_TE).^2; % [1]
phase_TM = mod(angle(completed_jobs.amplitude_TM),2*pi); % [rad]
transmission_TM = abs(completed_jobs.amplitude_TM).^2; % [1]
z_real = real(completed_jobs.amplitude_TE);
z_imag = imag(completed_jobs.amplitude_TE);
n_eff = completed_jobs.max_real_n_eff;

save_data = true;
if save_data
    save("../z_real.mat", 'z_real');
    save("../z_imag.mat", 'z_imag');
end

total_time = toc(total_tic);
fprintf('Total time elapsed: %g s\n', total_time);