clc;clear all;

addpath('./reticolo_allege');

input = struct;
input.angle_delta = 0;
input.angle_theta = 0;
input.substrate_material = 'air';
input.superstrate_material = 'air';
input.num_fourier = 10;
input.material_inside_circle = 'Si';
input.material_outside_circle = 'air';
input.cent_x = 0;
input.cent_y = 0;
input.patterned_layer_t_range = 5e-6; % [m]
input.period_range = 0.7e-6; % [m]
input.staircase_N = 10; 
input.diameter_range = linspace(0,input.period_range,100); % [m]
input.lambda_range = 0.76e-6; % [m]
input.lambda0 = 0.76e-6; % [m]
input.debug = false;
input.multiplex = true;

[all_meta_atoms,completed_jobs] = get_circular_meta_atom_GD(input);


%% Plot the amplitude obtained

phases = angle(completed_jobs.amplitude_TE); % [rad]
transmission = abs(completed_jobs.amplitude_TE).^2; % [1]

figure;

subplot(1,2,1);
plot(input.diameter_range*1e6,phases);
xlabel('Diameter [um]');
ylabel('Phase [rad]');
title('Small library phase');


subplot(1,2,2);
plot(input.diameter_range*1e6,100*transmission);
xlabel('Diameter [um]');
ylabel('Transmission [%]');
title('Small library transmission');

%% Estimate the amplitude based on the dominant bloch mode
neff = completed_jobs.max_real_n_eff;
k0L = (2*pi/input.lambda0)*input.patterned_layer_t_range;
n_substrate = get_refractive_index(input.substrate_material,input.lambda0*1e6);
n_superstrate = get_refractive_index(input.superstrate_material,input.lambda0*1e6);
est_transmission = trans_fn(neff,k0L,n_substrate,n_superstrate);
est_phases = angle(est_transmission); % [rad]
est_transmission = abs(est_transmission).^2; % [1]

figure('units','normalized','position',[0,0,1,1]);

subplot(1,2,1);
plot(input.diameter_range*1e6,est_phases,'-','linewidth',2);
xlabel('Diameter [um]');
ylabel('Phase [rad]');
title('Small library phase (dominant Bloch mode est.)');
set(gca,'Fontsize',16);
set(gca,'Linewidth',2);
box on;


subplot(1,2,2);
plot(input.diameter_range*1e6,100*est_transmission,'-','linewidth',2);
xlabel('Diameter [um]');
ylabel('Transmission [%]');
title('Small library transmission (dominant Bloch mode est.)');
set(gca,'Fontsize',16);
set(gca,'Linewidth',2);
box on;

%% Auxillary functions
function transmission = trans_fn(n1_eff,k0L,n_substrate,n_superstrate)
% Gives the complex transmission coefficient for a medium 
% of effective index n1_eff, vacuum optical path lens
% k0L, on a substrate of index n_substrate, and below
% a superstrate of index n_superstrate.

% params:
% -n1_eff [1]: effective index of patterned layer. can be vector.
% -k0L [rad]: vacuum optical path length of patterned layer. must be either
% a single value or the same size as n1_eff.
% -n_substrate [1]: substrate index. must be a single value 
% or the same size as n1_eff.
% -n_superstrate [1]: superstrate index. must be a single value 
% or the same size as n1_eff.

% returns:
% -transmission [1]: complex transmission amplitude in the superstrate so
% that |transmission|^2 is the efficiency.

size_neff = size(n1_eff);
assert(length(k0L)==1 || all(size_neff == size(k0L)),'k0L must be same size as n1_eff');
assert(length(n_substrate)==1 || all(size_neff == size(n_substrate)),...
    'n_substrate must be same size as n1_eff');
assert(length(n_superstrate)==1 ||all(size_neff == size(n_superstrate)),...
    'n_superstrate must be same size as n1_eff');
c_t12 = 2*n_substrate./(n_substrate + n1_eff);
c_t23 = 2*n1_eff./(n1_eff + n_superstrate);
c_r21 = (n_substrate - n1_eff)./(n_substrate + n1_eff);
c_r23 = (n_superstrate - n1_eff)./(n_superstrate + n1_eff);
delta = n1_eff.*k0L;
transmission = sqrt(n_superstrate./n_substrate).*c_t12.*c_t23.*exp(1i*delta)./(1-c_r21.*c_r23.*exp(2i*delta));
end