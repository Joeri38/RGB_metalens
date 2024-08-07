% Load phases
load('library/1.6um_320nm_2.5um.mat');
phase_16 = phases;
load('library/2.0um_320nm_2.5um.mat');
phase_20 = phases;
load('library/4.8um_320nm_2.5um.mat');
phase_48 = phases;

% Plot phases
figure;
color = linspace(0,0.5,100);
scatter3(phase_16/pi, phase_20/pi, phase_48/pi, 50, color);
colorbar;
xlabel('Phase 1.6um');
ylabel('Phase 2.0um');
zlabel('Phase 4.8um');

% Load neff
load('library/neff_1.6um_320nm_2.5um.mat');
neff_16 = n_eff;
load('library/neff_2.0um_320nm_2.5um.mat');
neff_20 = n_eff;
load('library/neff_4.8um_320nm_2.5um.mat');
neff_48 = n_eff;

%Plot neff
figure;
plot(linspace(0,0.5,100), neff_16);
hold on;
plot(linspace(0,0.5,100), neff_20);
plot(linspace(0,0.5,100), neff_48);
hold off;
legend('1.6um', '2.0um', '4.8um');