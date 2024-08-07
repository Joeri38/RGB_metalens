% This is a sample script to demonstrate the function of psf2mtf.
% Date created: 20200613
% Author: Daniel Lim
clear all;

image_filename = '0um_nobackground.tif';
% image_filename = '0um.tiff';

focal_plane_intensity = im2double(imread(image_filename));

pixel_size = 9.08e-6; % [m] for the Optikos MTF
magnification = 40; % using the 40x objective on the Optikos MTF

padding = 0; % [pixels], number of additional zeros to pad around the focal spot
             % window to increase frequency resolution of the MTF

[N_vert,N_horz] = size(focal_plane_intensity);

% we need to identify the pixel with the maximum intensity to be the center
% of the PSF
[max_row,max_col] = find(focal_plane_intensity == max(focal_plane_intensity(:)));
max_row = max_row(1);
max_col = max_col(1);

% Also define a window around the pixel with the maximum intensity
window_size = 200; 

% Plot the input values
figure('position',[0,0,1000,1000]);
subplot(2,2,1);
surf(focal_plane_intensity);
shading flat;
view([0,90]);
box on;
hold on
% plot3(max_col,max_row,1.1*max(focal_plane_intensity(:)),'rx','markersize',10);
plot3([max_col - window_size/2,max_col - window_size/2,max_col + window_size/2, max_col + window_size/2, max_col - window_size/2],...
      [max_row - window_size/2,max_row + window_size/2, max_row + window_size/2, max_row - window_size/2, max_row - window_size/2],...
      ones(5,1)*1.1*max(focal_plane_intensity(:)),'w-','linewidth',2);
hold off
colormap gray
axis tight equal
xlabel('X [pixel]');
ylabel('Y [pixel]');
set(gca,'Fontsize',16);
set(gca,'Linewidth',2);
set(gca,'TickDir','out');
cb = colorbar;
title('Focal plane intensity');

%% Also calculate the theoretical MTF for the diffraction limit
NA_system = 0.1; 
lambda0 = 0.633e-6; % [m]
lines_per_mm_range = linspace(1,350,1000); % [lp/mm]
phi_range = acos(lambda0*1e3*lines_per_mm_range/(2*NA_system)); % note we change the units of lambda0 to [mm]
MTF_diff_limit = (2/pi)*(phi_range - cos(phi_range).*sin(phi_range));
MTF_diff_limit(imag(MTF_diff_limit) ~= 0) = 0;

%% Load the inputs into the PSF2MTF function and get the results

input = struct;
input.focal_plane_intensity = focal_plane_intensity;
input.horz_center = max_col;
input.vert_center = max_row;
input.window_size = window_size;
input.padding = padding;
input.pixel_size = pixel_size/magnification;
nu_cutoff = 2*NA_system/(lambda0*1e3);

output = psf2mtf(input);

%% Plot the results
subplot(2,2,2);
plot(output.horz_lp_per_mm,output.horz_MTF,'-rx');
hold on
plot(output.vert_lp_per_mm,output.vert_MTF,'-bx');
plot(lines_per_mm_range,MTF_diff_limit,'k--');
hold off
set(gca,'Fontsize',16);
set(gca,'Linewidth',2);
xlabel('Line pairs per mm');
ylabel('MTF');
legend({'Horizontal MTF','Vertical MTF','Diffraction limit'});
legend boxoff
xlim([0,nu_cutoff*1.2]);
ylim([0,1]);
title('MTF from PSF');

subplot(2,2,3);
surf(output.MTF_2D_fx,output.MTF_2D_fy,output.MTF_2D);
shading flat;
view([0,90]);
xlabel('x line pairs per mm');
ylabel('y line pairs per mm');
cb = colorbar;
set(gca,'Fontsize',16);
set(gca,'Linewidth',2);
title('2D MTF');
set(gca,'TickDir','out');
axis tight equal
xlim([-nu_cutoff*1.1,nu_cutoff*1.2]);
ylim([-nu_cutoff*1.1,nu_cutoff*1.2]);