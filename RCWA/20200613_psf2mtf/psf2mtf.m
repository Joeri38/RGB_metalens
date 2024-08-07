function output = psf2mtf(input)
% PSF2MTF converts 1D or 2D focal intensity to an MTF
% Usage, with a 1D radial plot of the focal plane intensity: 
%
% Args:
% 	input: structure array with fields:
%       focal_plane_intensity [a.u.]: 2D matrix of size (N_vert,N_horz) for the
%           focal plane intensity pattern. This must be a real-valued
%           matrix and sampled with square pixels. 
%       pixel_size [m]: horizontal and vertical pixel size for
%           focal_plane_intensity
%       horz_center [pixel]: horizontal position of the focal point center
%       vert_center [pixel]: vertical position of the focal point center
%       window_size [pixel]: size of the window around the focal point to
%           consider. the full window is window_size x window_size. 
%       padding [1]: number of pixels of zero-padding to add around the
%           selected window. this increases the frequency resolution of the
%           MTF.
% 
% Returns:
%   output: structure array with fields:
%       horz_lp_per_mm [1/mm]: x-axis for the horizontal (tangential) MTF
%       vert_lp_per_mm [1/mm]: x-axis for the vertical (sagittal) MTF
%       horz_MTF [1]: horizontal (tangential) MTF
%       vert_MTF [1]: vertical (sagittal) MTF
%       MTF_2D [1]: two-dimensional MTF values
%       MTF_2D_fx [1/mm]: x-axis for MTF_2D
%       MTF_2D_fy [1/mm]: y-axis for MTF_2D

%% Do assertion checks
assert(isstruct(input),...
    sprintf("%s: input must be a structure array. It is currently a %s",...
    mfilename,class(input)));

% Verify that the input contains the required fields
field_names = {'focal_plane_intensity','pixel_size',...
               'horz_center','vert_center',...
               'window_size','padding'};
check_field_exists = isfield(input,field_names);
for i = 1:length(check_field_exists)
    assert(check_field_exists(i),...
        sprintf("%s: input is missing field %s",mfilename,field_names{i}));
end

focal_plane_intensity = input.focal_plane_intensity;

pixel_size = input.pixel_size;
window_size = input.window_size;
padding = input.padding;

[N_vert,N_horz] = size(focal_plane_intensity);

% ensure that the intensity matrix is strictly real
assert(isreal(focal_plane_intensity),...
    sprintf("%s: focal_plane_intensity is not strictly real",mfilename));

horz_center = input.horz_center;
vert_center = input.vert_center;
assert(horz_center > 1 && horz_center < N_horz,...
    sprintf("%s: horz_center is not within the bounds",mfilename));
assert(vert_center > 1 && vert_center < N_vert,...
    sprintf("%s: vert_center is not within the bounds",mfilename));

%% Window the focal plane intensity to be centered around (horz_center,vert_center)
% window_size = min([N_vert-vert_center,vert_center,N_horz-horz_center,horz_center]);
% half_window_size = floor(window_size/2);
half_window_size = floor(window_size/2);

focal_plane_intensity_windowed = focal_plane_intensity((vert_center-half_window_size):(vert_center+half_window_size),...
                                                       (horz_center-half_window_size):(horz_center+half_window_size));
% focal_plane_intensity_windowed = focal_plane_intensity_windowed/sum(focal_plane_intensity_windowed(:));

N_windowed = 2*half_window_size+1;

%% Add the padding
focal_plane_intensity_padded = zeros(N_windowed + 2*padding,N_windowed + 2*padding);
center = padding + half_window_size + 1;
focal_plane_intensity_padded((center-half_window_size):(center+half_window_size),(center-half_window_size):(center+half_window_size)) = focal_plane_intensity_windowed;

new_N = 2*(half_window_size + padding) + 1;
%% Perform the Fourier transforms
df_horz = 1/(new_N*pixel_size*1e3); % [1/mm], incremental change in frequency in horizontal direction
df_vert = 1/(new_N*pixel_size*1e3); % [1/mm], incremental change in frequency in vertical direction

MTF_2D_unshifted = abs(fft2((focal_plane_intensity_padded)));

% extract the horizontal and vertical components
horz_MTF = MTF_2D_unshifted(1,:);
vert_MTF = MTF_2D_unshifted(:,1);

MTF_2D_shifted = fftshift(MTF_2D_unshifted);

horz_lp_per_mm_full = (0:(new_N-1))*df_horz;
vert_lp_per_mm_full = (0:(new_N-1))*df_vert;

% only return frequencies up to the Nyquist frequency for the camera pixel
f_nyquist = 1/(2*pixel_size*1e3);

horz_MTF = horz_MTF(horz_lp_per_mm_full < f_nyquist);
horz_lp_per_mm = horz_lp_per_mm_full(horz_lp_per_mm_full < f_nyquist);
vert_MTF = vert_MTF(vert_lp_per_mm_full < f_nyquist);
vert_lp_per_mm = vert_lp_per_mm_full(vert_lp_per_mm_full < f_nyquist);

output = struct;
output.horz_MTF = horz_MTF/horz_MTF(1);
output.vert_MTF = vert_MTF/vert_MTF(1); 
output.horz_lp_per_mm = horz_lp_per_mm;
output.vert_lp_per_mm = vert_lp_per_mm;
output.MTF_2D = MTF_2D_shifted/MTF_2D_unshifted(1,1);
output.MTF_2D_fx = (-new_N/2:1:(new_N/2-1))*df_horz;
output.MTF_2D_fy = (-new_N/2:1:(new_N/2-1))*df_vert;
end

