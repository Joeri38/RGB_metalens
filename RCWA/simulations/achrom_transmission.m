function out = achrom_transmission(settings)
% Reticolo function to compute the complex zeroth order transmission
% coefficient through a hole or pillar array.

% The following geometries are available:
%   -hollow_square
%   -cross_mid_circle
%   -hollow_circle 
%   -two_hollow_squares
%   -two_hollow_circles
%   -supported_hollow_square
%   -four_squares
%   -square_cog

% Version 0: updated 20200725 (all SI units)
% Version 1: updated 20210622 (added supported_hollow_square, four_squares
%            and square_cog)

% settings: struct with the following keys:
%   -wavelength [m]: design vacuum wavelength
%   -angle_delta [deg]: rotation angle of the XZ plane of incidence 
%                       from the x-axis
%   -angle_theta [deg]: incidence angle from the normal (z-axis)
%   -n_incident_medium [1]: refractive index, incident medium
%   -n_substrate_medium [1]: refractive index, substrate medium
%   -n_superstrate_medium [1]: refractive index, superstrate medium
%   -n_transmitted_medium [1]: refractive index, transmitted medium
%   -num_fourier [1]: number of Fourier components in x and y directions
%                     each
%   -shape [string]: "hollow_square", "cross_mid_circle", "hollow_circle"
%                   "two_hollow_squares","two_hollow_circles"
%   -pillar_cent_x [m]: pillar center x-position
%   -pillar_cent_y [m]: pillar center y-position
%   -pillar_layer_t [m]: thickness of the meta-element
%   -d1 to d10 [m] as defined in 20200724_achromatic_holens_metaelements.pptx
%   -nbackground [1]: refractive index of Silicon layer at target
%                      wavelength
%   -ninclusion [1]: refractive index of hole in silicon layer
%   -period [m]: center-to-center spacing in x and y directions
%   -staircase_N [1]: Number of points for staircase approximation
%   -incident_layer_t [m]: thickness of incident layer
%   -substrate_layer_t [m]: thickness of substrate layer
%   -superstrate_layer_t [m]: thickness of superstrate layer
%   -transmitted_layer_t [m]: thickness of transmitted layer
%   -debug: either 0 or 1. automatically assumes 0 if not set. 

% Returns:
%   -out: struct with:
%       -amplitude_TM: complex TM transmission amplitude
%       -amplitude_TE: complex TE transmission amplitude
%       -n_eff: effective index for membrane layer

if isfield(settings,'debug')
    debug =  settings.debug;
else 
    debug = 0;
end
wavelength = settings.wavelength*1e6; % [m] to [um]
angle_delta = settings.angle_delta; % [deg], rotation angle of the XZ plane of incidence 
                 % from the x-axis
angle_theta = settings.angle_theta; % [deg], incidence angle from the normal (z-axis)

n_incident_medium = settings.n_incident_medium; % [1]
n_substrate_medium = settings.n_substrate_medium; % [1]
n_superstrate_medium = settings.n_superstrate_medium; % [1]
n_transmitted_medium = settings.n_transmitted_medium; % [1]
k_parallel = n_incident_medium*sin(angle_theta*pi/180);
nn = settings.num_fourier*[1,1]; % number of Fourier harmonics in the x and y directions

% set up the textures for the pillar array
nbackground = settings.nbackground; % [1], refractive index of Silicon at 1.55 um (using 
                      % Li 1980: n 1.2-14 um; 293 K on refractiveindex.info
ninclusion = settings.ninclusion; % fill the pillar with air                      
pillar_cent_x = settings.pillar_cent_x*1e6; % [m] to [um]
pillar_cent_y = settings.pillar_cent_y*1e6; % [m] to [um]
period = [1, 1]*settings.period*1e6; % [m] to [um], x direction, y direction
textures = {};
textures{1} = {n_incident_medium};
textures{2} = {n_substrate_medium};

% get the geometrical details and load it into textures(3)
shape = lower(settings.shape); % string
switch shape
    case 'hollow_square'
        d1 = settings.d1*1e6; % [um]
        d2 = settings.d2*1e6; % [um]
        d3 = settings.d3*1e6; % [um]
        d4 = settings.d4*1e6; % [um]
        assert(d3 < d1,"D3 must be smaller than D1");
        assert(d4 < d1,"D4 must be smaller than D2");
        assert(d1 > 0 && d2 > 0,"D1 and D2 must be positive");
        if (d3 > 0 && d4 > 0)
            % then we need to add a hollow core
            textures{3} = {nbackground,...
                            [pillar_cent_x,pillar_cent_y,...
                                d1,d2,...
                                ninclusion,1],...
                            [pillar_cent_x,pillar_cent_y,...
                                d3,d4,...
                                nbackground,1]};
        else
            % no hollow core needed
            textures{3} = {nbackground,[pillar_cent_x,pillar_cent_y,...
                            d1,d2,...
                            ninclusion,1]};
        end
    case 'two_hollow_squares'
        d1 = settings.d1*1e6; % [um]
        d2 = settings.d2*1e6; % [um]
        d3 = settings.d3*1e6; % [um]
        d4 = settings.d4*1e6; % [um]
        d5 = settings.d5*1e6; % [um]
        d6 = settings.d6*1e6; % [um]
        d7 = settings.d7*1e6; % [um]
        d8 = settings.d8*1e6; % [um]
        assert(d3 < d1,"D3 must be smaller than D1");
        assert(d4 < d1,"D4 must be smaller than D2");
        assert(d1 > 0 && d2 > 0,"D1 and D2 must be positive");
        assert(d5 < d3,"D5 must be smaller than D3");
        assert(d6 < d4,"D6 must be smaller than D4");
        assert(d7 < d5,"D7 must be smaller than D5");
        assert(d8 < d6,"D8 must be smaller than D6");
        assert(d7 > 0 && d8 > 0,"You must have two squares present");
        textures{3} = {nbackground,...
                        [pillar_cent_x,pillar_cent_y,...
                            d1,d2,...
                            ninclusion,1],...
                        [pillar_cent_x,pillar_cent_y,...
                            d3,d4,...
                            nbackground,1],...
                        [pillar_cent_x,pillar_cent_y,...
                            d5,d6,...
                            ninclusion,1],...
                        [pillar_cent_x,pillar_cent_y,...
                            d7,d8,...
                            nbackground,1]};
    case 'four_squares'
        d1 = settings.d1*1e6; % [um]
        d2 = settings.d2*1e6; % [um]
        d3 = settings.d3*1e6; % [um]
        d4 = settings.d4*1e6; % [um]
        assert(d3 < d1,"D3 must be smaller than D1");
        assert(d4 < d1,"D4 must be smaller than D2");
        assert(d1 > 0 && d2 > 0,"D1 and D2 must be positive");
        assert(d3 >= 0 && d4 >= 0,"D3 and D4 must be non-negative");
        cen_x = (d1+d3)/4;
        cen_y = (d2+d4)/4;
        L_x = (d1-d3)/2;
        L_y = (d2-d4)/2;
        textures{3} = {nbackground,...
                        [-cen_x,-cen_y,...
                            L_x,L_y,...
                            ninclusion,1],...
                        [cen_x,-cen_y,...
                            L_x,L_y,...
                            ninclusion,1],...
                        [-cen_x,cen_y,...
                            L_x,L_y,...
                            ninclusion,1],...
                        [cen_x,cen_y,...
                            L_x,L_y,...
                            ninclusion,1]};  
    case 'square_cog'
        d1 = settings.d1*1e6; % [um]
        d2 = settings.d2*1e6; % [um]
        d3 = settings.d3*1e6; % [um]
        d4 = settings.d4*1e6; % [um]
        d5 = settings.d5*1e6; % [um]
        d6 = settings.d6*1e6; % [um]
        assert(d3 < d1,"D3 must be smaller than D1");
        assert(d4 < d1,"D4 must be smaller than D2");
        assert(d1 > 0 && d2 > 0,"D1 and D2 must be positive");
        assert(d3 >= 0 && d4 >= 0,"D3 and D4 must be non-negative");
        assert(d5 < d3,"D5 must be smaller than D3");
        assert(d6 < d4,"D6 must be smaller than D4");

        textures{3} = {nbackground,...
                        [-((d1-d5)/4+d5/2),((d2-d6)/4+d6/2),...
                          (d1-d5)/2,(d2-d6)/2,...
                            ninclusion,1],...
                        [((d1-d5)/4+d5/2),((d2-d6)/4+d6/2),...
                          (d1-d5)/2,(d2-d6)/2,...
                            ninclusion,1],...
                        [-((d1-d5)/4+d5/2),-((d2-d6)/4+d6/2),...
                          (d1-d5)/2,(d2-d6)/2,...
                            ninclusion,1],...
                        [((d1-d5)/4+d5/2),-((d2-d6)/4+d6/2),...
                          (d1-d5)/2,(d2-d6)/2,...
                            ninclusion,1],...
                        [0,0,...
                         d3,d4,...
                         ninclusion,1] 
                          };     
    case 'two_hollow_circles'
        d1 = settings.d1*1e6; % [um]
        d2 = settings.d2*1e6; % [um]
        d3 = settings.d3*1e6; % [um]
        d4 = settings.d4*1e6; % [um]
        d5 = settings.d5*1e6; % [um]
        d6 = settings.d6*1e6; % [um]
        d7 = settings.d7*1e6; % [um]
        d8 = settings.d8*1e6; % [um]
        assert(d3 < d1,"D3 must be smaller than D1");
        assert(d4 < d1,"D4 must be smaller than D2");
        assert(d1 > 0 && d2 > 0,"D1 and D2 must be positive");
        assert(d5 < d3,"D5 must be smaller than D3");
        assert(d6 < d4,"D6 must be smaller than D4");
        assert(d7 < d5,"D7 must be smaller than D5");
        assert(d8 < d6,"D8 must be smaller than D6");
        assert(d7 > 0 && d8 > 0,"You must have two squares present");
        textures{3} = {nbackground,...
                        [pillar_cent_x,pillar_cent_y,...
                            d1,d2,...
                            ninclusion,settings.staircase_N],...
                        [pillar_cent_x,pillar_cent_y,...
                            d3,d4,...
                            nbackground,settings.staircase_N],...
                        [pillar_cent_x,pillar_cent_y,...
                            d5,d6,...
                            ninclusion,settings.staircase_N],...
                        [pillar_cent_x,pillar_cent_y,...
                            d7,d8,...
                            nbackground,settings.staircase_N]};
    case 'cross_mid_circle'
        d1 = settings.d1*1e6; % [um]
        d2 = settings.d2*1e6; % [um]
        d3 = settings.d3*1e6; % [um]
        d4 = settings.d4*1e6; % [um]
        d5 = settings.d5*1e6; % [um]
        d6 = settings.d6*1e6; % [um]
        assert(d1 > 0 && d2 > 0,"d1 and d2 must be positive");
        assert(d3 < d1/2 && d4 < d2/2,"d3 and d4 must be less than half the cross length");
        assert(d5 <= d1/2 && d6 <= d2/2,"d5 and d6 must be less than half the cross length");
        % check if the central ellipse protrudes out of the pattern
        theta_test = atan((d2/2.-d4)/(d1/2.-d3));
        ellipse_extent = sqrt((d5*cos(theta_test))^2 + (d6*sin(theta_test))^2);
        cross_extent = sqrt((d2/2.-d4)^2 + (d1/2.-d3)^2);
        if (ellipse_extent <= cross_extent || ellipse_extent == 0)
            % then we just have a pure cross
            textures{3} = {nbackground,...
                            [pillar_cent_x,pillar_cent_y,...
                                d1,d2-2*d4,...
                                ninclusion,1],...
                            [pillar_cent_x,pillar_cent_y,...
                                d1-2*d3,d4,...
                                ninclusion,1]};
        else
            % then we have a cross and a circle
            textures{3} = {nbackground,...
                            [pillar_cent_x,pillar_cent_y,...
                                d1,d2-2*d4,...
                                ninclusion,1],...
                            [pillar_cent_x,pillar_cent_y,...
                                d1-2*d3,d2,...
                                ninclusion,1],...
                            [pillar_cent_x,pillar_cent_y,...
                                2*d5,2*d6,...
                                ninclusion,settings.staircase_N]};
        end
            
    case 'hollow_circle'
        d1 = settings.d1*1e6; % [um]
        d2 = settings.d2*1e6; % [um]
        d3 = settings.d3*1e6; % [um]
        d4 = settings.d4*1e6; % [um]
        assert(d3 < d1 || d1 == 0,"D3 must be smaller than D1");
        assert(d4 < d2 || d2 == 0,"D4 must be smaller than D2");
%         assert(d1 > 0 && d2 > 0,"D1 and D2 must be positive");
        if (d3 > 0 && d4 > 0)
            % then we need to add a hollow core
            textures{3} = {nbackground,...
                            [pillar_cent_x,pillar_cent_y,...
                                d1,d2,...
                                ninclusion,settings.staircase_N],...
                            [pillar_cent_x,pillar_cent_y,...
                                d3,d4,...
                                nbackground,settings.staircase_N]};
        else
            % no hollow core needed
            textures{3} = {nbackground,[pillar_cent_x,pillar_cent_y,...
                            d1,d2,...
                            ninclusion,settings.staircase_N]};
        end
    case 'supported_hollow_square'
        d1 = settings.d1*1e6; % [um]
        d2 = settings.d2*1e6; % [um]
        d3 = settings.d3*1e6; % [um]
        d4 = settings.d4*1e6; % [um]
        d9 = settings.d9*1e9; % [um]
        d10 = settings.d10*1e9; % [um]
        assert(d3 < d1,"D3 must be smaller than D1");
        assert(d4 < d1,"D4 must be smaller than D2");
        assert(d1 > 0 && d2 > 0,"D1 and D2 must be positive");
        assert(d9 < d3,"D9 must be smaller than D3")
        assert(d10 < d4,"D10 must be smaller than D4")
        assert(d9 > 0 && d10 > 0,"Inner square must be supported");
        textures{3} = {nbackground,...
                        [pillar_cent_x,pillar_cent_y,...
                            d1,d2,...
                            ninclusion,1],...
                        [pillar_cent_x,pillar_cent_y,...
                            d3,d4,...
                            nbackground,1]};  
    otherwise
        assert(0,'Invalid shape');
end

textures{4} = {n_superstrate_medium};
textures{5} = {n_transmitted_medium};

% set up the profile
% the top layer and bottom layer are uniform media with n_incident_medium
% and sandwich the pillar array layer
transmitted_layer_t = settings.transmitted_layer_t*1e6; % [m] to [um]
incident_layer_t = settings.incident_layer_t*1e6; % [m] to [um]
substrate_layer_t = settings.substrate_layer_t*1e6; % [m] to [um]
superstrate_layer_t = settings.superstrate_layer_t*1e6; % [m] to [um]
pillar_layer_t = settings.pillar_layer_t*1e6; % [m] to [um]
if (substrate_layer_t > 0 && superstrate_layer_t > 0)
    profile = {[transmitted_layer_t,...
                superstrate_layer_t,...
                pillar_layer_t,...
                substrate_layer_t,...
                incident_layer_t],...
                [5,4,3,2,1]};
elseif (substrate_layer_t > 0 && superstrate_layer_t <= 0)
    profile = {[transmitted_layer_t,...
                pillar_layer_t,...
                substrate_layer_t,...
                incident_layer_t],...
                [5,3,2,1]};
elseif (substrate_layer_t <= 0 && superstrate_layer_t > 0)
    profile = {[transmitted_layer_t,...
                superstrate_layer_t,...
                pillar_layer_t,...
                incident_layer_t],...
                [5,4,3,1]};
else
    profile = {[transmitted_layer_t,...
                pillar_layer_t,...
                incident_layer_t],...
                [5,3,1]};
end

% Check for symmetries based on the geometry requested
if (angle_delta == 0 && angle_theta == 0 && pillar_cent_x == 0 && pillar_cent_y == 0)
   sym_x = 1;
   sym_y = 1;
else
   sym_x = 0;
   sym_y = 0;
end

%% Step 1: calculate the eigenmodes
parm = res0;
if sym_x
    parm.sym.x = 0; % put the symmetry axis through the origin 
end
if sym_y
    parm.sym.y = 0; % put the symmetry axis through the origin 
end
% check if the structure should be polarization-independent
if ((strcmp(shape,'hollow_square') && d1 == d2 && d3 == d4) || ...
    (strcmp(shape,'hollow_circle') && d1 == d2 && d3 == d4) || ...
    (strcmp(shape,'cross_mid_circle') && d1 == d2 && d3 == d4 && d5 == d6) || ...
    (strcmp(shape,'two_hollow_squares') && d1==d2 && d3==d4 && d5==d6 && d7==d8) || ...
    (strcmp(shape,'two_hollow_circles') && d1==d2 && d3==d4 && d5==d6 && d7==d8) || ...
    (strcmp(shape,'four_squares') && d1==d2 && d3==d4) || ...
    (strcmp(shape,'square_cog') && d1==d2 && d3==d4 && d5==d6))
    pol_independent = 1;
else
    pol_independent = 0;
end



if pol_independent
    parm.sym.pol = 1; % only compute TE 
end

if debug
    parm.res1.trace = 1;
end
[aa,n_eff] = res1(wavelength,period,textures,nn,k_parallel,angle_delta,parm);

%% Step 2: compute the diffracted waves
result = res2(aa,profile);

% Get the information for propagated transmitted waves from the bottom
if pol_independent
    res_TE = result.TEinc_bottom_transmitted;
    % Extract the phase information
    % Get the complex zeroth order amplitude
    N_orders = size(res_TE.order,1);
    if N_orders > 1
        order_ind = ceil(N_orders/2);
    else 
        order_ind = 1;
    end
    assert(res_TE.order(order_ind,1) == 0);
    assert(res_TE.order(order_ind,2) == 0);

    amplitude_TM = 0;
    amplitude_TE = res_TE.amplitude_TE(order_ind);
else
    res_TM = result.TMinc_bottom_transmitted;
    res_TE = result.TEinc_bottom_transmitted;
    % Extract the phase information
    % Get the complex zeroth order amplitude
    N_orders = size(res_TM.order,1);
    if N_orders > 1
        order_ind = ceil(N_orders/2);
    else 
        order_ind = 1;
    end
    assert(res_TM.order(order_ind,1) == 0);
    assert(res_TM.order(order_ind,2) == 0);
    assert(res_TE.order(order_ind,1) == 0);
    assert(res_TE.order(order_ind,2) == 0);

    amplitude_TM = res_TM.amplitude_TM(order_ind);
    amplitude_TE = res_TE.amplitude_TE(order_ind);
end


out.amplitude_TM = amplitude_TM;
out.amplitude_TE = amplitude_TE;
out.n_eff = n_eff{3};
out.result = result; % store the full result too
retio;
end

