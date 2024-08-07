function [all_meta_atoms,completed_jobs] = get_circular_meta_atom_GD(input)
%GET_META_ATOM_GD performs RETICOLO simulations for the circular pillar
% or holey meta-atoms as a full parameter sweep in all dimensions, then
% returns the group delay information as a Table.

% Date created: 20211209
% Author: Soon Wei Daniel Lim

% params:
%   input: struct with the following fields:
%       angle_delta [deg]
%       angle_theta [deg] 
%       substrate_material [string]: see allowed materials in get_refractive_index
%       superstrate_material [string]: see allowed materials in get_refractive_index
%       num_fourier [int]: total number of fourier modes is
%           (2*num_fourier+1)^2
%       material_inside_circle [string]: see allowed materials in get_refractive_index
%       material_outside_circle [string]: see allowed materials in get_refractive_index
%       cent_x [m]: pillar position in unit cell
%       cent_y [m]: pillar position in unit cell
%       patterned_layer_t_range [m, vector]: vector of pillar layer
%           thicknesses
%       period_range [m, vector]: vector of periodicities
%       staircase_N [int]: staircase approximation for arcs
%       diameter_range [m, vector]: vector of circle diameters. you are
%           responsible to make sure that this is compatible with the
%           periodicity.
%       lambda_range [m,vector]: vector of wavelengths
%       lambda0 [m]: central wavelength, must be in lambda_range
%       i [bool]: 1 or 0
%       multiplex [bool]: 1 or 0. if 1, then each combination of patterned
%           layer thickness, period, and diameter is sampled. otherwise,
%           only the specified thickness, period, and diameter is computed.

% returns:
%   all_meta_atoms: table of all meta-atom information, including fitted GD
%       and GDD values. 
%   completed_jobs: table of all raw information from the run. 

fs = 1e-15;
c = 299792458; % [m/s]
lambda0 = input.lambda0; % [m] central wavelength
omega0 = c*2*pi/lambda0; % [rad/s] central frequency
assert(lambda0 <= max(input.lambda_range));
assert(lambda0 >= min(input.lambda_range));

angle_delta = input.angle_delta; % [deg]
angle_theta = input.angle_theta; % [deg]
lambda_range = input.lambda_range; % [m]
N_lambda = length(lambda_range);
[n_incident,k_incident] = get_refractive_index(input.substrate_material,...
                                        lambda_range*1e6); % [1]
% n_incident_medium_range = n_incident + 1i*k_incident;
n_incident_medium_range = n_incident;
[n_substrate,k_substrate] = get_refractive_index(input.substrate_material,...
                                        lambda_range*1e6); % [1]
% n_substrate_medium_range = n_substrate + 1i*k_substrate;
n_substrate_medium_range = n_substrate;
[n_superstrate,k_superstrate] = get_refractive_index(input.superstrate_material,...
                                        lambda_range*1e6); % [1]
% n_superstrate_medium_range = n_superstrate + 1i*k_superstrate;
n_superstrate_medium_range = n_superstrate;
[n_transmitted,k_transmitted] = get_refractive_index(input.superstrate_material,...
                                        lambda_range*1e6); % [1]
% n_transmitted_medium_range = n_transmitted + 1i*k_transmitted;
n_transmitted_medium_range = n_transmitted;

num_fourier = input.num_fourier;
[n_inclusion,k_inclusion] = get_refractive_index(input.material_inside_circle,...
                                        lambda_range*1e6); % index of the circular region
% ninclusion_range = n_inclusion + 1i*k_inclusion;                              
ninclusion_range = n_inclusion;                              
[n_background,k_background] = get_refractive_index(input.material_outside_circle,...
                                        lambda_range*1e6); % index outside circular region
k_background = 0; % set absorption to 0
% nbackground_range = n_background + 1i*k_background;
nbackground_range = n_background;
pillar_cent_x = input.cent_x; % [m]
pillar_cent_y = input.cent_y; % [m]
pillar_layer_t_range = input.patterned_layer_t_range; % [m]
N_thickness = length(pillar_layer_t_range);
period_range = input.period_range; % [m]
N_period = length(period_range);
staircase_N = input.staircase_N; 
incident_layer_t = 0e-6; % [m]
substrate_layer_t = 2e-6; % [m]
superstrate_layer_t = 0e-6; % [m]
transmitted_layer_t = 0e-6; % [m]
diameter1_range = input.diameter1_range;
diameter2_range = input.diameter2_range;
%diameter3_range = input.diameter3_range;
%diameter4_range = input.diameter4_range;
N_diameter1 = length(diameter1_range);
N_diameter2 = length(diameter2_range);
%N_diameter3 = length(diameter3_range);
%N_diameter4 = length(diameter4_range);
debug = input.debug;

%% Set up the containers

d1_list = [];
d2_list = [];
d3_list = [];
d4_list = []; 
d5_list = [];
d6_list = [];
d7_list = [];
d8_list = [];
shape_list = {};
lambda_list = [];
n_incident_medium_list = [];
n_substrate_medium_list = [];
n_superstrate_medium_list = [];
n_transmitted_medium_list = [];
ninclusion_list = [];
nbackground_list = [];
pillar_layer_t_list = [];
period_list = [];

%% Generate the queue

switch input.multiplex
    case 1
        if debug
            N_jobs = N_lambda*N_thickness*N_period*N_diameter;
            fprintf('Number of jobs (multiplexed): %d\n',N_jobs);
        end
        % here, d1 is the diameter. d3=d4=0 because this is a solid circle

        for i_thickness = 1:N_thickness
            this_thickness = pillar_layer_t_range(i_thickness);
            for i_period = 1:N_period
                this_period = period_range(i_period);
                for i_diameter = 1:N_diameter1
                    this_diameter1 = diameter1_range(i_diameter);
                    for j_diameter = 1:N_diameter2
                        this_diameter2 = diameter2_range(j_diameter);
                        %if j_diameter >= i_diameter
                            for i_lambda = 1:N_lambda
                                    this_lambda = lambda_range(i_lambda);
                                    this_n_incident_medium = n_incident_medium_range(i_lambda);
                                    this_n_substrate_medium = n_substrate_medium_range(i_lambda);
                                    this_n_superstrate_medium = n_superstrate_medium_range(i_lambda);
                                    this_n_transmitted_medium = n_transmitted_medium_range(i_lambda);
                                    this_ninclusion = ninclusion_range(i_lambda);
                                    this_nbackground = nbackground_range(i_lambda);
                                    d1 = this_diameter1;
                                    d2 = d1;
                                    d3 = this_diameter2;
                                    d4 = d3;
                                    d5 = this_diameter2;
                                    d6 = d5;
                                    d7 = this_diameter1;
                                    d8 = d7;
                                    shape_list = vertcat(shape_list,{'four_circles'});
                                    d1_list = vertcat(d1_list,d1);
                                    d2_list = vertcat(d2_list,d2);
                                    d3_list = vertcat(d3_list,d3);
                                    d4_list = vertcat(d4_list,d4);
                                    d5_list = vertcat(d5_list,d5);
                                    d6_list = vertcat(d6_list,d6);
                                    d7_list = vertcat(d7_list,d7);
                                    d8_list = vertcat(d8_list,d8);
                                    lambda_list = vertcat(lambda_list,this_lambda);
                                    n_incident_medium_list = vertcat(n_incident_medium_list,this_n_incident_medium);
                                    n_substrate_medium_list = vertcat(n_substrate_medium_list,this_n_substrate_medium);
                                    n_superstrate_medium_list = vertcat(n_superstrate_medium_list,this_n_superstrate_medium);
                                    n_transmitted_medium_list = vertcat(n_transmitted_medium_list,this_n_transmitted_medium);
                                    ninclusion_list = vertcat(ninclusion_list,this_ninclusion);
                                    nbackground_list = vertcat(nbackground_list,this_nbackground);
                                    pillar_layer_t_list = vertcat(pillar_layer_t_list,this_thickness);
                                    period_list = vertcat(period_list,this_period);
                            end
                        %end                          
                    end
                end
            end
        end
    case 0

        % check that the variable vectors have the same size
        assert(N_period==N_thickness,"N_period not equal to N_thickness");
        assert(N_period==N_diameter,"N_period not equal to N_diameter");
        
        for i_period = 1:N_period
            this_thickness = pillar_layer_t_range(i_period);
            this_period = period_range(i_period);
            this_diameter = diameter_range(i_period);
            for i_lambda = 1:N_lambda
                this_lambda = lambda_range(i_lambda);
                this_n_incident_medium = n_incident_medium_range(i_lambda);
                this_n_substrate_medium = n_substrate_medium_range(i_lambda);
                this_n_superstrate_medium = n_superstrate_medium_range(i_lambda);
                this_n_transmitted_medium = n_transmitted_medium_range(i_lambda);
                this_ninclusion = ninclusion_range(i_lambda);
                this_nbackground = nbackground_range(i_lambda);
                d1 = this_diameter;
                d2 = d1;
                d3 = 0;
                d4 = d3;
                d5 = 0;
                d6 = 0;
                d7 = 0;
                d8 = 0;
                shape_list = vertcat(shape_list,{'hollow_circle'});
                d1_list = vertcat(d1_list,d1);
                d2_list = vertcat(d2_list,d2);
                d3_list = vertcat(d3_list,d3);
                d4_list = vertcat(d4_list,d4);
                d5_list = vertcat(d5_list,d5);
                d6_list = vertcat(d6_list,d6);
                d7_list = vertcat(d7_list,d7);
                d8_list = vertcat(d8_list,d8);
                lambda_list = vertcat(lambda_list,this_lambda);
                n_incident_medium_list = vertcat(n_incident_medium_list,this_n_incident_medium);
                n_substrate_medium_list = vertcat(n_substrate_medium_list,this_n_substrate_medium);
                n_superstrate_medium_list = vertcat(n_superstrate_medium_list,this_n_superstrate_medium);
                n_transmitted_medium_list = vertcat(n_transmitted_medium_list,this_n_transmitted_medium);
                ninclusion_list = vertcat(ninclusion_list,this_ninclusion);
                nbackground_list = vertcat(nbackground_list,this_nbackground);
                pillar_layer_t_list = vertcat(pillar_layer_t_list,this_thickness);
                period_list = vertcat(period_list,this_period);
            end
        end
        if debug
            fprintf('Number of jobs (not multiplexed): %d\n',N_lambda*N_period);
        end
    otherwise
        assert(false,'input.multiplex is either 1 or 0');
end

N_jobs = size(lambda_list,1);

d1 = d1_list;
d2 = d2_list;
d3 = d3_list;
d4 = d4_list;
d5 = d5_list; 
d6 = d6_list;
d7 = d7_list;
d8 = d8_list;
shape = shape_list;
lambda = lambda_list;
nbackground = nbackground_list;

angle_delta = angle_delta*ones([N_jobs,1]);
angle_theta = angle_theta*ones([N_jobs,1]);
n_incident_medium = n_incident_medium_list;
n_substrate_medium = n_substrate_medium_list;
n_superstrate_medium = n_superstrate_medium_list;
n_transmitted_medium = n_transmitted_medium_list;
num_fourier = num_fourier*ones([N_jobs,1]);
ninclusion = ninclusion_list;
pillar_cent_x = pillar_cent_x*ones([N_jobs,1]);
pillar_cent_y = pillar_cent_y*ones([N_jobs,1]);
pillar_layer_t = pillar_layer_t_list;
period = period_list;
staircase_N = staircase_N*ones([N_jobs,1]);
incident_layer_t = incident_layer_t*ones([N_jobs,1]);
substrate_layer_t = substrate_layer_t*ones([N_jobs,1]);
superstrate_layer_t = superstrate_layer_t*ones([N_jobs,1]);
transmitted_layer_t = transmitted_layer_t*ones([N_jobs,1]);

all_jobs = table(shape,d1,d2,d3,d4,d5,d6,d7,d8,...
                lambda,nbackground,...
                angle_delta,angle_theta,n_incident_medium,...
                n_substrate_medium,n_superstrate_medium,n_transmitted_medium,...
                num_fourier,ninclusion,pillar_cent_x,pillar_cent_y,...
                pillar_layer_t,period,staircase_N,incident_layer_t,...
                substrate_layer_t,superstrate_layer_t,...
                transmitted_layer_t);

%% Run the main sequence
if isempty(gcp('nocreate'))
    parpool;
end
parforloop_tic = tic;
completed_jobs = run_table_parfor(all_jobs);
parfor_loop_total_time = toc(parforloop_tic);
if debug
    fprintf('Parfor loop total time elapsed: %g s, time per job: %g s\n',parfor_loop_total_time,parfor_loop_total_time/N_jobs);
end
% runfile = mfilename;

%% Now we need to get the chromatic properties (GD,GDD)
N_jobs = size(completed_jobs,1);
all_meta_atoms = table();

for row_idx = 1:N_jobs
    row = completed_jobs(row_idx,:);
    n_eff = row.n_eff{1};
    real_n_eff = real(n_eff(imag(n_eff)<1e-9));
    max_real_n_eff = max(real_n_eff);
    completed_jobs.max_real_n_eff(row_idx) = max_real_n_eff;
end

current_row = 1; % current row being examined

while current_row <= N_jobs
    ref_row = completed_jobs(current_row,:);
    all_meta_atoms = vertcat(all_meta_atoms,ref_row); % this row is assumed to be a new unique meta-atom
    temp_max_real_n_eff = ref_row.max_real_n_eff;
    temp_lambda = ref_row.lambda;
    temp_amplitude_TE = ref_row.amplitude_TE;
    % also compute the expected phase based on effective medium
    % theory (just use the largest real index)
    temp_n_substrate_medium = ref_row.n_substrate_medium;
    temp_n_superstrate_medium = ref_row.n_superstrate_medium;
    temp_pillar_layer_t = ref_row.pillar_layer_t;
    temp_amplitude_TE_from_neff = estimate_complex_amplitude(temp_n_substrate_medium(:),...
                                                            temp_max_real_n_eff(:),...
                                                            temp_n_superstrate_medium(:),...
                                                            temp_pillar_layer_t(:),...
                                                            temp_lambda(:));
    % now check the following rows for duplicate geometries which indicate
    % a different wavelength
    found_different_row = 0;
    while ~found_different_row
        current_row = current_row + 1;
        if current_row > N_jobs
            break;
        end
        test_row = completed_jobs(current_row,:);
        if (strcmp(ref_row.shape{1},test_row.shape{1}) && ...
            ref_row.d1 == test_row.d1 && ...
            ref_row.d2 == test_row.d2 && ...
            ref_row.d3 == test_row.d3 && ...
            ref_row.d4 == test_row.d4 && ...
            ref_row.d5 == test_row.d5 && ...
            ref_row.d6 == test_row.d6 && ...
            ref_row.d7 == test_row.d7 && ...
            ref_row.d8 == test_row.d8 && ...
            ref_row.period == test_row.period && ...
            ref_row.pillar_layer_t == test_row.pillar_layer_t)
            % then we have found a similar geometry but at a different
            % wavelength
            temp_max_real_n_eff = vertcat(temp_max_real_n_eff,test_row.max_real_n_eff);
            temp_lambda = vertcat(temp_lambda,test_row.lambda);
            temp_amplitude_TE = vertcat(temp_amplitude_TE,test_row.amplitude_TE);
            % also compute the expected phase based on effective medium
            % theory (just use the largest real index)
            temp_n_substrate_medium = test_row.n_substrate_medium;
            temp_n_superstrate_medium = test_row.n_superstrate_medium;
            temp_pillar_layer_t = test_row.pillar_layer_t;
            estimated_complex_amplitude = estimate_complex_amplitude(temp_n_substrate_medium(:),...
                                                            test_row.max_real_n_eff(:),...
                                                            temp_n_superstrate_medium(:),...
                                                            temp_pillar_layer_t(:),...
                                                            test_row.lambda(:));
            temp_amplitude_TE_from_neff = vertcat(temp_amplitude_TE_from_neff,...
                                                  estimated_complex_amplitude(:));
        else
            % then this row is different
            found_different_row = 1;
        end
    end
    % now extract the statistics 
    [temp_lambda,I] = sort(temp_lambda);
    temp_max_real_n_eff = temp_max_real_n_eff(I);
    temp_amplitude_TE = temp_amplitude_TE(I);
    temp_omega = fs*c*2*pi./temp_lambda; % [rad/fs]
    temp_phase = unwrap(angle(temp_amplitude_TE)); % [rad]
    if length(temp_phase) > 3
%         [p,S] = polyfit(temp_omega,temp_phase,2); % quadratic fit
        % also try out getting confidence intervals
        [f,gof] = fit(temp_omega,temp_phase,'poly2');
        p = [f.p1,f.p2,f.p3];
        confidence_intervals = confint(f);
        p1_error = 0.5*(confidence_intervals(2,1) - confidence_intervals(1,1));
        p2_error = 0.5*(confidence_intervals(2,2) - confidence_intervals(1,2));
        p3_error = 0.5*(confidence_intervals(2,3) - confidence_intervals(1,3));
        GDD = 2*p(1); % [fs^2]
        GDD_error = 2*p1_error;
        omega0_fs = omega0*fs; % [rad/fs]
        GD = p(2)+ GDD*omega0_fs; % [fs]
        GD_error = sqrt(p2_error^2 + (GDD_error*omega0_fs)^2);
        phase = p(3) + GD*omega0_fs - 0.5*GDD*omega0_fs^2; % [rad]
        phase_error = sqrt(p3_error^3 + (GD_error*omega0_fs)^2 + (0.5*GDD_error*omega0_fs^2)^2);
        % check the R^2 value
%         r2 = 1 - (S.normr/norm(temp_phase - mean(temp_phase)))^2;
        r2 = gof.rsquare;
        % also do the fitting for the complex amplitude from the effective
        % index
        temp_phase_from_n_eff = unwrap(angle(temp_amplitude_TE_from_neff));
        [f2,gof2] = fit(temp_omega(:),temp_phase_from_n_eff(:),'poly2');
        p2 = [f2.p1,f2.p2,f2.p3];
        GDD2 = 2*p2(1); % [fs^2]
        GD2 = p2(2)+ GDD2*omega0_fs; % [fs]
        phase2 = p2(3) + GD2*omega0_fs - 0.5*GDD2*omega0_fs^2; % [rad]
        r2_from_neff = gof2.rsquare;
    else
        phase = NaN;
        GDD = NaN;
        GD = NaN;
        r2 = NaN;
        phase_error = NaN;
        GD_error = NaN;
        GDD_error = NaN;
        phase2 = NaN;
        GDD2 = NaN;
        GD2 = NaN;
        r2_from_neff = NaN;
    end
    
    atom_idx = size(all_meta_atoms,1);
    mean_n_eff(atom_idx,1) = mean(temp_max_real_n_eff);
    min_n_eff(atom_idx,1) = min(temp_max_real_n_eff);
    max_n_eff(atom_idx,1) = max(temp_max_real_n_eff);
    delta_n_eff(atom_idx,1) = max(temp_max_real_n_eff) - min(temp_max_real_n_eff);
    red_n_eff(atom_idx,1) = temp_max_real_n_eff(end); % last wavelength is the longest
    blue_n_eff(atom_idx,1) = temp_max_real_n_eff(1);
    phase_list(atom_idx,1) = phase;
    GD_list(atom_idx,1) = GD;
    GDD_list(atom_idx,1) = GDD;
    r2_list(atom_idx,1) = r2;
    phase_error_list(atom_idx,1) = phase_error;
    GD_error_list(atom_idx,1) = GD_error;
    GDD_error_list(atom_idx,1) = GDD_error;
    phase_from_neff_list(atom_idx,1) = phase2;
    GD_from_neff_list(atom_idx,1) = GD2;
    GDD_from_neff_list(atom_idx,1) = GDD2;
    r2_from_neff_list(atom_idx,1) = r2_from_neff;
%     
%     % stop the loop (temporary)
%     current_row = N_jobs+1;
end

all_meta_atoms.mean_n_eff = mean_n_eff;
all_meta_atoms.min_n_eff = min_n_eff;
all_meta_atoms.max_n_eff = max_n_eff;
all_meta_atoms.delta_n_eff = delta_n_eff;
all_meta_atoms.red_n_eff = red_n_eff;
all_meta_atoms.blue_n_eff = blue_n_eff;
all_meta_atoms.phase = phase_list;
all_meta_atoms.GD = GD_list;
all_meta_atoms.GDD = GDD_list;
all_meta_atoms.r2 = r2_list;
all_meta_atoms.phase_error = phase_error_list;
all_meta_atoms.GD_error = GD_error_list;
all_meta_atoms.GDD_error = GDD_error_list;
all_meta_atoms.phase_from_neff = phase_from_neff_list;
all_meta_atoms.GD_from_neff = GD_from_neff_list;
all_meta_atoms.GDD_from_neff = GDD_from_neff_list;
all_meta_atoms.r2_from_neff = r2_from_neff_list;
all_meta_atoms = removevars(all_meta_atoms,{'lambda','nbackground','amplitude_TE','amplitude_TM','n_eff','max_real_n_eff'});


end

function completed_jobs = run_table(all_jobs)
    completed_jobs = all_jobs;
    N_jobs = size(all_jobs,1);
    for i = 1:N_jobs
%         fprintf('Running job %d of %d...',i,N_jobs);
        row = all_jobs(i,:);
        out = run_table_row(row);
        completed_jobs.amplitude_TM(i) = out.amplitude_TM;
        completed_jobs.amplitude_TE(i) = out.amplitude_TE;
        completed_jobs.n_eff{i} = out.n_eff;
    end
end

function completed_jobs = run_table_parfor(all_jobs)
    
    completed_jobs = all_jobs;
    N_jobs = size(all_jobs,1);
    output_all = cell(N_jobs,1);
    parfor i = 1:N_jobs
    %for i = 1:N_jobs
        row = all_jobs(i,:);
        tic;
        out = run_table_row(row);
        duration = toc;
        out_struct = struct;
        out_struct.amplitude_TM = out.amplitude_TM;
        out_struct.amplitude_TE = out.amplitude_TE;
        out_struct.n_eff = out.n_eff;
        output_all{i} = out_struct;
        fprintf('Job number %d of %d done! Time = %g s\n',i,N_jobs,duration);
    end

    % process the results
    for i = 1:N_jobs
        completed_jobs.amplitude_TM(i) = output_all{i}.amplitude_TM;
        completed_jobs.amplitude_TE(i) = output_all{i}.amplitude_TE;
        completed_jobs.n_eff{i} = output_all{i}.n_eff;
    end
end

function out = run_table_row(row)
    settings.shape = row.shape{1};
    settings.d1 = row.d1;
    settings.d2 = row.d2;
    settings.d3 = row.d3;
    settings.d4 = row.d4;
    settings.d5 = row.d5;
    settings.d6 = row.d6;
    settings.d7 = row.d7;
    settings.d8 = row.d8;
    settings.wavelength = row.lambda;
    settings.nbackground = row.nbackground;
    settings.angle_delta = row.angle_delta;
    settings.angle_theta = row.angle_theta;
    settings.n_incident_medium = row.n_incident_medium;
    settings.n_substrate_medium = row.n_substrate_medium;
    settings.n_superstrate_medium = row.n_superstrate_medium;
    settings.n_transmitted_medium = row.n_transmitted_medium;
    settings.num_fourier = row.num_fourier;
    settings.ninclusion = row.ninclusion;
    settings.pillar_cent_x = row.pillar_cent_x;
    settings.pillar_cent_y = row.pillar_cent_y;
    settings.pillar_layer_t = row.pillar_layer_t;
    settings.period = row.period;
    settings.staircase_N = row.staircase_N;
    settings.incident_layer_t = row.incident_layer_t;
    settings.substrate_layer_t = row.substrate_layer_t;
    settings.superstrate_layer_t = row.superstrate_layer_t;
    settings.transmitted_layer_t = row.transmitted_layer_t;
    settings.debug = false;
    out = achrom_transmission(settings);
end

function t = estimate_complex_amplitude(n_substrate,n_eff,n_superstrate,t,lambda0)
% Use the analytic Fabry-Perot formula to estimate the complex transmission
% amplitude
k0 = 2*pi/lambda0; % [1/m]
k0L = k0*t; % [1]
c_t12 = 2*n_substrate./(n_substrate + n_eff);
c_t23 = 2*n_eff./(n_eff + n_superstrate);
c_r21 = (n_eff - n_substrate)./(n_eff + n_substrate);
c_r23 = (n_eff - n_superstrate)./(n_eff + n_superstrate);
t = sqrt(n_superstrate./n_substrate).*c_t12.*c_t23.*exp(1i*n_eff.*k0L)./(1-c_r21.*c_r23.*exp(2i*n_eff.*k0L));
end