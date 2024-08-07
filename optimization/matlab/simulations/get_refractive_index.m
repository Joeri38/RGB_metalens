function [n,k] = get_refractive_index(material,wavelengths)
%GET_REFRACTIVE_INDEX returns a list of refractive indices for requested
%wavelengths. interpolation is linear. extrapolation will return NaN. 

% params:
%   -material [string]: either 'SiO2','Si','Au','YIG','TiO2_ALD','TiO2_LGA'
%                       or 'Ag','PMMA','air','vacuum', 'Pt'
%   -wavelengths [um]: vector or matrix of requested wavelengths. 

% returns:
%   -n: refractive index [1], same shape as "wavelengths"
%   -k: absorption [1], same shape as "wavelengths"

% Sources:
%   'Au': https://refractiveindex.info/?shelf=main&book=Au&page=McPeak
%         https://doi.org/10.1021/ph5004237
%   'Si': https://refractiveindex.info/?shelf=main&book=Si&page=Green-2008
%         https://doi.org/10.1016/j.solmat.2008.06.009
%         https://refractiveindex.info/?shelf=main&book=Si&page=Li-293K
%         https://doi.org/10.1063/1.555624
%   'SiO2': https://refractiveindex.info/?shelf=main&book=SiO2&page=Malitson
%           https://doi.org/10.1364/JOSA.55.001205
%   'YIG': https://doi.org/10.3390/nano8050355 (Table 2, n2 column)
%   'TiO2_ALD': Capasso group internal experimental data. ALD-1 TiO2. 
%   'TiO2_LGA': Maryna's measurements of LGA thin films deposited TiO2
%   using e-beam evaporation. 
%   'Ag': https://refractiveindex.info/?shelf=main&book=Ag&page=Johnson
%         https://doi.org/10.1103/PhysRevB.6.4370
%   'PMMA': https://refractiveindex.info/?shelf=organic&book=poly%28methyl_methacrylate%29&page=Zhang-Mitsubishi
%           https://doi.org/10.1364/AO.383831
%           https://doi.org/10.1016/j.jqsrt.2020.107063
%   'Pt': https://doi.org/10.1063/1.3243762
%         https://refractiveindex.info/tmp/data/main/Pt/Werner.txt
%   'air': n=1,k=0
%   'vacuum': n=1,k=0

interp_method = 'linear';
extrap_val = NaN;

load('./refractive_index_info.mat');
% load('./refractive_index_info.mat');

material = upper(material);

switch material
    case 'AU'
        n_table = Au_n;
        k_table = Au_k;
    case 'AG'
        n_table = Ag_n;
        k_table = Ag_k;
    case 'SI'
        n_table = Si_n;
        k_table = Si_k;
    case 'SIO2'
        n_table = SiO2_n;
        k_table = SiO2_k;
    case 'YIG'
        n_table = YIG_n;
        k_table = YIG_k;
    case 'TIO2_ALD'
        n_table = TiO2_ALD_n;
        k_table = TiO2_ALD_k;
    case 'TIO2_LGA'
        n_table = TiO2_LGA_n;
        k_table = TiO2_LGA_k;
    case 'PMMA'
        n_table = PMMA_n;
        k_table = PMMA_k;
    case 'AIR'
        wl = [min(wavelengths)-1;max(wavelengths)+1];
        n = [1;1];
        k = [0;0];
        n_table = table(wl,n);
        k_table = table(wl,k);
    case 'VACUUM'
        wl = [min(wavelengths)-1;max(wavelengths)+1];
        n = [1;1];
        k = [0;0];
        n_table = table(wl,n);
        k_table = table(wl,k);
    case 'PT'
        n_table = Pt_n;
        k_table = Pt_k;
    otherwise
        assert(false,sprintf('%s: material not available',mfilename));
end

% perform interpolation for n
x = n_table.wl; % [um]
v = n_table.n; % [1]
n = interp1(x,v,wavelengths,interp_method,extrap_val); % [1]

% perform interpolation for k
x = k_table.wl; % [um]
v = k_table.k; % [1]
k = interp1(x,v,wavelengths,interp_method,extrap_val); % [1]

end

