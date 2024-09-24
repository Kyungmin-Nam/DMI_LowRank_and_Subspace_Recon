% demo_step4_bart_generate_AMARES_fiting_noise_v1.m
% Written by Kyungmin Nam
% Email: K.M.Nam@umcutrecht.nl
% Started: 07/19/2024, Last modified: 09/24/2024

%% Clean slate
close all; clear variables; clc;

%% Set source paths
package_path0 = 'E:\mfiles_kmin\MRM_dSPICE\demo\src';
package_path1 = 'E:\mfiles_kmin\MRM_dSPICE\demo\thirdparty';

%% Add source paths to search path
addpath(genpath(package_path0));
addpath(genpath(package_path1));

%% Define input files
input_path  = 'E:\mfiles_kmin\MRM_dSPICE\demo\data\output\noise\';
output_path = [input_path,'AMARES\output\']; mkdir(output_path);

%% Define a path to BART
bart_path = '/home/kmin/bart';

%% Start a stopwatch timer
start_time = tic;

%% Set up BART commands
%--------------------------------------------------------------------------
% Define a BART command
%--------------------------------------------------------------------------
if ispc
    command_prefix = 'wsl';
else
    command_prefix = '';
end
bart_command = sprintf('%s %s/bart', command_prefix, bart_path);

%--------------------------------------------------------------------------
% Translate from a Windows path to a WSL path 
%--------------------------------------------------------------------------
if ispc
    bart_output_path = strrep(output_path, '\', '/');
    bart_output_path = sprintf('/mnt/%s/%s/', lower(bart_output_path(1)), bart_output_path(4:end));
else
    bart_output_path = sprintf('%s/', output_path);
end

%% Define simulation parameters
SimParam.BW           = 5000;      % spectral bandwidth [Hz]
SimParam.NP           = 1024;      % number of spectral points
SimParam.B0           = 7;         % main field strength [T]
SimParam.GyromagRatio = 6.53569e6; % gyromagnetic ratio [Hz/T]

%% Calculate a time axis [sec]
SimParam.time = (0:SimParam.NP-1) / SimParam.BW; % [sec]
SimParam.Freq = SimParam.B0 * SimParam.GyromagRatio;
SimParam.ppmwindow = 109.2900;

%% Define metabolites
carrier_frequency_ppm = 4.7; % [ppm]

%% Calculate a ppm axis [ppm]
df = SimParam.BW / SimParam.NP; % [Hz]
freq_axis = (-floor(SimParam.NP/2):ceil(SimParam.NP/2)-1).' * df; % [Hz]
SimParam.ppm_axis = freq_axis / (SimParam.GyromagRatio * SimParam.B0 * 1e-6) +...
    carrier_frequency_ppm; % [Hz] / [Hz/ppm] => [ppm]

%--------------------------------------------------------------------------
% ksp_lowres (N1 x N2 x N3 x 1 x 1 x M)
%--------------------------------------------------------------------------
ksp_filename = 'ksp_chcomb_nslvl0_nsitr0';

%% Read a .cfl file
%--------------------------------------------------------------------------
% ksp (N1 x N2 x N3 x Nc x 1 x M)
%--------------------------------------------------------------------------
cfl_file = fullfile(input_path, ksp_filename);
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
ksp = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
[N1,N2,N3,~,~,M] = size(ksp);

%% Fit using AMARES
SimPhantom = calculate_FFT_AMARES_Fitting_v1(SimParam,ksp);

%% Write a .mat file
%--------------------------------------------------------------------------
% AMARES_FFT_nslvl%d_nsitr%d
%--------------------------------------------------------------------------
ksp_filename = 'AMARES_FFT_nslvl0_nsitr0';
mat_file = fullfile(output_path, ksp_filename);
tstart = tic; fprintf('%s: Writing a .mat file: %s... ', datetime, cfl_file);
save(mat_file, 'SimPhantom');
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Set noise parameters
num_noise = 4;  % 4 noise level
itr_noise = 30; % 30 iteration

for idx_nslvl = 3%:num_noise
    for idx_itr = 1%:itr_noise
        
    ksp_filename = sprintf('ksp_chcomb_nslvl%d_nsitr%d', idx_nslvl, idx_itr);
    
    %% Read a .cfl file
    %--------------------------------------------------------------------------
    % ksp (N1 x N2 x N3 x Nc x 1 x M)
    %--------------------------------------------------------------------------
    cfl_file = fullfile(input_path, ksp_filename);
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    ksp = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    [N1,N2,N3,Nc,~,M] = size(ksp);
    
    %% Fit using AMARES
    SimPhantom = calculate_FFT_AMARES_Fitting_v1(SimParam,ksp);

    %% Write a .mat file
    %--------------------------------------------------------------------------
    % AMARES_FFT_nslvl%d_nsitr%d
    %--------------------------------------------------------------------------
    ksp_filename = sprintf('AMARES_FFT_nslvl%d_nsitr%d', idx_nslvl, idx_itr);
    mat_file = fullfile(output_path, ksp_filename);
    tstart = tic; fprintf('%s: Writing a .mat file: %s... ', datetime, cfl_file);
    save(mat_file, 'SimPhantom');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    end
end
%%
AMARES_amp = (SimPhantom.AMARESmaps.Amplitude);

%% Display
figure('color','w'); 
imagesc(AMARES_amp(:,:,7,4)); axis image off; colormap(gray); 
colorbar; title('HDO map');

figure('color','w'); 
imagesc(AMARES_amp(:,:,7,3)); axis image off; colormap(gray); 
colorbar; title('Glx map');

figure('color','w'); 
imagesc(AMARES_amp(:,:,7,2)); axis image off; colormap(gray); 
colorbar; title('Glc map');

figure('color','w'); 
imagesc(AMARES_amp(:,:,7,1)); axis image off; colormap(gray); 
colorbar; title('Lipid map');



































