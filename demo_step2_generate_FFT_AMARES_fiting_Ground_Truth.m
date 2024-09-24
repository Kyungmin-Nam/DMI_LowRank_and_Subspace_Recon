%--------------------------------------------------------------------------
% demo_step2_generate_FFT_AMARES_fiting_Ground_Truth.m
% Written by Kyungmin Nam
% Email: K.M.Nam@umcutrecht.nl
% Started: 09/16/2024, Last modified: 
%--------------------------------------------------------------------------

%% Clean slate
close all; clear variables; clc;

%% Set source paths
package_path = 'E:\mfiles_kmin\MRM_dSPICE\demo\thirdparty';

%% Add source paths to search path
addpath(genpath(package_path));

%% Define input files
main_path   = 'E:\mfiles_kmin\MRM_dSPICE\demo'; 
input_path  = [main_path, '\data\output\'];
output_path = [input_path,'AMARES\output\'];

%% Start a stopwatch timer
start_time = tic;

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
img_filename = 'img_12x18x15';

%% Read a .cfl file
%--------------------------------------------------------------------------
% ksp (N1 x N2 x N3 x Nc x 1 x M)
%--------------------------------------------------------------------------
cfl_file = fullfile(input_path, img_filename);
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
img_t = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
[N1,N2,N3,~,~,M] = size(img_t);

%% Fit using AMARES
SimPhantom = generate_FFT_AMARES_Fitting(main_path,SimParam,img_t);
AMARES_amp = SimPhantom.AMARESmaps.Amplitude;

%% Display metabolite maps
figure('Color','w'); imagesc(AMARES_amp(:,:,8,1)); 
axis image off; colormap(gray); colorbar;
title('Ground truth, Lipid map');

figure('Color','w'); imagesc(AMARES_amp(:,:,8,2)); 
axis image off; colormap(gray); colorbar;
title('Ground truth, Glx map');

figure('Color','w'); imagesc(AMARES_amp(:,:,8,3)); 
axis image off; colormap(gray); colorbar;
title('Ground truth, Glc map');

figure('Color','w'); imagesc(AMARES_amp(:,:,8,4)); 
axis image off; colormap(gray); colorbar;
title('Ground truth, HDO map');

%% Write a .mat file
%--------------------------------------------------------------------------
% AMARES_FFT_nslvl%d_nsitr%d
%--------------------------------------------------------------------------
img_filename = 'AMARES_FFT_nslvl0_nsitr0';
mat_file = fullfile(output_path, img_filename);
tstart = tic; fprintf('%s: Writing a .mat file: %s... ', datetime, cfl_file);
save(mat_file, 'SimPhantom');
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

