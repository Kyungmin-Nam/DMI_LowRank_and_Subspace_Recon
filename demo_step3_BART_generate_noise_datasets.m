% demo_step3_BART_generate_noise_datasets.m
% Written by Kyungmin Nam
% Email: K.M.Nam@umcutrecht.nl
% Started: 09/16/2024, Last modified: 

%% Clean slate
close all; clear variables; clc;

%% Set source paths
package_path = 'E:\mfiles_kmin\MRM_dSPICE\demo\';

%% Add source paths to search path
addpath(genpath(package_path));

%% Define input files
data_path   = 'E:\mfiles_kmin\MRM_dSPICE\demo\data\';
input_path  = [data_path, 'input\'];
noise_path  = [data_path, 'output\noise\'];

%% Make an output path
mkdir(noise_path);

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
    bart_output_path = strrep(input_path, '\', '/');
    bart_output_path = sprintf('/mnt/%s/%s/', lower(bart_output_path(1)), bart_output_path(4:end));
else
    bart_output_path = sprintf('%s/', input_path);
end

%% Set the matrix size for a low-resolution dataset
N1_lowres = 12;
N2_lowres = 18;
N3_lowres = 15;

%% Define simulation parameters
spectral_bandwidth    = 5000;      % spectral bandwidth [Hz]
nr_acquisition_points = 1024;      % number of acquisition points
nr_spectral_points    = 1024;      % number of spectral points
B0                    = 7;         % main field strength [T]
gamma                 = 6.53569e6; % gyromagnetic ratio [Hz/T]

%% Calculate a time axis [sec]
time_axis = (0:nr_acquisition_points-1).' / spectral_bandwidth; % [sec]

%% Define metabolites
carrier_frequency_ppm = 4.7; % [ppm]

%% Calculate a ppm axis [ppm]
df = spectral_bandwidth / nr_spectral_points; % [Hz]
freq_axis = (-floor(nr_spectral_points/2):ceil(nr_spectral_points/2)-1).' * df; % [Hz]
ppm_axis = freq_axis / (gamma * B0 * 1e-6) + carrier_frequency_ppm; % [Hz] / [Hz/ppm] => [ppm]

%--------------------------------------------------------------------------
% ksp_lowres (N1 x N2 x N3 x 1 x 1 x M)
%--------------------------------------------------------------------------
ksp_filename = 'ksp_ori';

%% Read a .cfl file
%--------------------------------------------------------------------------
% ksp (N1 x N2 x N3 x Nc x 1 x M)
%--------------------------------------------------------------------------
cfl_file = fullfile(input_path, ksp_filename);
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
ksp = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
[N1,N2,N3,~,~,M] = size(ksp);

%% Calculate low-resolution signle-channel spectra
tstart = tic; fprintf('%s: Calculating low-resolution single-channel spectra... ', datetime);

%--------------------------------------------------------------------------
% image space (x-y-z) <=> k-space (kx-ky-kz)
%--------------------------------------------------------------------------
img_lowres = ksp;
for dim = 1:3
    img_lowres = sqrt(size(img_lowres,dim)) * ...
        fftshift(ifft(ifftshift(img_lowres, dim), [], dim), dim);
end

%--------------------------------------------------------------------------
% kf(=t) => f
%--------------------------------------------------------------------------
img_lowres = 1 / sqrt(size(img_lowres,6)) * fftshift(fft(img_lowres, [], 6), 6);

fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read a .cfl file
%--------------------------------------------------------------------------
% sens (N1 x N2 x N3 x Nc)
%--------------------------------------------------------------------------
cfl_file = fullfile(input_path, 'sens');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
sens = readcfl(cfl_file);
nr_channels = size(sens,4);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read a .mat file
%--------------------------------------------------------------------------
% mask_wavg (N1 x N2 x N3 x 1 x 1 x M x ... x Nc(15th))
%--------------------------------------------------------------------------
cfl_file = fullfile(input_path, 'mask_wavg');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
mask_wavg = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Add gaussian random noise
num_noise = 4;  % 4 noise level
itr_noise = 30; % 30 iteration
averages  = 4;  % weighted averages

%--------------------------------------------------------------------------
% Genreate SNR table (Add noise)
% 60.0477+-7.2879
% 30.4202+-3.8472
% 15.0509+-2.0765
% 05.1114+-1.1892
%--------------------------------------------------------------------------
noise_std = [0.923,1.830,3.700,10.900]; 
NC = 4;

for idx_nslvl = 3%1:num_noise
    for idx_itr = 1%1:itr_noise
        img_noise_ave = zeros(N1_lowres,N2_lowres,N3_lowres,NC,1,nr_spectral_points,1,1,1,1,1,1,1,1,averages);
        for idx_nc = 1:NC   
            for idx_nsa = 1:averages
                % add ground trunth with noise 
                img_noise_ = img_lowres + noise_std(idx_nslvl) * ...
                    (1/sqrt(2)*(randn(N1,N2,N3,1,1,M)+1i*randn(N1,N2,N3,1,1,M)));
                img_noise_ave(:,:,:,idx_nc,1,:,1,1,1,1,1,1,1,1,idx_nsa) = img_noise_;           
            end
        end
        %--------------------------------------------------------------------------
        % image space (x-y-z) => k-space (kx-ky-kz)
        %--------------------------------------------------------------------------
        ksp_noise = img_noise_ave;
        for dim = 1:3
            ksp_noise = 1/sqrt(size(ksp_noise,dim)) * ...
                fftshift(fft(ifftshift(ksp_noise, dim), [], dim), dim);
        end
        %--------------------------------------------------------------------------
        % kf(=t) <= f
        %--------------------------------------------------------------------------
        ksp_noise = sqrt(size(ksp_noise,6)) * ifft(ifftshift(ksp_noise,6),[],6);
            
        % Multiply hamming weighting acquisition mask
        ksp_noise = bsxfun(@times, mask_wavg, ksp_noise);

        % Calculate the weighted average of k-space
        mask_wavg_ = sum(mask_wavg,15);
        ksp_noise = sum(ksp_noise,15) ./ mask_wavg_;

        %--------------------------------------------------------------------------
        % image space (x-y-z) <= k-space (kx-ky-kz)
        %--------------------------------------------------------------------------
        img_noise = ksp_noise;
        for dim = 1:3
            img_noise = sqrt(size(img_noise,dim)) * ...
                fftshift(ifft(ifftshift(img_noise, dim), [], dim), dim);
        end
                
        % Calculate low-resolution multi-channel noise spectra (x-y-z-kf(=t))
        img_chans_noise = bsxfun(@times, sens, img_noise);               
        
        %--------------------------------------------------------------------------
        %image space (x-y-z) <=> k-space (kx-ky-kz)
        %--------------------------------------------------------------------------
        ksp_chans_noise = img_chans_noise;
        for dim = 1:3
            ksp_chans_noise = 1 / sqrt(size(ksp_chans_noise,dim)) * ...
                fftshift(fft(ifftshift(ksp_chans_noise, dim), [], dim), dim);
        end
        
        %% Write a .cfl file
        %--------------------------------------------------------------------------
        %ksp_chans_noise_nslvl%d_itr%d (N1 x N2 x N3 x Nc x 1 x M)
        %--------------------------------------------------------------------------
        ksp_filename = sprintf('ksp_chans_nslvl%d_nsitr%d',idx_nslvl,idx_itr);
        cfl_file = fullfile(noise_path, ksp_filename);
        tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
        writecfl(cfl_file, ksp_chans_noise);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

        %% Create a sampling mask
        ksp_file = strcat(bart_output_path, ksp_filename);  % N1 x N2 x N3 x Nc x 1 x M
        pat_file = strcat(bart_output_path, 'mask');        % N1 x N2 x N3 x  1 x 1 x M

        command = sprintf('%s pattern %s %s', bart_command, ksp_file, pat_file);
        tstart = tic; fprintf('%s:[BART] Creating a sampling mask:\n%s\n', datetime, command);
        [status_pattern,result_pattern] = system(command);
        fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));
    
        %% Read a .cfl file
        %--------------------------------------------------------------------------
        % ksp (N1 x N2 x N3 x Nc x 1 x M)
        %--------------------------------------------------------------------------
        cfl_file = fullfile(noise_path, ksp_filename);
        tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
        ksp_chans_noise = readcfl(cfl_file);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
        [N1,N2,N3,Nc,~,M] = size(ksp_chans_noise);

        %% Calculate low-resolution multi-channel spectra
        tstart = tic; fprintf('%s: Calculating low-resolution multi-channel spectra... ', datetime);

        %--------------------------------------------------------------------------
        % image space (x-y-z) <=> k-space (kx-ky-kz)
        %--------------------------------------------------------------------------
        img_chans_noise = ksp_chans_noise;
        for dim = 1:3
            img_chans_noise = sqrt(size(img_chans_noise,dim)) * ...
                fftshift(ifft(ifftshift(img_chans_noise, dim), [], dim), dim);
        end
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
        
        %% Perform optimal coil combination
        img_chcomb_noise = sum(bsxfun(@times, conj(sens), img_chans_noise), 4); % x-y-z-f
        
        %% Calculate low-resolution multi-channel k-space data (kx-ky-kz-kf(=t))
        %--------------------------------------------------------------------------
        % BART: image space (x-y-z) <=> k-space (kx-ky-kz)
        %--------------------------------------------------------------------------
        ksp_chcomb_noise = img_chcomb_noise;
        for dim = 1:3
            ksp_chcomb_noise = 1 / sqrt(size(ksp_chcomb_noise,dim)) * ...
                fftshift(fft(ifftshift(ksp_chcomb_noise, dim), [], dim), dim);
        end    
        %% Write a .cfl file
        %--------------------------------------------------------------------------
        % ksp_lowres (N1 x N2 x N3 x Nc x 1 x M)
        %--------------------------------------------------------------------------
        ksp_filename = sprintf('ksp_chcomb_nslvl%d_nsitr%d', idx_nslvl, idx_itr);
        cfl_file = fullfile(noise_path, ksp_filename);
        tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
        writecfl(cfl_file, ksp_chcomb_noise);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
        
        %--------------------------------------------------------------------------
        % kf(=t) => f
        %--------------------------------------------------------------------------
        img_chcomb_noise = 1 / sqrt(size(img_chcomb_noise,6)) * fftshift(fft(img_chcomb_noise, [], 6), 6);
            
        %% Display coil-combined spectra
        max_val = max(abs(img_chcomb_noise(:)));

        img_display = squeeze(img_chcomb_noise);

        c1 = floor(N1_lowres/2) + 1;
        c2 = floor(N2_lowres/2) + 1;
        c3 = floor(N3_lowres/2) + 1;

        figure('color', 'w', 'Position', [1063 50 855 513]);
        ColorOrder = get(gca, 'ColorOrder');
        for idx1 = 1:N1_lowres
            for idx2 = 1:N2_lowres
                spectrum = squeeze(real(img_display(idx1,idx2,c3,:)));
                hAxe = axes('Parent', gcf, 'linewidth', 1, 'XTick', [], 'YTick', [], ...
                    'Box', 'on', 'Position', [(idx2 - 1) / N2_lowres, 1 - (idx1 / N1_lowres), 1 / N2_lowres, 1 / N1_lowres], ...
                    'XDir', 'reverse', 'XTickLabel', {''; ''}, 'YTickLabel', {''; ''}, ...
                    'TickDir', 'out', 'YLim', [-0.2 * max_val max_val], 'XLim', [0 10]);
                hLine = line('Parent', hAxe, 'Linewidth', 0.5, 'XData', ppm_axis, 'YData', spectrum, 'Color', ColorOrder(1,:));
            end
        end
    end
end