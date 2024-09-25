% demo_step6_bart_pics_lowrank_subspace_v1.m
% Initialized by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Modified by Kyungmin Nam at UMC Utrecht
% Email: K.M.Nam@umcutrecht.nl
% Started: 06/02/2024, Last modified: 09/25/2024

%% Clean slate
close all; clear variables; clc;

%% Set source paths
package_path = 'E:\mfiles_kmin\MRM_dSPICE\demo';

%% Add source paths to search path
addpath(genpath(package_path));

%% Define input files
data_path   = 'E:\mfiles_kmin\MRM_dSPICE\demo\';
sens_path   = [data_path,   'data\input\'];
input_path  = [data_path,   'data\output\noise\'];
basis_path  = [data_path,   'data\output\basis\'];
output_path = [input_path,  'LRSM\'];
mask_path   = [output_path, 'mask\'];
debug_path  = [output_path, 'debug_msg\'];
figure_path = [output_path, 'figs\'];

%% Make output path
mkdir(output_path);
mkdir(debug_path);
mkdir(mask_path);
mkdir(figure_path);

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
    bart_input_path = strrep(input_path, '\', '/');
    bart_input_path = sprintf('/mnt/%s/%s/', lower(bart_input_path(1)), bart_input_path(4:end));   

    bart_sens_path = strrep(sens_path, '\', '/');
    bart_sens_path = sprintf('/mnt/%s/%s/', lower(bart_sens_path(1)), bart_sens_path(4:end));   
    
    bart_basis_path = strrep(basis_path, '\', '/');
    bart_basis_path = sprintf('/mnt/%s/%s/', lower(bart_basis_path(1)), bart_basis_path(4:end));   

    bart_mask_path = strrep(mask_path, '\', '/');
    bart_mask_path = sprintf('/mnt/%s/%s/', lower(bart_mask_path(1)), bart_mask_path(4:end));   

    bart_output_path = strrep(output_path, '\', '/');
    bart_output_path = sprintf('/mnt/%s/%s/', lower(bart_output_path(1)), bart_output_path(4:end));   
else
    bart_input_path  = sprintf('%s/', input_path);
    bart_sens_path   = sprintf('%s/', sens_path);
    bart_basis_path  = sprintf('%s/', basis_path);
    bart_mask_path   = sprintf('%s/', basis_path);    
    bart_output_path = sprintf('%s/', output_path);
end

%% Define sequence parameters
maxiter = 30;
Lmax    = 10;

%% Calculate a ppm axis [ppm]
spectral_bandwidth    = 5000;      % spectral bandwidth [Hz]
nr_acquisition_points = 1024;      % number of acquisition points
nr_spectral_points    = 1024;      % number of spectral points
B0                    = 7;         % main field strength [T]
gamma                 = 6.53569e6; % gyromagnetic ratio [Hz/T]

carrier_frequency_ppm = 4.7; % [ppm]

spectral_resolution = spectral_bandwidth / nr_spectral_points; % [Hz]
freq_axis = (-floor(nr_spectral_points/2):ceil(nr_spectral_points/2)-1).' * spectral_resolution; % [Hz]
ppm_axis  = freq_axis / (gamma * B0 * 1e-6) + carrier_frequency_ppm; % [Hz] / [Hz/ppm] => [ppm]

%%
num_noise = 4;  % 4 noise level
itr_noise = 30; % 30 iteration
lambda    = [0.001 0.01 0.03 0.05];
D1_matrix = [3 5 3];
R         = [1.00 1.10 1.30];
%%

for idx_D1 = 1
for idx_lambda = 1%:length(lambda)    
for idx_nslvl = 3%:num_noise
for idx_nsitr = 1%:itr_noise
for idx_R = 1%:length(R)

    %% Define the size of navigator data (kx-ky-kz-kf(=t))
    D1_size = D1_matrix(idx_D1,:);

    %% Perform DMI SPICE reconstruction
for L = 5%4:Lmax
    %% Perform a low-rank and subspace model-based reconstruction
    pat_file = strcat(bart_mask_path, ...
    sprintf('vd_poisson_mask_%dx%dx%d_R%1.2f', D1_size(1),D1_size(2),D1_size(3),R(idx_R)));  % N1 x N2 x N3 x  1 x 1 x M
    ksp_file = strcat(bart_input_path, sprintf('ksp_chans_nslvl%d_nsitr%d',idx_nslvl,idx_nsitr));   % N1 x N2 x N3 x Nc x 1 x M
    basis_file = strcat(bart_basis_path, sprintf('basis_%dx%dx%d_L%d_nslvl%d_nsitr%d',...
        D1_size(1), D1_size(2), D1_size(3), L, idx_nslvl, idx_nsitr));  % 1 x 1 x 1 x 1 x 1 x M (TE) x L (COEFF)
    sens_file = strcat(bart_sens_path, 'sens'); % N1 x N2 x N3 x Nc

    cfimg_filename = sprintf('cfimg_tv_i%d_l%g_%dx%dx%d_L%d_nslvl%d_nsitr%d_R%1.1f', ...
        maxiter, lambda(idx_lambda), D1_size(1), D1_size(2), D1_size(3), L, idx_nslvl, idx_nsitr, R(idx_R));
    img_filename   = sprintf('img_tv_i%d_l%g_%dx%dx%d_L%d_nslvl%d_nsitr%d_R%1.1f', ...
        maxiter, lambda(idx_lambda), D1_size(1), D1_size(2), D1_size(3), L, idx_nslvl, idx_nsitr, R(idx_R));

    cfimg_file = strcat(bart_output_path, cfimg_filename); % N1 x N2 x N3 x  1 x 1 x 1 x L

    command = sprintf('%s pics -e -d5 -S -i%d -R T:7:0:%d -p %s -B %s %s %s %s', ...
        bart_command, maxiter, lambda(idx_lambda), pat_file, basis_file, ksp_file, sens_file, cfimg_file);

    tstart = tic; fprintf('%s:(L=%d/%d)[BART] Performing DMI SPICE reconstruction:\n%s\n',...
        datetime, L, Lmax, command);
    [status_pics,result_pics] = system(command);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %% Save the output as a .txt file
    txt_file = fullfile(debug_path, sprintf('%s_debug.txt', cfimg_filename));
    [fid,message] = fopen(txt_file, 'w');
    fwrite(fid, result_pics, 'char');
    fclose(fid);

    if 1
        %% Read a .cfl file
        %----------------------------------------------------------------------
        % cfimg (N1 x N2 x N3 x 1 x 1 x 1 x L) (x-y-z-kf(=t))
        %----------------------------------------------------------------------
        cfl_file = fullfile(output_path, cfimg_filename);
        tstart = tic; fprintf('%s:(L=%d/%d) Reading a .cfl file: %s... ', datetime, L, Lmax, cfl_file);
        cfimg = readcfl(cfl_file);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
        [N1,N2,N3,~,~,~,L] = size(cfimg);

        %% Display spatial coefficients
        cfimg_montage = complex(zeros(N1 * N3, N2 * L, 'single'));
        for idx1 = 1:N3
            idx1_range = (1:N1).' + (idx1 - 1) * N1;
            cfimg_montage(idx1_range,:) = reshape(cfimg(:,:,idx1,1,1,1,:), [N1 N2 * L]);
        end

        figure('Color', 'w', 'Position', [91 14 526 964]);
        imagesc(abs(cfimg_montage));
        axis image;
        saveas(gcf, fullfile(figure_path, [cfimg_filename,'.fig']));
        close gcf;
    end

    %% Multiply the temporal basis with the spatial coefficients to reconstrct the Casorati matrix
    img_file = strcat(bart_output_path, img_filename);

    % bart bitmask 6 => 64
    command = sprintf('%s fmac -s 64 %s %s %s', bart_command, basis_file, cfimg_file, img_file);
    tstart = tic; fprintf('%s:(L=%d/%d)[BART] Multiplying the temporal basis with the spatial coefficients:\n%s\n', datetime, L, Lmax, command);
    [status_fmac,result_fmac] = system(command);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    if 1
        %% Read a .cfl file
        %----------------------------------------------------------------------
        % img (N1 x N2 x N3 x 1 x 1 x M) (x-y-z-kf(=t))
        %----------------------------------------------------------------------
        cfl_file = fullfile(output_path, img_filename);
        tstart = tic; fprintf('%s:(L=%d/%d) Reading a .cfl file: %s... ', datetime, L, Lmax, cfl_file);
        img = readcfl(cfl_file);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

        %% Apply a 1D forward Fourier transform along the CSHIFT dimension
        %----------------------------------------------------------------------
        % kf(=t) => f
        %----------------------------------------------------------------------
        img = 1 / sqrt(size(img,6)) * fftshift(fft(img, [], 6), 6);
        
        %% Display spectra (x-y-z-f)
        img_display = squeeze(img) .* exp(1i*pi/2);

        [N1_,N2_,N3_,~] = size(img_display);

        max_val = max(abs(img_display(:)));

        c1 = floor(N1_/2) + 1;
        c2 = floor(N2_/2) + 1;
        c3 = floor(N3_/2) + 1;

        figure('color', 'w', 'Position', [1063 2 855 513]);
        ColorOrder = get(gca, 'ColorOrder');
        for idx1 = 1:N1_
            for idx2 = 1:N2_
                spectrum = squeeze(real(img_display(idx1,idx2,c3,:)));
                hAxe = axes('Parent', gcf, 'linewidth', 1, 'XTick', [], 'YTick', [], ...
                    'Box', 'on', 'Position', [(idx2 - 1) / N2_, 1 - (idx1 / N1_), 1 / N2_, 1 / N1_], ...
                    'XDir', 'reverse', 'XTickLabel', {''; ''}, 'YTickLabel', {''; ''}, ...
                    'TickDir', 'out', 'YLim', [-0.2 * max_val max_val], 'XLim', [0 10]);
                hLine = line('Parent', hAxe, 'Linewidth', 0.5, 'XData', ppm_axis, 'YData', spectrum, 'Color', ColorOrder(1,:));
            end
        end
        saveas(gcf, fullfile(figure_path, [img_filename,'.tif']));
        close gcf;       
    end
    
end    
end
end
end
end
end
%%

