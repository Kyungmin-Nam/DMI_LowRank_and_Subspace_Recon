%--------------------------------------------------------------------------
% demo_display_DMI_spectra.m
% Written by Kyungmin Nam
% Email: K.M.Nam@umcutrecht.nl
% Started: [mm/dd/yyyy] 09/16/2024, Last modified: 
% Inspired to use the BART Toolbox by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
%--------------------------------------------------------------------------

%% Clean slate
close all; clear variables; clc;

%% Set source paths
package_path = 'E:\mfiles_kmin\MRM_dSPICE\demo\thirdparty\';

%% Add source paths to search path
addpath(genpath(package_path));

%% Define input files
data_path   = 'E:\mfiles_kmin\MRM_dSPICE\demo\data\';
input_path  = [data_path, 'input\'];
output_path = [data_path, 'output'];

%% Make an output path
mkdir(output_path);

%% Define a path to BART (i.e., WSL)
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

%% Read a .cfl file
%--------------------------------------------------------------------------
% ksp_ori (N1 x N2 x N3 x 1 x 1 x M) (kx-ky-kz-1-1-kf(=t))
%--------------------------------------------------------------------------
ksp_filename = 'ksp_ori';
cfl_file = fullfile(input_path, ksp_filename);
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
ksp_lowres = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Calculate low-resolution spectra (x-y-z-kf(=t))
%--------------------------------------------------------------------------
% BART: image space (x-y-z) <= k-space (kx-ky-kz)
%--------------------------------------------------------------------------
img_lowres = ksp_lowres;
for dim = 1:3
    img_lowres = sqrt(size(img_lowres,dim)) * fftshift(ifft(ifftshift(img_lowres, dim), [], dim), dim);
end

%% Display ground truth spectra (x-y-z-f)
img_display = squeeze(img_lowres);

%--------------------------------------------------------------------------
% kf(=t) => f
%--------------------------------------------------------------------------
img_display = 1/ sqrt(size(img_display, 4)) * fftshift(fft(img_display, [], 4), 4);

%% Figure 1B
img_low_res_map = sum(abs(img_display),4); 
figure('color', 'w', 'Position', [600 50 100 100]);
imagesc(img_low_res_map(:,:,12)); 
axis tight; axis equal; xticks([]); yticks([]);

figure('color', 'w', 'Position', [800 50 100 100]);
imagesc(img_low_res_map(:,:,2)); 
axis tight; axis equal; xticks([]); yticks([]);

figure('color', 'w', 'Position', [1000 50 100 100]);
imagesc(img_low_res_map(:,:,1)); 
axis tight; axis equal; xticks([]); yticks([]);

%% Define simulation parameters
spectral_bandwidth    = 5000;      % spectral bandwidth [Hz]
nr_acquisition_points = 1024;      % number of acquisition points
nr_spectral_points    = 1024;      % number of spectral points
B0                    = 7;         % main field strength [T]
gamma                 = 6.53569e6; % gyromagnetic ratio [Hz/T]
carrier_frequency_ppm = 4.7;       % [ppm]

%% Calculate a time axis [sec]
time_axis = (0:nr_acquisition_points-1).' / spectral_bandwidth; % [sec]

%% Calculate a ppm axis [ppm]
df = spectral_bandwidth / nr_spectral_points; % [Hz]
freq_axis = (-floor(nr_spectral_points/2):ceil(nr_spectral_points/2)-1).' * df; % [Hz]
ppm_axis = freq_axis / (gamma * B0 * 1e-6) + carrier_frequency_ppm; % [Hz] / [Hz/ppm] => [ppm]

%% Figure 1C
max_val = max(real(img_display(:)));
[N1_lowres,N2_lowres,N3_lowres,~] = size(img_display);
c1 = floor(N1_lowres/2) + 1;
c2 = floor(N2_lowres/2) + 1;
c3 = 1;

figure('color', 'w', 'Position', [1063 50 855 513]);
ColorOrder = get(gca, 'ColorOrder');
for idx1 = 1:N1_lowres
    for idx2 = 1:N2_lowres
        spectrum = squeeze(real(img_display(idx1,idx2,c3,:)));
        hAxe = axes('Parent', gcf, 'linewidth', 1, 'XTick', [], 'YTick', [], 'Box', 'on', ...
            'Position', [(idx2 - 1) / N2_lowres, 1 - (idx1 / N1_lowres), 1 / N2_lowres, 1 / N1_lowres], ...
            'XDir', 'reverse', 'XTickLabel', {''; ''}, 'YTickLabel', {''; ''}, ...
            'TickDir', 'out', 'YLim', [-0.2 * max_val max_val], 'XLim', [0 10]);
        hLine = line('Parent', hAxe, 'Linewidth', 1.5, 'XData', ppm_axis, 'YData', spectrum,...
            'Color', 'k');
    end
end

%% Figure 1D, Liver
figure('color', 'w', 'Position', [1063 50 855 513]);
cnt = 0;
for idx_row = 1:3
    for idx_col = 1:3
        cnt = cnt + 1;
        subplot(3,3,cnt); plot(ppm_axis, squeeze(real(img_display(4+idx_row,5+idx_col,1,:)./max_val)));
        xlim([0 10]); ylim([-0.2 1.2]); box off; 
        set(gca,'Xdir','reverse'); grid on;
    end
end

%% Figure 1E, Stomach
figure('color', 'w', 'Position', [1063 50 300 513]);
cnt = 0;
for idx_row = 1:2
    for idx_col = 1
        cnt = cnt + 1;
        subplot(2,1,cnt); plot(ppm_axis, squeeze(real(img_display(5+idx_row,12+idx_col,1,:)./max_val)));
        xlim([0 10]); ylim([-0.2 1.2]); box off; 
        set(gca,'Xdir','reverse'); grid on;
    end
end

%% Figure 1F, Body boundary
figure('color', 'w', 'Position', [1063 50 855 513]);
cnt = 0;
for idx_row = 1:2
    for idx_col = 1:2
        cnt = cnt + 1;
        subplot(2,2,cnt); plot(ppm_axis, squeeze(real(img_display(9+idx_row,9+idx_col,1,:)./max_val)));
        xlim([0 10]); ylim([-0.2 1.2]); box off; 
        set(gca,'Xdir','reverse'); grid on;
    end
end

%% Read a .cfl file
%--------------------------------------------------------------------------
% sens (N1 x N2 x N3 x Nc)
%--------------------------------------------------------------------------
cfl_file = fullfile(input_path, 'sens');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
sens = readcfl(cfl_file);
nr_channels = size(sens,4);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Calculate low-resolution multi-channel spectra (x-y-z-kf(=t))
cimg_lowres = bsxfun(@times, sens, img_lowres);

%% Calculate low-resolution multi-channel k-space data (kx-ky-kz-kf(=t))
%--------------------------------------------------------------------------
% BART: image space (x-y-z) <=> k-space (kx-ky-kz)
%--------------------------------------------------------------------------
ksp_lowres = cimg_lowres;
for dim = 1:3
    ksp_lowres = 1 / sqrt(size(ksp_lowres,dim)) * fftshift(fft(ifftshift(ksp_lowres, dim), [], dim), dim);
end

%% Write a .cfl file
%--------------------------------------------------------------------------
% ksp_lowres (N1 x N2 x N3 x Nc x 1 x M)
%--------------------------------------------------------------------------
ksp_filename = sprintf('ksp_%dx%dx%d', N1_lowres, N2_lowres, N3_lowres);
cfl_file = fullfile(output_path, ksp_filename);
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, ksp_lowres);
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
cfl_file = fullfile(output_path, ksp_filename);
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
ksp_lowres = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
[N1,N2,N3,Nc,~,M] = size(ksp_lowres);

%--------------------------------------------------------------------------
% mask (N1 x N2 x N3 x Nc x 1 x M)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'mask');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
mask = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Calculate low-resolution multi-channel spectra
tstart = tic; fprintf('%s: Calculating low-resolution multi-channel spectra... ', datetime);

%--------------------------------------------------------------------------
% image space (x-y-z) <=> k-space (kx-ky-kz)
%--------------------------------------------------------------------------
cimg_lowres = ksp_lowres;
for dim = 1:3
    cimg_lowres = sqrt(size(cimg_lowres,dim)) * fftshift(ifft(ifftshift(cimg_lowres, dim), [], dim), dim);
end

fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Perform optimal coil combination
img_lowres = sum(bsxfun(@times, conj(sens), cimg_lowres), 4); % x-y-z-f

%--------------------------------------------------------------------------
% kf(=t) => f
%--------------------------------------------------------------------------
img = 1 / sqrt(size(img_lowres,6)) * fftshift(fft(img_lowres, [], 6), 6);

%% sanity check
figure('color','w'); plot(ppm_axis, squeeze(real(img(7,7,7,1,1,:))), 'Color', 'k'); 
box off; xlim([0 10]); set(gca, 'Xdir', 'reverse');

%% Write a .cfl file
%--------------------------------------------------------------------------
% img_lowres (N1 x N2 x N3 x 1 x 1 x M)
%--------------------------------------------------------------------------
img_filename = sprintf('img_%dx%dx%d', N1_lowres, N2_lowres, N3_lowres);
cfl_file = fullfile(output_path, img_filename);
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, img_lowres);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% save
save(fullfile(output_path, 'img_12x18x15.mat'), 'img_lowres', 'ppm_axis');

%% Display coil-combined spectra
max_val = max(abs(img(:)));

img_display = squeeze(img_lowres);

c1 = floor(N1_lowres/2) + 1;
c2 = floor(N2_lowres/2) + 1;
c3 = floor(N3_lowres/2) + 1;

figure('color', 'w', 'Position', [1063 50 855 513]);
for idx1 = 1:N1_lowres
    for idx2 = 1:N2_lowres
        spectrum = squeeze(real(img_display(idx1,idx2,c3,:)));
        hAxe = axes('Parent', gcf, 'linewidth', 1, 'XTick', [], 'YTick', [], 'Box', 'on', ...
            'Position', [(idx2 - 1) / N2_lowres, 1 - (idx1 / N1_lowres), 1 / N2_lowres, 1 / N1_lowres], ...
            'XDir', 'reverse', 'XTickLabel', {''; ''}, 'YTickLabel', {''; ''}, ...
            'TickDir', 'out', 'YLim', [-0.2 * max_val max_val], 'XLim', [0 10]);
        hLine = line('Parent', hAxe, 'Linewidth', 0.5, 'XData', ppm_axis, 'YData', spectrum, 'Color', ColorOrder(1,:));
    end
end
