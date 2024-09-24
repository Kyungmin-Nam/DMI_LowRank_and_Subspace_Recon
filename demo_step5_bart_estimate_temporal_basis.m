% demo_step5_bart_estimate_temporal_basis.m
% Written by Nam G. Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Modified by Kyungmin Nam
% Email: K.M.Nam@umcutrecht.nl
% Started: 06/01/2024, Last modified: 09/24/2024

%% Clean slate
close all; clear variables; clc;

%% Set source paths
package_path = 'E:\mfiles_kmin\MRM_dSPICE\demo\thirdparty';

%% Add source paths to search path
addpath(genpath(package_path));

%% Define input files
data_path   = 'E:\mfiles_kmin\MRM_dSPICE\demo';
input_path  = [data_path, '\data\input\'];
output_path = [data_path, '\data\output\'];
basis_path  = [output_path, 'basis\']

%% Define a path to BART
bart_path = '/home/kmin/bart'; % Please install WSL in your MS Windows

%% Start a stopwatch timer
start_time = tic;

%% Make an output path
mkdir(output_path);

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
% sens (N1 x N2 x N3 x Nc)
%--------------------------------------------------------------------------
cfl_file = fullfile(input_path, 'sens');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
sens = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Apply forward Fourier transforms (spatial domain only)
%--------------------------------------------------------------------------
% BART: image space (x-y-z) => k-space (kx-ky-kz)
%--------------------------------------------------------------------------
ksp_sens = sens;
for dim = 1:3
    ksp_sens = 1 / sqrt(size(ksp_sens,dim)) * fftshift(fft(ifftshift(ksp_sens, dim), [], dim), dim);
end

%%
D1_matrix = [3 5 3]; % D1 size
num_noise = 4;       % 4 noise level
itr_noise = 30;      % 30 iteration

for idx_D1size = 1
for idx_nslvl = 3%:num_noise
for idx_nsitr = 1%:itr_noise

%% Define the size of navigator data (kx-ky-kz-kf(=t))
D1_size = D1_matrix(idx_D1size,:);

%--------------------------------------------------------------------------
% ksp (N1 x N2 x N3 x Nc x 1 x M)
%--------------------------------------------------------------------------
cfl_file = fullfile([output_path,'noise\'], sprintf('ksp_chans_nslvl%d_nsitr%d',idx_nslvl,idx_nsitr));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
ksp = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
[N1,N2,N3,Nc,~,M] = size(ksp);

%% Calculate low-resolution coil sensitivity maps
idx1_range = (-floor(D1_size(1)/2):ceil(D1_size(1)/2)-1).' + floor(N1/2) + 1;
idx2_range = (-floor(D1_size(2)/2):ceil(D1_size(2)/2)-1).' + floor(N2/2) + 1;
idx3_range = (-floor(D1_size(3)/2):ceil(D1_size(3)/2)-1).' + floor(N3/2) + 1;

ksp_sens1 = ksp_sens(idx1_range, idx2_range, idx3_range, :);

%% Apply inverse Fourier transforms (spatial domain only)
%--------------------------------------------------------------------------
% BART: image space (x-y-z) <=> k-space (kx-ky-kz)
%--------------------------------------------------------------------------
sens1 = ksp_sens1;
for dim = 1:3
    sens1 = sqrt(size(sens1,dim)) * fftshift(ifft(ifftshift(sens1, dim), [], dim), dim);
end

%% Normalize low-resolution coil sensitivity maps
%sens1 = sens1 ./ sqrt(sum(abs(sens1).^2,4));

%% Extract a limited region of central k-space
idx1_range = (-floor(D1_size(1)/2):ceil(D1_size(1)/2)-1).' + floor(N1/2) + 1;
idx2_range = (-floor(D1_size(2)/2):ceil(D1_size(2)/2)-1).' + floor(N2/2) + 1;
idx3_range = (-floor(D1_size(3)/2):ceil(D1_size(3)/2)-1).' + floor(N3/2) + 1;

ksp1 = ksp(idx1_range, idx2_range, idx3_range, :, 1, :);

%% Apply inverse Fourier transforms (spatial domain only)
%--------------------------------------------------------------------------
% BART: image space (x-y-z) <= k-space (kx-ky-kz)
%--------------------------------------------------------------------------
cimg1 = ksp1;
for dim = 1:3
    cimg1 = sqrt(size(cimg1,dim)) * fftshift(ifft(ifftshift(cimg1, dim), [], dim), dim);
end

%% Perform optimal coil combination
img1 = sum(bsxfun(@times, conj(sens1), cimg1), 4); % x-y-z-kf(=t)

%% Define simulation parameters
spectral_bandwidth    = 5000;      % spectral bandwidth [Hz]
nr_acquisition_points = 1024;      % number of acquisition points
nr_spectral_points    = 1024;      % number of spectral points
B0                    = 7;         % main field strength [T]
gamma                 = 6.53569e6; % gyromagnetic ratio [Hz/T]

%% Calculate a ppm axis [ppm]
carrier_frequency_ppm = 4.7; % [ppm]

df = spectral_bandwidth / nr_spectral_points; % [Hz]
freq_axis = (-floor(nr_spectral_points/2):ceil(nr_spectral_points/2)-1).' * df; % [Hz]
ppm_axis = freq_axis / (gamma * B0 * 1e-6) + carrier_frequency_ppm; % [Hz] / [Hz/ppm] => [ppm]

%% Calculate the Casorati matrix of dataset "D1"
Casorati = reshape(img1, [prod(D1_size) M]);

%% Perform the SVD of navigator data
tstart = tic; fprintf('%s: Computing the SVD of D... ', datetime);
[U,S,V] = svd(Casorati,0); % U: N x N, S: N x M, V: M x M
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Set the maximum rank of a temporal subspace
Lmax = min(size(Casorati));

%% Write a .cfl file
for L = 4:Lmax
    %----------------------------------------------------------------------
    % basis (1 x 1 x 1 x 1 x 1 x M x L)
    %----------------------------------------------------------------------
    cfl_file = fullfile(basis_path, sprintf('basis_%dx%dx%d_L%d_nslvl%d_nsitr%d',...
        D1_size(1), D1_size(2), D1_size(3), L, idx_nslvl, idx_nsitr));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, reshape(conj(V(:,1:L)), [1 1 1 1 1 M L]));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

if 1
    %% Display for sanity check
    figure('Color', 'w', 'Position', [5 349 1917 629]); 
    subplot(2,2,1);
    imagesc(abs(Casorati));
    axis image; title('Casorati matrix', 'FontWeight', 'normal');

    subplot(2,2,3);
    imagesc(abs(V(:,1:min(size(V)))'));
    axis image; title('Right singular vectors, V^H', 'FontWeight', 'normal');

    subplot(2,2,2);
    plot(diag(S) / max(S(:)));
    title('Singular values', 'FontWeight', 'normal'); grid on; grid minor;

    subplot(2,2,4);
    plot(abs(V(:,1:L)));
    title(sprintf('%d principal right singular vectors', L), 'FontWeight', 'normal');
    grid on; grid minor;

    %% Save figure
    exportgraphics(gcf, fullfile(output_path, ...
        sprintf('subspace_estimation_%dx%dx%d_nslvl%d_nsitr%d.tif', ...
        D1_size(1), D1_size(2), D1_size(3), idx_nslvl, idx_nsitr)));

end
end
end
end

