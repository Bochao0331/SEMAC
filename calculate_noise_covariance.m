function [Psi,inv_L] = calculate_noise_covariance(noise)
%   Calculae noise covariance marix from noise samples
%   Input:
%       - noise: Ncoil x Nsamples x Nrepetitions

%   Output:
%       - Psi: scaled noise coavariance matrix (Ncoil x Ncoil)
%       - inv_L: inverse sqaure roof of the noise covariance (Ncoil x Ncoil)

%   Reference:
%   1. Kellman P, McVeigh ER. Magn Reson Med. 2005;54(6):1439-1447.

%   Author: Bochao Li
%   Email: bochaoli@usc.edu

[nr_channels,nr_samples,nr_repetitions] = size(noise);
tic; fprintf('Calculating the receiver noise matrix... ');
Nc = nr_channels;
Ns = nr_samples * nr_repetitions;
eta = reshape(noise, [Nc Ns]);
Psi = eta * eta' / Ns; % Equation B1
fprintf('done! (%6.4f sec)\n', toc);

% ---- covairaince matrix scaled by noise equivalent bandwidth (Eq[2]) ---- 
tic; fprintf('Calculating noise equivalent bandwidth... ');
noise_vec = permute(reshape(noise, [Nc Ns]),[2,1]);
for k =1:size(noise_vec,2)
    [p(:,k),f(:,k)] = pspectrum(noise_vec(:,k).','FrequencyLimits',[-pi pi]);
end
flat_mag_noise = mean(p(910:3115,:),1); % flat portion of spectrum
all_mag_noise = mean(p,1);
noise_ratio = mean(all_mag_noise./flat_mag_noise);
Psi = Psi/noise_ratio;
fprintf('done! (%6.4f sec)\n', toc);

%% Calculate the Cholesky decomposition of the receiver noise matrix
tic; fprintf('Calculating the Cholesky decomposition of the receiver noise matrix... ');
L = chol(Psi, 'lower'); % Equation B4
inv_L = sqrt(2)*inv(L);
fprintf('done! (%6.4f sec)\n', toc);
end