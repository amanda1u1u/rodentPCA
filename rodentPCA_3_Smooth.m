% Smoothing by Gaussian kernels, std = 100 ms
% Author:  Zhouxiao Lu
% Date: Mar. 28, 2022
% Last modified on: Apr. 24, 2022

clear; clc; close all;
p1 = fullfile('toolbox');
addpath(genpath(p1));
DataFolder = '../Datasets';

% Open file dialog box to select a .mat file
[file,path] = uigetfile('*.mat','Select the .mat file');
if isequal(file,0)
    disp('User selected Cancel');
else
    % Load the selected .mat file
    load(fullfile(path,file));  
end

%%
% Define Gaussian kernel parameters
sigma = 2; % *50ms = standard deviation
window = 5*sigma; % window size
x = -window:window; % domain of the kernel
kernel = exp(-x.^2/(2*sigma^2)) / (sigma*sqrt(2*pi)); % Gaussian kernel

% Smooth spike firing counts with Gaussian kernel
for i = 1:size(Data.trial_neuron_timestamps,1)
    for j = 1:size(Data.trial_neuron_timestamps,2)
        spike_counts = Data.binned{i,j};
        smooth_counts = conv(spike_counts, kernel, 'same');
        Data.smoothed{i,j} = smooth_counts;
    end
end
%%
% Plot original and smoothed spike firing counts
figure;
plot(Data.binned{1,1}, 'ko-', 'LineWidth', 1.5);
hold on;
plot(Data.smoothed{1,1}, 'r-', 'LineWidth', 1.5);
legend('Original', 'Smoothed');
xlabel('Time (ms)');
ylabel('Spike Count');
title('Spike Firing Counts, trial 1');

%% Save smoothed data
mat_name = fullfile(path, file);
save(mat_name,'Data')