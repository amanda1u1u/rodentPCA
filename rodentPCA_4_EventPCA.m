% PCA: Event-wise
% Author:  Zhouxiao Lu
% Date: Apr. 24, 2022
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
% Define time window for event
bin_size = Data.bin_size;
extra_time = Data.ExtraMarginTime;
t_peri = Data.PeriTime;
t_start = extra_time/bin_size;
t_end = t_start+(t_peri/bin_size);

for i = 1:size(Data.trial_neuron_timestamps,1)
    for j = 1:size(Data.trial_neuron_timestamps,2)
        trial_data = Data.smoothed{i,j};
        n = length(trial_data);
        sample_data = trial_data(t_start:t_end);
        response_data = trial_data(n-t_end:n-t_start);
        Data.sample_neuron{i,j} = sample_data;
        Data.response_neuron{i,j} = response_data;
    end
end

%% Groups
groups = {'LEFT','SAMPLE','SUCCESS';'LEFT','SAMPLE','FAILURE';
    'LEFT','RESPONSE','SUCCESS';'LEFT','RESPONSE','FAILURE';
    'RIGHT','SAMPLE','SUCCESS';'RIGHT','SAMPLE','FAILURE';
    'RIGHT','RESPONSE','SUCCESS';'RIGHT','RESPONSE','FAILURE'};

for i = 1:length(groups)
    str_position = groups(i,1);
    str_phase = groups(i,2);
    str_type = groups(i,3);
    if strcmp(str_phase,'SAMPLE');
        phase = Data.SamplePosition;
        data = Data.sample_neuron;
    else strcmp(str_phase,'RESPONSE');
        phase = Data.ResponsePosition;
        data = Data.response_neuron;
    end
    type = Data.TrialType;
    ind_position = find(strcmp(phase,str_position));
    ind_type = find(strcmp(type,str_type));
    ind = intersect(ind_position,ind_type);

    m = zeros(size(Data.smoothed,1),t_end-t_start+1);
    for j = 1:length(ind)
        k = ind(j);
        d = cell2mat(data(:,k));
        m = m + d;
    end
    Data.group_neuron{i} = m/j; % average among events for each group
end
Data.groups = groups;

%% Center each average spike time series vector by subtracting its mean.
for i = 1:length(groups)
    centered_group{i} = Data.group_neuron{1,i} - mean(Data.group_neuron{1,i});
end

%% covariance
covariances = cell(1, length(groups));
for i = 1:length(groups)
    covariances{i} = cov(centered_group{i}');
end

%%
% Compute eigenvectors and eigenvalues for each covariance matrix
eigenvectors = cell(1, length(groups));
eigenvalues = cell(1, length(groups));
for i = 1:length(groups)
    [eigenvectors{i}, eigenvalues{i}] = eig(covariances{i});
end

% Sort eigenvectors by descending eigenvalues for each behavior type
for i = 1:length(groups)
    [sorted_eigenvalues{i}, idx] = sort(diag(eigenvalues{i}), 'descend');
    var_explained{i} = sorted_eigenvalues{i}/sum(sorted_eigenvalues{i});
    var_explained{i} = cumsum(var_explained{i});
    eigenvectors{i} = eigenvectors{i}(:, idx);
    % find top k PCs that explain more than 90% of variance
    id = find(var_explained{i}>=0.9);
    idx_var{i} = id(1);
end

% Select top k eigenvectors
for i = 1:length(groups)
    k = idx_var{i};
    projection = eigenvectors{i}(:,1:k);
    % Project spike time series matrix onto projection matrix
    spike_reduced = centered_group{i}' * projection;
    Data.pca_event{i} = spike_reduced;
end

mat_name = fullfile(path, file);
save(mat_name,'Data')

