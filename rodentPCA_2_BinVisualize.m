% Binning spiking, bin size = 50ms
% This plots the raster plots for each trial in one DNMS session
% blue represents LEFT lever, red represents RIGHT lever
% Author:  Zhouxiao Lu
% Date: Mar. 27, 2023
% Last modified on: Apr. 24, 2023

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
tic;

% Set the bin size in seconds
bin_size = 0.05;
Data.bin_size = bin_size;

matches = regexp(file, '(\d+)_(\d+)\.mat', 'tokens');
AnimalID = str2double(matches{1}{1});
SessionID = str2double(matches{1}{2});
%%
trial_data = Data.trial_neuron_timestamps;

% Define the time resolution for the raster plot
dt = 1/Data.Frequency; % seconds
t_extra = Data.ExtraMarginTime;


%% Loop over each trial
save_name = strcat(DataFolder,'/',num2str(AnimalID),'/2_RasterPlots/',num2str(AnimalID),'_',num2str(SessionID));
if ~exist(save_name,'dir')
    mkdir(save_name);
end

for i = 1:size(trial_data, 2)
    figure('Visible','off');
    % Get the timestamps for the current trial
    trial_spikes = trial_data(:, i);
    sample_type = Data.SamplePosition{i,1};
    response_type = Data.ResponsePosition{i,1};

    % Define the time window for plotting
    t_start = Data.trials_timestamps(i,1)-t_extra;
    t_end = Data.trials_timestamps(i,2)+t_extra;

    % Create time bins
    bins = t_start:bin_size:t_end;
    num_bins = length(bins)-1;
    t_window = [bins(1),bins(end)];

    % Allocate memory for the binned spike counts
    binned_spikes = zeros(1,num_bins);

    % Loop over each neuron
    for j = 1:size(trial_data, 1)
        % Get the timestamps for the current neuron in the current trial
        neuron_spikes = trial_spikes{j,1}';
        % Loop through each time bin
        for k = 1:num_bins
            % Identify spikes that occurred within the current bin
            spikes_in_bin = neuron_spikes >= bins(k) & neuron_spikes < bins(k+1);
            
            % Count the number of spikes in the bin
            binned_spikes(k) = sum(spikes_in_bin);
        end
        Data.binned{j,i} = binned_spikes;
        %neuron_spikes = binned_spikes;
        % Plot a vertical line for each spike
        subplot(2,1,1)
        spikes = t_start + find(binned_spikes)*bin_size;
        y_vals = j * ones(size(spikes));
        line([spikes; spikes], [y_vals-0.5; y_vals+0.5], 'Color', 'k');
        xx = (bins(1:end-1) + bins(2:end)) / 2;
        subplot(2,1,2)
        plot(xx,binned_spikes+j);
        hold on
    end
    hold off

    for m = 1:2
        subplot(2,1,m)
        % Set the axis limits and labels
        xlim(t_window);
        ylim([0, size(trial_data, 1)+1]);
        ylabel('Neuron');
        xlabel('Time (s)');
        % Draw reactangles  
        % draw the SAMPLE rectangle
        t_peri = Data.PeriTime;
        rect_x = t_start + t_extra;
        rect_y = 0;
        rect_width = 2*t_peri;
        rect_height = size(trial_data,1);
        if strcmp(sample_type,'LEFT')
            rectangle('Position', [rect_x, rect_y, rect_width, rect_height], ...
                      'FaceColor', [0, 0, 1, 0.5],'EdgeColor', 'none');
        else
            rectangle('Position', [rect_x, rect_y, rect_width, rect_height], ...
                      'FaceColor', [1, 0, 0, 0.5],'EdgeColor', 'none');
        end
    
        % draw the RESPONSE rectangle
        rect_x = t_end - t_extra - 2*t_peri;
        rect_y = 0;
        rect_width = 2*t_peri;
        rect_height = size(trial_data,1);
        if strcmp(response_type,'LEFT')
            rectangle('Position', [rect_x, rect_y, rect_width, rect_height], ...
                      'FaceColor', [0, 0, 1, 0.5],'EdgeColor', 'none');
        else
            rectangle('Position', [rect_x, rect_y, rect_width, rect_height], ...
                      'FaceColor', [1, 0, 0, 0.5],'EdgeColor', 'none');
        end
    end
    
    
    % Set the title of the subplot to the trial number
    sgtitle(sprintf('Animal %d, Session %d, Trial %d, Bin size %d', AnimalID, SessionID, i,bin_size));
    fig_name = sprintf('trial%d.png', i);
    file_path = fullfile(save_name, fig_name);

    % Save the figure in PNG/EPS format with high resolution
    saveas(gcf,file_path)
    % saveas(gcf,file_path,'epsc')
end

mat_name = fullfile(path, file);
save(mat_name,'Data')

elapsed_time = toc;
fprintf('Elapsed time: %.2f seconds\n', elapsed_time);