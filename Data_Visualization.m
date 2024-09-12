%% Data inspection and evaluation
clear; close all;clc
load trainFD001.mat
%load data/RUL_DATA.mat RULFD001
% Summarize and visualize the data

table_data = trainFD001{:,:};

% Separate the data to different matrices for clearance
%operational_settings = table_data(:,3:5); % Operational settings
time = table_data(:,2); % Time series data Integers
data = table_data(:,6:end); % Sensor data
iteration = table_data(:,1); % group/iterationn
%RUL = RULFD001{:,:}; % Not needed yet

clear table_data

format shortG
clc; close all
num_sensors = size(data, 2);
sensor_names = [1:num_sensors];
fprintf('Total sensors: %d\n', num_sensors);

% Calculating summary statistics for each sensor
mean_vals = mean(data);   
std_vals = std(data);            
min_vals = min(data);            
max_vals = max(data);
unique_vals = zeros(1,num_sensors);
for i = 1:num_sensors
    unique_vals(i) = length(unique(data(:,i)));
end


fprintf('Basic Statistics for each sensor:\n');
disp(table(sensor_names', mean_vals', std_vals', min_vals', max_vals', unique_vals', ...
    'VariableNames', {'Sensor No.','Mean', 'STD', 'Min', 'Max','Unique'}));


%Total sensor data (not separated by iterations)
%Uncomment to see all the data on all runs
%{
figure;
for i = 1:num_sensors
    subplot(7, 3, i);
    plot(data(:, i));
    title(sprintf('Sensor %d', i));
    xlabel('Time');
    ylabel('Value');
end
%}
[normal_data, mu] = zscore(data);

figure;
boxplot(normal_data,'Labels',sensor_names)

%Histogram for each sensor
figure;
for i = 1:num_sensors
    subplot(3, 7, i);
    histogram(data(:, i));
    title(sprintf('Sensor %d', i));
    xlabel('Value');
    ylabel('Count');
end



% Loop over each sensor to compute and plot mean ± std over time
figure;
for sensor_idx = 1:num_sensors
    subplot(7, 3, sensor_idx);
    hold on;
    title(sprintf('Sensor %d: Mean ± 2*STD Over Time', sensor_idx));
    xlabel('Time');
    ylabel(sprintf('Sensor %d Value', sensor_idx));
    
    % Get unique time points
    unique_times = unique(time);
    
    % Preallocate arrays to store mean and std for each time point
    mean_values = zeros(length(unique_times), 1);
    std_values = zeros(length(unique_times), 1);
    
    % Calculate mean and std for each time point across all iterations
    for t = 1:length(unique_times)
        time_indices = time == unique_times(t);  % Find rows corresponding to the current time
        sensor_data_at_time = data(time_indices, sensor_idx);  % Extract sensor data at this time point
        mean_values(t) = mean(sensor_data_at_time);  % Mean value at this time point
        std_values(t) = std(sensor_data_at_time);    % Standard deviation at this time point
     
    end
    
    % Plot the mean line
    plot(unique_times, mean_values, 'b', 'LineWidth', 1.5, 'DisplayName', 'Mean');
    
    % Plot shaded area for mean ± 2 standard deviation
    fill([unique_times; flipud(unique_times)], ...
         [mean_values - 2*std_values; flipud(mean_values + 2*std_values)], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Mean ± 1 STD');
    
    hold off;
end
%% Heatmaps for the correlations between the sensors
% I have taken correlation matrix on every iteration and plotted the mean
% of the correlations and the standard deviation of the iterations to see
% the mean correlation and variation of correlation between variables
% between iterations

% Initialize variables
iterations = unique(iteration);  % Get the unique iterations/groups
num_iterations = length(iterations);  % Number of unique iterations
num_sensors = size(data, 2);  % Number of sensors
all_correlations = zeros(num_sensors, num_sensors, num_iterations);  % Store correlation matrices

% Loop over each iteration/group
for i = 1:num_iterations
    % Extract the data for the current iteration
    iteration_data = data(iteration == iterations(i), :);
    
    % Compute the correlation matrix for this iteration
    corr_matrix = corr(iteration_data);
    
    % Store the correlation matrix
    all_correlations(:, :, i) = corr_matrix;
end

% Calculate the mean correlation matrix across iterations
mean_corr_matrix = round(mean(all_correlations, 3),2);

% Calculate the std correlation matrix across iterations
std_corr_matrix = std(all_correlations, 0, 3);

% Visualize the mean correlation matrix as a heatmap
figure;
heatmap(mean_corr_matrix, 'Colormap', autumn, 'ColorbarVisible', 'on','FontSize',15);
title('Mean Correlation Matrix Across Iterations');

% Visualize the std correlation matrix as a heatmap
figure;
heatmap(std_corr_matrix, 'Colormap', hot, 'ColorbarVisible', 'on','FontSize',15);
title('STD of Correlation Matrix Across Iterations');
