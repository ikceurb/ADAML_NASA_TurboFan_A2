clear; close all; clc
clear all

% Load the data
load trainFD001.mat

% Extract the data into a matrix
table_data = trainFD001{:,:};

% Define index names
index_names = {'engine', 'cycle'};

% Define setting names
setting_names = {'setting_1', 'setting_2', 'setting_3'};

% sensor names
sensor_names = cell(1, 21);  % Preallocate a cell array for 21 sensors
for i = 1:21
    sensor_names{i} = ['sensor ', num2str(i)];
end


% Concatenate all the names into one cell array
col_names = [index_names, setting_names, sensor_names];
disp(col_names);



% Convert the table_data back to a table and assign column names
dftrain = array2table(table_data, 'VariableNames', col_names);

% Separate the data to different matrices for clearance
%operational_settings = table_data(:,3:5); % Operational settings
time = dftrain{:, 2}; % Time series data Integers
data = dftrain{:, 6:end}; % Sensor data
iteration = dftrain{:, 1}; % group/iteration (engine number)

format shortG


num_sensors = size(data, 2);
fprintf('Total sensors: %d\n', num_sensors);

% Calculating summary statistics for each sensor
mean_vals = mean(data);
std_vals = std(data);
min_vals = min(data);
max_vals = max(data);
unique_vals = zeros(1, num_sensors);

for i = 1:num_sensors
    unique_vals(i) = length(unique(data(:, i)));
end

fprintf('Basic Statistics for each sensor:\n');
disp(table(sensor_names', mean_vals', std_vals', min_vals', max_vals', unique_vals', ...
    'VariableNames', {'Sensor No.', 'Mean', 'STD', 'Min', 'Max', 'Unique'}));



boxplot(data,'Labels',sensor_names)
figure;

%Histogram for each sensor

for i = 1:num_sensors
    subplot(3, 7, i);
    histogram(data(:, i));
    title(sprintf('Sensor %d', i));
    xlabel('Value');
    ylabel('Count');
end

% Calculate maximum time cycles for each unit number
units = unique(dftrain.engine);
max_time_cycles = zeros(length(units), 1);
for i = 1:length(units)
    max_time_cycles(i) = max(time(iteration == units(i)));
end

% Create figure and adjust the size for a long narrow plot
figure('Position', [100, 100, 800, 1200]); % Width 800, Height 1200 (tall figure)

% Plot horizontal bar graph
barh(units, max_time_cycles, 'FaceColor', [0.2 0.6 0.8]);

% Customize the plot
ylim([0, 100]); % Set y-axis limits (0 to 100)
xlabel('Number of Cycles', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Engine Number', 'FontWeight', 'bold', 'FontSize', 12);
title('Turbofan Engine Maximum Cycles', 'FontWeight', 'bold', 'FontSize', 14);


% Adjust yticks and ylim for better spacing (1, 5, 10, ..., 100 intervals)
ytick_interval = [1, 5:5:100];  % 1, 5, 10, ..., 100
yticks(ytick_interval);
ylim([0 length(units) + 1]);

% Grid and layout adjustments
grid on;
set(gca, 'YDir', 'reverse'); % Optional: to have Engine 1 at the top
box on;

% Plot histogram with count on y-axis
figure;

% Histogram with 20 bins (default normalization is 'Count')
histogram(max_time_cycles, 20, 'FaceColor', [0.2 0.6 0.8]);
hold on;

% Kernel Density Estimate (KDE) (optional)
[f, xi] = ksdensity(max_time_cycles);
yyaxis right; % Create second y-axis for KDE curve
plot(xi, f, 'r', 'LineWidth', 2);
ylabel('Probability Density', 'FontWeight', 'bold', 'FontSize', 15);

% Labels and title
xlabel('Max Time Cycle', 'FontWeight', 'bold', 'FontSize', 15);
yyaxis left; % Switch back to the left y-axis
ylabel('Count', 'FontWeight', 'bold', 'FontSize', 15);
title('Distribution of Maximum Time Cycles', 'FontWeight', 'bold', 'FontSize', 20);

% Grid and hold off
grid on;
hold off;


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



%Function to add RUL column
function df = add_RUL_column(df)
    % Group by 'unit_number' and calculate the maximum 'time_cycles' for each unit
    [unitss, ~, unit_indices] = unique(df.engine); 
    max_time_cycles = accumarray(unit_indices, df.cycle, [], @max);
    
    % Add a 'max_time_cycle' column to the table
    max_time_cycles_per_unit = max_time_cycles(unit_indices);
    df.max_time_cycle = max_time_cycles_per_unit;
    
    % Calculate the Remaining Useful Life (RUL)
    df.RUL = df.max_time_cycle - df.cycle;
    
    % Remove the 'max_time_cycle' column
    df.max_time_cycle = [];
end

% Call the function to add the RUL column
dftrain = add_RUL_column(dftrain);

% Display a preview of the unit_number and RUL columns
disp(dftrain(1:10, :));



% Remove columns with constant sensor values
sens_const_values = {};
sens_const_values_other = sensor_names; 

for i = 1:length(sensor_names)
    feature = sensor_names{i}; 
    try
        if min(dftrain.(feature)) == max(dftrain.(feature))
            sens_const_values{end+1} = feature; % Add the sensor name
            sens_const_values_other(strcmp(sens_const_values_other, feature)) = []; 
        end
    catch
        % Do nothing if an error occurs
    end
end

disp(sens_const_values); % Display sensors with constant values
dftrain(:, sens_const_values) = []; % Remove them from dftrain
%dftrain(:, {"cycle", 'engine','setting_1', 'setting_2', 'setting_3'} ) = [];


% Randomly shuffle the units
rng('shuffle'); % For reproducibility
shuffled_units = units(randperm(length(units)));

% Calculate split index
split_idx = round(0.8 * length(shuffled_units));

% Split units into training and testing
train_units = shuffled_units(1:split_idx);
test_units = shuffled_units(split_idx + 1:end);

% Create training and testing datasets
train_data = dftrain(ismember(dftrain.engine, train_units), :);
test_data = dftrain(ismember(dftrain.engine, test_units), :);

% Count missing values in each column
missing_dftrain = dftrain
missing_dftrain(:, {'cycle', 'engine', 'setting_1', 'setting_2', 'setting_3'}) = [];

missing_counts = sum(ismissing(missing_dftrain), 1); % Sum missing values per column

% Create a table with the desired row and column names
missing_table = table(missing_counts', 'VariableNames', {'MissingCount'}, ...
    'RowNames', missing_dftrain.Properties.VariableNames');

% Display the table
disp(missing_table);


% Example selecting specific columns 
std_train = train_data{:, 6:(width(dftrain)-1)};  
[std_train, mu, sigma] = zscore(std_train);  % Apply zscore

% Standardizing the test data (using mu and sigma from the training data)
std_test = (test_data{:, 6:(width(test_data)-1)} - mu) ./ sigma;

figure;
boxplot(std_train,'Labels',sens_const_values_other)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
