% Load the data from the file
data = load('Jf_fast_2.dat');
data2 = load('J.dat');

% Separate the data into x and y components
x = data(:, 1);
y = data(:, 2);

x2 = data2(:, 1);
y2 = data2(:, 2);

% Filter out zero values
non_zero_indices = y ~= 0;
x_non_zero = x(non_zero_indices);
y_non_zero = y(non_zero_indices);

non_zero_indices2 = y2 ~= 0;
x2_non_zero = x2(non_zero_indices2);
y2_non_zero = y2(non_zero_indices2);

% Create a figure
figure;

% Plot the data with a thinner line
plot(x_non_zero, y_non_zero, 'b-o', 'MarkerSize', 2, 'LineWidth', 0.15); % Setting the LineWidth to 0.15 makes the line thinner

hold on

plot(x2_non_zero, y2_non_zero, 'r-o', 'MarkerSize', 2, 'LineWidth', 0.15);

% Add labels and title
xlabel('X axis');
ylabel('Y axis');
title('J ');

% Add grid for better readability
grid on;
