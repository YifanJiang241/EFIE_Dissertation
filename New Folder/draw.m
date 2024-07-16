% Load the data from the file
data = load('Jf_fast_2.dat');
data2 = load('Ef_fast_2.dat')

% Separate the data into x and y components
x = data(:, 1);
y = data(:, 2);

x2 = data2(:,1);
y2 = data2(:,2);

% Create a figure
figure;

% Plot the data with a thinner line
plot(x, y, 'b-o','MarkerSize',2, 'LineWidth', 0.15); % Setting the LineWidth to 0.25 makes the line thinner

hold on

%plot(x2, y2, 'r-o','MarkerSize',2, 'LineWidth', 0.15);
% Add labels and title
xlabel('X axis');
ylabel('Y axis');
title('J ');

% Add grid for better readability
grid on;
