% 加载数据，指定无标题行且使用空白符分隔
opts = detectImportOptions(data_path, 'Delimiter', 'whitespace', 'NumHeaderLines', 0);
signal_data = readtable(data_path, opts);

% 转换距离单位从米到公里
signal_data.Var1 = signal_data.Var1 / 1000; % 假设第一列是距离

% 重命名列以方便引用
signal_data.Properties.VariableNames = {'Distance', 'Signal_Strength'};

% 定义参数
f = 970;  % 频率 (MHz)
hte = 442;  % 发射天线高度 (m)
hre = 1.5;  % 接收天线高度 (m)
d = linspace(0, 0.7, 100);  % 距离 (km), 从0开始

% 计算a(hre)
a = (1.1 * log10(f) - 0.7) * hre - (1.56 * log10(f) - 0.8);

% 基本的城市模型损失
L_urban = 69.55 + 26.16 * log10(f) - 13.82 * log10(hte) - a + (44.9 - 6.55 * log10(hte)) * log10(d + 0.1);

% 郊区和开放区域模型损失
L_suburban = L_urban - 2 * (log10(f/28).^2) - 5.4;
L_open_area = L_urban - 4.78 * (log10(f).^2) + 18.33 * log10(f) - 40.94;

% 翻转路径损失曲线并调整至信号强度起始高度
offset_urban = signal_data.Signal_Strength(1) + L_urban(1);
offset_suburban = signal_data.Signal_Strength(1) + L_suburban(1);
offset_open_area = signal_data.Signal_Strength(1) + L_open_area(1);

L_urban_flipped = -L_urban + offset_urban;
L_suburban_flipped = -L_suburban + offset_suburban;
L_open_area_flipped = -L_open_area + offset_open_area;

% Urban model plot
figure;
plot(d, L_urban_flipped, 'b-', 'DisplayName', 'Urban Model (Adjusted)');
hold on;
scatter(signal_data.Distance, signal_data.Signal_Strength, 'red', 'DisplayName', 'Signal Strength Data');
title('Urban Model Path Loss and Signal Strength Comparison');
xlabel('Distance (km)');
ylabel('Signal Strength/Path Loss (dB)');
legend('show');
grid on;
saveas(gcf, 'F:/tcd/s1/dissertation/pathlossmodel/Urban_Model.png');
close;

% Suburban model plot
figure;
plot(d, L_suburban_flipped, 'k-', 'DisplayName', 'Suburban Model (Adjusted)');
hold on;
scatter(signal_data.Distance, signal_data.Signal_Strength, 'red', 'DisplayName', 'Signal Strength Data');
title('Suburban Model Path Loss and Signal Strength Comparison');
xlabel('Distance (km)');
ylabel('Signal Strength/Path Loss (dB)');
legend('show');
grid on;
saveas(gcf, 'F:/tcd/s1/dissertation/pathlossmodel/Suburban_Model.png');
close;

% Open area model plot
figure;
plot(d, L_open_area_flipped, 'g-', 'DisplayName', 'Open Area Model (Adjusted)');
hold on;
scatter(signal_data.Distance, signal_data.Signal_Strength, 'red', 'DisplayName', 'Signal Strength Data');
title('Open Area Model Path Loss and Signal Strength Comparison');
xlabel('Distance (km)');
ylabel('Signal Strength/Path Loss (dB)');
legend('show');
grid on;
saveas(gcf, 'F:/tcd/s1/dissertation/pathlossmodel/Open_Area_Model.png');
close;
