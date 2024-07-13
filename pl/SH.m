% 参数定义
f = 970;  % 频率 (MHz)
hte = 442;  % 发射天线高度 (m)
hre = 1.5;  % 接收天线高度 (m)
d = linspace(0.1, 0.7, 100);  % 距离 (km), 避免0公里开始

% Okumura-Hata模型常数
a = (1.1*log10(f) - 0.7)*hre - (1.56*log10(f) - 0.8);
c = 0;
if f < 150
    c = 8.29 * (log10(1.54*hre))^2 - 1.1;
elseif f > 1500
    c = 3.2 * (log10(11.75*hre))^2 - 4.97;
end

% 城市环境损耗公式
L = 69.55 + 26.16*log10(f) - 13.82*log10(hte) - a + (44.9 - 6.55*log10(hte)).*log10(d) + c;

% 绘制图像
figure;
plot(d, L);
title('路径损失 (dB) vs 距离 (km)');
xlabel('距离 (km)');
ylabel('路径损失 (dB)');
grid on;
