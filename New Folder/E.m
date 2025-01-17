% 初始化全局变量
global Eamp Epsilon_0 Epsilon_d Mu_0 c GrossStep f Lambda DeltaX Omega;
global Beta_0 Eta_0 TOL NoLinesubs Xsource Ysource I X Y;
%高度442 = 52 + 390
% 常量定义
Eamp = 1.0;
Epsilon_0 = 8.854e-12;
Epsilon_d = 2.47912e-11;
Mu_0 = 12.566370614e-7;
c = 1.0 / sqrt(Mu_0 * Epsilon_0); % 光速
GrossStep = 10.0;
f = 970e6; % 频率
Lambda = c / f; % 波长
DeltaX = Lambda / 4.0; % 空间步长
Omega = 2.0 * pi * f; % 角频率
Beta_0 = Omega * sqrt(Mu_0 * Epsilon_0); % 波数
Eta_0 = sqrt(Mu_0 / Epsilon_0); % 本征阻抗
TOL = 10e-15; % 容忍误差
NoLinesubs = floor((GrossStep * 70) / DeltaX); % 根据实际情况调整
Xsource = 0.0; % 源位置X坐标
Ysource = 442.0; % 源位置Y坐标
I = 1.0; % 源电流

cst = 4.8
stp = 5
% cst = 1
% stp = 1

% 地形数据初始化
X = zeros(1, 385); % 假定有385个数据点
Y = zeros(1, 385);

% 请根据实际情况调整以上变量的值






    
    % 初始化J数组
    J = complex(zeros(1, NoLinesubs)); % 预分配复数数组
    
    % 读取地形
    % 逐行读取
     % 打开文件
    data = load("X.04");
    X = data(:, 1)';
    Y = data(:, 2)';
    
    % Forward Scattering
    %% 
    
    J(1) = EiRad(R_source_p(1), 1) / Zself(1);
    for p = 1:stp:NoLinesubs %numofgroups = nolinesubs / stp
        SUM = 0 + 0i; 
        for q = 1:stp:p-1 
           
            SUM = SUM + R_p_q(q, q+1) * Z(p, q) * J(q) * cst;
        end
        p
        J(p) = (EiRad(R_source_p(p), p) - SUM) / Zself(p); % Equation (6)
    end
    
    % 


    %% backscattering

    SUM = 0 + 0i; 
    for p = NoLinesubs :stp: 1
        SUM = 0 + 0i;
        for q = NoLinesubs :stp: p+1
            SUM = SUM + R_p_q(q, q+1) * Z(p, q) * J(q) * cst;
        end
        J(p) = J(p) -SUM/Zself(p);
    end


    %% 


    % 打开文件，准备写入，获取文件句柄
fileID = fopen('J_fast_stp5_matlab_fb.dat', 'w');

% 检查文件是否成功打开
if fileID == -1
    error('Failed to open file.');
end

% 遍历J数组，写入x坐标和J的模长(abs(J))到文件
for index = 1:NoLinesubs % MATLAB索引从1开始
    % 注意x函数需要根据您的实现来调用
    % 假设x坐标可以直接用x(index)获取，这里可能需要调整
    fprintf(fileID, '%f %f\n', x(index), abs(J(index)));
end

% 关闭文件
fclose(fileID);
%% 


% 打开文件，准备写入
fileID = fopen('E_fast_stp5_matlab_fb.dat', 'w');

% 检查文件是否成功打开
if fileID == -1
    error('Failed to open file.');
end

Et = complex(zeros(1, NoLinesubs)); 

for index = 1:stp:NoLinesubs
    Et(index) = complex(0.0, 0.0); 
    index
    for n = 1:stp:index 
        Et(index) = Et(index) + (J(n) * R_p_q(n, n+1) * Zcomplex(R_surf_obs(n, index))) * cst;
    end

    Et(index) = EiRad(R_source_obs(index), index) - Et(index); 

    % 写入x坐标和Et的20*log10(abs(Et)/sqrt(R_source_obs(index)))到文件
    fprintf(fileID, '%f %f\n', x(index), 20.0 * log10(abs(Et(index)) / sqrt(R_source_obs(index))));
end

% 关闭文件
fclose(fileID);
%% 








% Function Definitions





function R = R_source_p(p)
    global Xsource Ysource X Y DeltaX;
    R = sqrt((Xsource - x(p))^2 + (Ysource - y(p))^2);
end



function result = x(a)
    global DeltaX;
    % 假设DeltaX已在全局变量中定义
    result = a * DeltaX; % 直接计算x坐标
end

function R = R_p_q(p, q)
    % 使用之前定义的x和y函数来获取p和q点的坐标
    xp = x(p);
    yp = y(p);
    xq = x(q);
    yq = y(q);
    
    % 计算两点间的欧几里得距离
    R = sqrt((xq - xp)^2 + (yq - yp)^2);
end


function s = y(a)
    global Y DeltaX GrossStep;
    % 计算临时变量Temp，表示a对应的相对位置
    Temp = (a * DeltaX) / GrossStep;
    % MATLAB的索引从1开始，所以需要对计算得到的索引进行调整
    Index = floor(Temp) + 1;
    % 计算比例Prop
    Prop = Temp - floor(Temp);
    
    % 当Index为数组的最后一个元素时，需要特别处理以避免索引越界
  
    
    % 使用线性插值计算s
    s = Y(Index) + Prop * (Y(Index + 1) - Y(Index));
end


function E = EiRad(dist, p)
    global Beta_0 Omega Epsilon_0;
    % 使用MATLAB的复数运算计算E
    E = -((Beta_0^2) / (4.0 * Omega * Epsilon_0)) * H02(Beta_0 * dist);
end

function H = H02(Arg)
    % 使用MATLAB内置的贝塞尔函数j0和y0计算Hankel函数的值
    H = besselj(0, Arg) - 1i * bessely(0, Arg);
end



function H = Zself(i)
    global Beta_0 Omega Epsilon_0 PI;
    % 确保所需的全局变量已经定义和初始化

    Linesubln = R_p_q(i, i + 1); % 这应该返回一个标量值

    % 计算复数H的实部和虚部
    realPart = ((Beta_0^2) / (4.0 * Omega * Epsilon_0)) * Linesubln;
    imagPart = -((Beta_0^2) / (4.0 * Omega * Epsilon_0)) * ((2.0 * Linesubln) / pi) * log((1.781 * Beta_0 * Linesubln) / (4.0 * exp(1)));

    H = complex(realPart, imagPart); % 使用complex创建复数
end

function H = Z(p, q)
    global Beta_0 Omega Epsilon_0;
    
    % 计算点p和点q之间的距离
    R = R_p_q(p, q);
    
    % 使用H02函数计算Hankel函数的值
    % 并计算复数阻抗H
    H = ((Beta_0^2) / (4.0 * Omega * Epsilon_0)) * H02(Beta_0 * R);
end


function R = R_source_obs(p)
    global Xsource DeltaX;
    % 计算R_source_obs值
    R = sqrt(((Xsource - x(p))^2) + ((442.0 - y(p) - 2.4)^2));
end


function R = R_surf_obs(p, q)
    % 计算R_surf_obs值
    R = sqrt(((x(q) - x(p))^2) + ((y(q) + 2.4 - y(p))^2));
end


function H = Zcomplex(R)
    global Beta_0 Omega Epsilon_0;

    % 假设H02函数计算与Bessel函数相关的Hankel函数的第二类
    % MATLAB中的besselh函数用于计算Hankel函数，这里假设H02等价于besselh(0, 2, Beta_0*R)
    % 第一个参数0表示计算的是第一类的Hankel函数，2表示第二种类型（Hankel函数的第二类），Beta_0*R是函数的参数

    H = ((Beta_0^2) / (4.0 * Omega * Epsilon_0)) * H02(Beta_0 * R);
end

