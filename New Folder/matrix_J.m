% 初始化全局变量
global Eamp Epsilon_0 Epsilon_d Mu_0 c GrossStep f Lambda DeltaX Omega;
global Beta_0 Eta_0 TOL NoLinesubs Xsource Ysource I X Y;

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

    %% 


% 调用precomputeInfluenceMatrix函数来计算影响矩阵
InfluenceMatrix = precomputeInfluenceMatrix(NoLinesubs);



    %% 

% 假设EiRadVec是所有点的EiRad值的向量
EiRadVec = arrayfun(@(p) EiRad(R_source_p(p), p), 1:NoLinesubs);

% 假设ZselfVec是所有点的Zself值的向量
ZselfVec = arrayfun(@(p) Zself(p), 1:NoLinesubs);

% 创建方程组的左侧矩阵A，包含自身阻抗和影响矩阵的组合
A = diag(ZselfVec) - [zeros(1, NoLinesubs); tril(InfluenceMatrix(2:end, :), -1)];

% 创建方程组的右侧向量B，包含所有点的外部影响
B = EiRadVec.';

% 解线性方程组A * J = B以找到J
J = B \ A;

    

    %% 

fileID = fopen('Jf_m.dat', 'w');

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

function InfluenceMatrix = precomputeInfluenceMatrix(NoLinesubs)
    global Beta_0 Omega Epsilon_0;
    
    % 初始化影响矩阵为复数零矩阵
    InfluenceMatrix = complex(zeros(NoLinesubs, NoLinesubs));
    
    % 预计算所有点的x和y坐标
    allX = x(1:NoLinesubs);
    allY = y(1:NoLinesubs);
    
    for p = 1:NoLinesubs
        p
        for q = 1:p-1

            % 计算p和q之间的距离
            R = sqrt((allX(p) - allX(q))^2 + (allY(p) - allY(q))^2);
            
            % 计算Hankel函数的值并更新影响矩阵
            H = H02(Beta_0 * R);
            InfluenceMatrix(p, q) = ((Beta_0^2) / (4.0 * Omega * Epsilon_0)) * H;
        end
    end
end




function result = x(a)
    global DeltaX;
    result = a * DeltaX; % 这里已经可以处理向量输入
end



function s = y(a)
    global Y DeltaX GrossStep;
    % 计算临时变量Temp，表示a对应的相对位置
    Temp = (a * DeltaX) / GrossStep;
    % MATLAB的索引从1开始，所以需要对计算得到的索引进行调整
    Index = floor(Temp) + 1;
    % 计算比例Prop
    Prop = Temp - floor(Temp);
    
    % 确保Index数组中的最后一个元素不超过Y的长度
    Index(Index >= length(Y)) = length(Y) - 1;
    
    % 使用线性插值计算s，适用于向量
    s = Y(Index) + Prop .* (Y(Index + 1) - Y(Index));
end



function E = EiRad(dist, p)
    global Beta_0 Omega Epsilon_0;
    % 针对向量dist进行操作，利用MATLAB的向量化特性
    E = -((Beta_0^2) / (4.0 * Omega * Epsilon_0)) .* H02(Beta_0 .* dist);
end




function H = Zself(i)
    global Beta_0 Omega Epsilon_0 PI;
    % 计算R_p_q向量，这里我们需要一个方法来计算连续元素间的距离
    Linesubln = R_p_q(i, i + 1); % 修改为向量化计算

    % 对于向量Linesubln，进行向量化计算
    realPart = ((Beta_0^2) / (4.0 * Omega * Epsilon_0)) .* Linesubln;
    imagPart = -((Beta_0^2) / (4.0 * Omega * Epsilon_0)) .* ((2.0 .* Linesubln) / pi) .* log((1.781 .* Beta_0 .* Linesubln) / (4.0 * exp(1)));

    H = complex(realPart, imagPart);
end

function R = R_source_p(p)
    global Xsource Ysource X Y;
    % 计算所有p点的x和y坐标
    xp = x(p);
    yp = y(p);
    % 向量化计算源点到每个点的距离
    R = sqrt((Xsource - xp).^2 + (Ysource - yp).^2);
end


function R = R_p_q(p, q)
    % 计算p和q点的x和y坐标
    xp = x(p);
    yp = y(p);
    xq = x(q);
    yq = y(q);
    % 向量化计算p和q之间的距离
    R = sqrt((xq - xp).^2 + (yq - yp).^2);
end

function H = H02(Arg)
    % 直接使用MATLAB内置的向量化贝塞尔函数
    H = besselj(0, Arg) - 1i * bessely(0, Arg);
end

