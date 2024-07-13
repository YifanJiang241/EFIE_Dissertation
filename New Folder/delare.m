% Constants
EXP = exp(1);
PI = 3.14159265358979323846;
CONJ = -1; % You may not need this in MATLAB since conjugate operations are built-in
NO_CONJ = 1; % Same as above
TRANS = 1; % If this is used for matrix transpose, MATLAB uses ' for conjugate transpose and .' for non-conjugate transpose

% No need to define complex type, MATLAB has built-in support for complex numbers

% Global Variables
Eamp = 1.0;
Epsilon_0 = 8.854e-12;
Epsilon_d = 2.47912e-11;
Mu_0 = 12.56637061e-7;
c = 1.0 / sqrt(Mu_0 * Epsilon_0);
GrossStep = 10.0;
f = 970e6;
Lambda = c / f;
DeltaX = Lambda / 4.0;
Omega = 2.0 * PI * f;

j = 1j; % MATLAB's imaginary unit is 1i or 1j by default
GrossNoSteps = 70;
Beta_0 = Omega * sqrt(Mu_0 * Epsilon_0);
Eta_0 = sqrt(Mu_0 / Epsilon_0);
TOL = 10e-15;
NoLinesubs = floor((GrossStep * GrossNoSteps) / DeltaX);
Xsource = 0.0;
Ysource = 442.0;
I = 1.0;

% In MATLAB, arrays are typically 1-indexed, so adjust if necessary
X = zeros(1, 385);
Y = zeros(1, 385);




% MATLAB handles file I/O differently, so we don't declare file pointers.
% Instead, we will use fopen, fprintf/fscanf, and fclose directly in the file operations.

% Integer variables
i = 0; p = 0; q = 0; n = 0; index = 1;

% Calculate NoLinesubs
NoLinesubs = floor((GrossStep * GrossNoSteps) / DeltaX);

% In MATLAB, you don't need to specify 'new' for allocation. 
% Also, MATLAB arrays are 1-indexed and dynamically typed.

% Double arrays
ModJ = zeros(1, NoLinesubs);
ModEt = zeros(1, NoLinesubs);

% Complex arrays
J = complex(zeros(1, NoLinesubs));
Et = complex(zeros(1, NoLinesubs));
Sigma = complex(zeros(1, NoLinesubs));

% Single complex variable
SUM = complex(0, 0);















function main
    % Load Terrain Profile Data from X.04
    fileID = fopen('X.04', 'r');
    % 逐行读取
    for i = 1:385
        data = fscanf(fileID, '%f %f', [1 2]);  % 每次读取一行，每行有两个浮点数
        if ~isempty(data)  % 确保读取到了数据
            X(i) = data(1);
            Y(i) = data(2);
        else
            break;  % 如果没有读取到数据，则跳出循环
        end
    end
    fclose(fileID);
    X,Y
    % Forward scattering calculation for J
    test0 = R_p_q(1,2)
    test1 = R_source_p(1)
    test2 = EiRad(test1, 1)
    test3 = Zself(1)
    test4 = x_func(1)
    test5 = y_func(1)
    disp(R_source_p(1));
    disp(EiRad(R_source_p(1), 1));
    disp(Zself(1));
    
    test = EiRad(R_source_p(1), 1) / Zself(1); 
    J(1) = test;
    for p = 1:NoLinesubs
        SUM = complex(0.0, 0.0);
        for q = 1:(p-1)
            SUM = SUM + R_p_q(q, q+1) * Z(p, q) * J(q);
        end
        J(p) = (EiRad(R_source_p(p), p) - SUM) / Zself(p);  % Equation (6)
    end
    
    Backscattering calculation for J
    for p = NoLinesubs:-1:2
        SUM = complex(0.0, 0.0);    
        for q = NoLinesubs:-1:(p+1)
            SUM = SUM + R_p_q(q, q+1) * Z(p, q) * J(q);
        end
        J(p) = J(p) + (-1.0 * SUM) / Zself(p);  % Equation (6)
    end
    
    % Write J to file J.dat
    fileID = fopen('J.dat', 'w');
    for index = 1:NoLinesubs
        fprintf(fileID, '%f  %e\n', x_func(index), abs(J(index)));
    end
    fclose(fileID);
    
    % ... Rest of your main function code ...

end



















function E = EiRad(dist, p)
    global Beta_0 Omega Epsilon_0;
    E = -((Beta_0^2) / (4.0 * Omega * Epsilon_0)) * H02(Beta_0 * dist);
end











function H = H02(Arg)
    % You need to define the complex Hankel function H02
    % MATLAB has bessel functions, e.g., besselj and bessely for first and second kind
    H = complex(besselj(0, Arg), -bessely(0, Arg));
end




function H = Zself(i)
    global Beta_0 Omega Epsilon_0 EXP;
    Linesubln = R_p_q(i, i + 1);
    H = complex(((Beta_0^2) / (4.0 * Omega * Epsilon_0)) * Linesubln, ...
                -((Beta_0^2) / (4.0 * Omega * Epsilon_0)) * ((2.0 * Linesubln) / pi) * log((1.781 * Beta_0 * Linesubln) / (4.0 * EXP)));
end








function R = R_p_q(p, q)
    R = sqrt((x_func(q) - x_func(p))^2 + (y_func(q) - y_func(p))^2);
end


function H = Z(p, q)
    global Beta_0 Omega Epsilon_0;
    H = ((Beta_0^2) / (4.0 * Omega * Epsilon_0)) * H02(Beta_0 * R_p_q(p, q));
end



% function result = x(a)
%     global DeltaX;
%     result = a * DeltaX;
% end

% function s = y(a)
%     global DeltaX GrossStep Y;
%     Temp = a * DeltaX / GrossStep;
%     Index = floor(Temp);
%     Prop = Temp - Index;
%     s = Y(Index + 1) + (Prop * (Y(Index + 2) - Y(Index + 1)));
% end

function R = R_source_p(p)
    % 根据提供的信息，计算源点和第p个点之间的距离
    % 假设 x_func 和 y_func 是计算坐标的函数
    global Xsource  Ysource;
    R = sqrt((Xsource - x_func(p))^2 + (Ysource - y_func(p))^2);
end


function result = x_func(a)
    global DeltaX;
    result = a  * DeltaX;
end

function result = y_func(a)
    global DeltaX GrossStep Y;
    % 补全 y_func 函数，a 为传入的索引，DeltaX 为步长，GrossStep 为总步数，Y 为数组
    Temp = a  * DeltaX / GrossStep;  % 计算 Temp
    Index = floor(Temp);  % 取整数部分
    Prop = Temp - double(Index);  % 计算 Prop
    % if Index == 0 % 防止索引为 0
    %     Index = 1;
    % end
    % if Index >= length(Y)-1 % 防止索引超过数组长度
    %     Index = length(Y) - 2;
    % end
    result = Y(Index) + Prop * (Y(Index + 1) - Y(Index));
end