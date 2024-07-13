function write(J)
% 打开文件，准备写入，获取文件句柄
fileID = fopen('Jf.dat', 'w');

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
end