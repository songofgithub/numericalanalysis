function [Tmk_array] = RombergIntegral(fn_str,integral_lower_limit,integral_upper_limit,precision)
%ROMBERGINTEGRAL 使用龙贝格算法求函数积分
%   fn_str：函数表达式
%   integral_lower_limit: 积分下限
%   integral_upper_limit：积分上限
%   precision：误差要求
%   Tmk_array：迭代的所有值，共四列，第一列为梯形迭代次数，第二列为加速次数，第三列为积分值，第四列为精度
%   Anhui University S.L. Xia 2024年1月11日

    f = str2func(['@(x)' vectorize(fn_str)]); % 获取函数句柄
    % 数据初始化
    a = integral_lower_limit;
    b = integral_upper_limit;
    max_iterators = 1000; % 最大初始化数值行数
    Tmk_array = zeros(max_iterators,4); % 初始化迭代积分值Tmk，k=0,1,2...，存储迭代次数，加速次数，积分结果，精度
    k = 0; % 初始化迭代次数为0    

    % 第一步：计算初值T1        
    T2k = ((b - a) / 2) * (f(a) + f(b)); % k=0    
    Tmk_array(1,1) = k;
    Tmk_array(1,2) = 0;
    Tmk_array(1,3) = T2k;

    temp_precision = inf; % 将初始精度设为无穷大
    % 第四步：判断精度是否符合要求
    while(temp_precision >= precision)
        k = k + 1;
        %if(k > 23)
        %    throw(MException('MATLAB:TooManyInputs',['已迭代' num2str(k) '次，过多迭代次数会导致计算十分耗时']));
        %end
        % 第二步：计算梯形迭代的结果
        sum = 0;
        for i = 1 : 2 ^ (k - 1)            
            sum = sum + f(a + (2 * i - 1) * ((b - a) / (2 ^ k))); % 迭代公式中的求和部分        
        end
        T2k = 1 / 2 * T2k + (b - a) / (2 ^ k) * sum; % 迭代计算出新的T2^k             
        T2k_index = ((1 + k) * k) / 2 + 1; % 每行元素按k元素等差数列递增，每次新计算的梯形迭代值应存放的行索引
        T2k_index_pre = ((1 + (k - 1)) * (k - 1)) / 2 + 1; % 上一个梯形迭代值存放的行索引        
        Tmk_array(T2k_index,1) = k;
        Tmk_array(T2k_index,2) = 0;
        Tmk_array(T2k_index,3) = T2k;        

        % 第三步：求加速值
        for j = 1 : k
            Tmk_array(T2k_index + j,1) = k - j;
            Tmk_array(T2k_index + j,2) = j;            
            Tmk_array(T2k_index + j,3) = ((4 ^ j) * Tmk_array(T2k_index + j - 1,3) - Tmk_array(T2k_index_pre + j - 1,3)) / ((4 ^ j) - 1);                        
        end          

        % 计算精度
        temp_precision = abs(Tmk_array(T2k_index + k,3) - Tmk_array(T2k_index_pre + k - 1,3));        
        Tmk_array(T2k_index + k,4) = temp_precision;
    end    

    rows = ((1 + (k + 1)) * (k + 1)) / 2;
    Tmk_array = Tmk_array(1:rows,:); % 去除尾部无用的零
end