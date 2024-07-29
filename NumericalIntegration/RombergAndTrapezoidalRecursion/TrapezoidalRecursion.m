function[T2k_array] = TrapezoidalRecursion(fn_str,integral_lower_limit,integral_upper_limit,precision)
%TRAPEZOIDALRECURSION 梯形迭代法求函数积分
%   fn_str：函数表达式
%   integral_lower_limit: 积分下限
%   integral_upper_limit：积分上限
%   precision：误差要求，计算时取误差限为3倍的设置值
%   T2k_array：迭代的所有值，共三列，第一列为迭代次数，第二列为积分值，第三列为误差
%   Anhui University S.L. Xia 2024年1月11日

    f = str2func(['@(x)' vectorize(fn_str)]); % 获取函数句柄
    % 数据初始化
    a = integral_lower_limit;
    b = integral_upper_limit;
    max_iterators = 1000; % 最大迭代次数
    T2k_array = zeros(max_iterators,3); % 初始化迭代积分值T2^k，k=0,1,2...，存储迭代次数，迭代结果，迭代误差
    k = 0; % 初始化梯形迭代次数为0
    
    % 第一步：计算初值T1
    tmp = 0; % 记录前一次的结果    
    T2k = ((b - a) / 2) * (f(a) + f(b)); % k=0
    T2k_array(1,1) = k;  
    T2k_array(1,2) = T2k;
    T2k_array(1,3) = abs(T2k - tmp);
    tmp = T2k; % 更新tmp的值为当次迭代结果
    k = k + 1; % 迭代次数加1

    % 第二步：计算新的梯形值T2^k
    while(k < max_iterators + 1 && T2k_array(k,3) >= 3 * precision)  % 此处判断是否超出迭代次数和第三步：判断是否达到精度要求
        sum = 0;
        for i = 1 : 2 ^ (k - 1)            
            sum = sum + f(a + (2 * i - 1) * ((b - a) / (2 ^ k))); % 迭代公式中的求和部分        
        end
        
        T2k = 1 / 2 * tmp + (b - a) / (2 ^ k) * sum; % 迭代计算出新的T2^k        
        T2k_array(k + 1,1) = k;  
        T2k_array(k + 1,2) = T2k;
        T2k_array(k + 1,3) = abs(T2k - tmp);    
        tmp = T2k; % 更新tmp的值为当次迭代结果
        k = k + 1; % 迭代次数加1
        %if(k > 23)
        %    throw(MException('MATLAB:TooManyInputs',['已迭代' num2str(k) '次，过多迭代次数会导致计算十分耗时']));
        %end
    end

    T2k_array = T2k_array(1:k,:); % 去除尾部无用的零
end