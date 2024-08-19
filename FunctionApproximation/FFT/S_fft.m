function [Xk] = S_fft(X,N)
%S_FFT 自定义fft算法（按时域奇偶抽取，基2-DIT-fft）
%   X:输入序列
%   N:运算点数(必须为2的幂)
    xSize = size(X);
    if(xSize(1)~=1)
        X = X'; % 统一为行向量
    end
    xSize = size(X);

    if(xSize(2)<N)
        X = [X,zeros(1,N-xSize(2))]; % 输入数据长度不够时补0
    end

    L = log2(N); % 蝶形算法层数（蝶形次数）
    Xk = zeros(1,N); % 初始化输出    
    for j=0:N-1
        Xk(j+1) = X(bin2dec(reverse(dec2bin(j,L)))+1); % 按位倒序重排输入序列
    end
    WN = get_W(N); % 获取Wnk的值
    for i=1:L % 按需要计算的总蝶形层数每层计算更新输出Xk的内容
        s = 2^(i-1); % 蝶形结对偶节点的距离（作为计算地址的偏移量）
        s1 = s*2; % 跳转到同一层下一个蝶形结的距离
        s2 = N/s1; % 每层需要计算的蝶形结个数       
        for j=0:s2-1                        
            kStart = j*s1; % 每层单个蝶形结的起始地址           
            for k=0:s-1
                addr = kStart+k+1;
                addrDuality = addr + s;                
                addrData = Xk(addr);
                addrDualityData = Xk(addrDuality);
                r = addr*(2^(L-i));
                WNr = WN(mod(r,N));      

                % Xm(k) = Xm-1(k) + Xm-1(k+2^(m-1))WN(r)
                % Xm(k) = Xm-1(k) - Xm-1(k+2^(m-1))WN(r)
                Xk(addr) = addrData+addrDualityData*WNr;
                Xk(addrDuality) = addrData-addrDualityData*WNr;
            end
        end
    end
end

