function [wOut] = get_W(N)
%GET_W 获取N点DFT的旋转因子W的值
%   N:DFT的点数
unitAngle = 2*pi/N;
wOut = zeros(1,N);
for loop=0:N-1
    wAngle = loop*unitAngle;
    wOut(loop+1) = cos(wAngle)-sin(wAngle)*1i;
end
end

