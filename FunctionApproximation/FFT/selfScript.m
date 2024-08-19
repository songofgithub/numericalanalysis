%% 测试自带fft函数
close all;
N = 1024;
X = [1,2,3,4,5,6,7];
plot([-2*pi:2*pi/N:2*pi-2*pi/N],[real(fft(X,N)),real(fft(X,N))]);
xlim([-2*pi 2*pi]);
set(gca,'XTick',[-2*pi:pi/2:2*pi]);
set(gca,'xtickLabel',{'-2π','-3π/2','-π','-π/2','0','π/2','π','3π/2','2π'});

%% 测试get_W函数
close all;
N = 1024;
y = get_W(N);
plot([0:2*pi/N:2*pi-2*pi/N],imag(y));
xlim([0 2*pi]);
set(gca,'XTick',[0:pi/2:2*pi]);
set(gca,'xtickLabel',{'0','π/2','π','3π/2','2π'});

%% 测试S_fft函数
close all;
N = 1024;
X = [1,2,3,4,5,6,7];
y = S_fft(X,N);
plot([-2*pi:2*pi/N:2*pi-2*pi/N],[real(y),real(y)]);
xlim([-2*pi 2*pi]);
set(gca,'XTick',[-2*pi:pi/2:2*pi]);
set(gca,'xtickLabel',{'-2π','-3π/2','-π','-π/2','0','π/2','π','3π/2','2π'});