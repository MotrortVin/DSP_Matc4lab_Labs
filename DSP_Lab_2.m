clc,clear
figure(1);
r = 0.98;
p = r^(1/12);
num = [1,-(exp((pi/6)*1i)+exp((-pi/6)*1i))*p,p*p,0,0,0,0,0,0,0,0,0,-1,(exp((pi/6)*1i)+exp((-pi/6)*1i))*p,-p*p];
den = [1,-(exp((pi/6)*1i)+exp((-pi/6)*1i)),1,0,0,0,0,0,0,0,0,0,-r,(exp((pi/6)*1i)+exp((-pi/6)*1i))*r,-r];
subplot(1,2,1)
zplane(num,den);
title('零极点图');
w = 0:pi;
[h,w] = freqz(num,den,1024*1024);
H = abs(h);
subplot(1,2,2);
plot(w/pi,H);
xlabel('Frequence(Hz)')
title('Magnitude Response');
fs = 600;
N = 1024*1024;
n = 1:N;
t = n/fs;
x = 1*sin(2*pi*50*t)+0.45*sin(4*pi*50*t)-0.6*sin(6*pi*50*t)+1.2*sin(8*pi*50*t)-0.5*sin(10*pi*50*t)+0.8;
T = 0:1/fs:6000/fs;
figure(2);
subplot(2,2,1);
plot(T,x(1:6001));
xlabel('Time(s)');
title('输入信号时域曲线');
subplot(2,2,3);
plot(T,x(1:6001));
xlabel('Time(s)');
title('放大');
axis([5.0,5.05,-1.5,3])
Xtmp = fft(x,N);
X = abs(Xtmp)*2/N; 
f = n*fs/N;
subplot(2,2,[2 4]);
plot(f,X);
xlabel('Frequence(Hz)')
title('输入信号频域曲线')
y=filter(num,den,x);
figure(3);
subplot(2,2,1);
plot(T,y(1:6001));
xlabel('Time(s)');
title('输出信号时域曲线');
subplot(2,2,3);
plot(T,y(1:6001));
xlabel('Time(s)');
title('放大');
axis([5.0,5.05,-2,2])
Ytmp=fft(y,N);
Y = abs(Ytmp)*2/N; 
subplot(2,2,[2 4]);
plot(f,Y);
xlabel('Frequence(Hz)')
title('输出信号频域曲线')
