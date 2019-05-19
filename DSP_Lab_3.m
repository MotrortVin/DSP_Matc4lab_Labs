L0=100;
L1=200;
n0=0:L0-1;
n1=0:L1-1;
f0=50;
fs=1000;
x_Rect1=cos(2*pi*f0/fs*n0);
x_Rect2=cos(2*pi*f0/fs*n1);
x_Hamm1=(0.54-0.46*cos(2*pi*n0/(L0-1))).*cos(2*pi*f0/fs*n0);
x_Hamm2=(0.54-0.46*cos(2*pi*n1/(L1-1))).*cos(2*pi*f0/fs*n1);
% Fig.9.1.5
figure(1);
subplot(2,1,1);
plot(n0,x_Rect1);
axis([0 200 -2 2]);
xlabel('time samples n(L=100)');
title('Rectangular Window');
subplot(2,1,2);
plot(n1,x_Rect2);
axis([0 200 -2 2]);
xlabel('time samples n(L=200)');
title('Rectangular Window');

%Fig.9.1.6
figure(2);
subplot(2,1,1);
plot(n0,x_Hamm1);
axis([0 200 -2 2]);
xlabel('time samples n(L=100)');
title('Hamming Window');
subplot(2,1,2);
plot(n1,x_Hamm2);
axis([0 200 -2 2]);
xlabel('time samples n(L=200)');
title('Hamming Window');

%Fig.9.1.7
N=2000;
k=0:1999;
x_Rect1=fft(x_Rect1,N);
x_Hamm1=fft(x_Hamm1,N);
w=k/N;
figure(3);
subplot(1,2,1);
plot(w*2,abs(x_Rect1));
axis([0 0.2 0 100]);
xlabel('w in units of pi');
title('Magnitude Spectra, L=100');
hold on;
plot(w*2,abs(x_Hamm1));
axis([0 0.2 0 100]);

x_Rect2=fft(x_Rect2,N);
x_Hamm2=fft(x_Hamm2,N);
subplot(1,2,2);
plot(w*2,abs(x_Rect2));
axis([0 0.2 0 100]);
xlabel('w in units of pi');
title('Magnitude Spectra, L=200');
hold on;
plot(w*2,abs(x_Hamm2));
axis([0 0.2 0 100]);

f1=2000;
f2=2500;
f3=3000;
fs=10000;
L_10=10;
L_20=20;
L_40=40;
L_100=100;
a0=0:L_10-1;
a1=0:L_20-1;
a2=0:L_40-1;
a3=0:L_100-1;

%Fig.9.1.8
x_R10=cos(2*pi*f1/fs*a0)+cos(2*pi*f2/fs*a0)+cos(2*pi*f3/fs*a0);
X_R10=fft(x_R10,N);
figure(4);
subplot(2,2,1);
plot(w,abs(X_R10));
axis([0 1 0 50]);
xlabel('f/fs');
title('Rectangular Window, L=10');

x_H10=(cos(2*pi*f1/fs*a0)+cos(2*pi*f2/fs*a0)+cos(2*pi*f3/fs*a0)).*(0.54-0.46*cos(2*pi*a0/(L_10-1)));
X_H10=fft(x_H10,N);
subplot(2,2,2);
plot(w,abs(X_H10));
axis([0 1 0 50]);
xlabel('f/fs');
title('Hamming Window, L=10');

x_R20=cos(2*pi*f1/fs*a1)+cos(2*pi*f2/fs*a1)+cos(2*pi*f3/fs*a1);
X_R20=fft(x_R20,N);
subplot(2,2,3);
plot(w,abs(X_R20));
axis([0 1 0 50]);
xlabel('f/fs');
title('Rectangular Window, L=20');

x_H20=(cos(2*pi*f1/fs*a1)+cos(2*pi*f2/fs*a1)+cos(2*pi*f3/fs*a1)).*(0.54-0.46*cos(2*pi*a1/(L_20-1)));
X_H20=fft(x_H20,N);
subplot(2,2,4);
plot(w,abs(X_H20));
axis([0 1 0 50]);
xlabel('f/fs');
title('Hamming Window, L=20');

%Fig.9.1.9
x_R40=cos(2*pi*f1/fs*a2)+cos(2*pi*f2/fs*a2)+cos(2*pi*f3/fs*a2);
X_R40=fft(x_R40,N);
figure(5);
subplot(2,2,1);
plot(w,abs(X_R40));
axis([0 1 0 50]);
xlabel('f/fs');
title('Rectangular Window, L=40');

x_H40=(cos(2*pi*f1/fs*a2)+cos(2*pi*f2/fs*a2)+cos(2*pi*f3/fs*a2)).*(0.54-0.46*cos(2*pi*a2/(L_40-1)));
X_H40=fft(x_H40,N);
subplot(2,2,2);
plot(w,abs(X_H40));
axis([0 1 0 50]);
xlabel('f/fs');
title('Hamming Window, L=40');

x_R100=cos(2*pi*f1/fs*a3)+cos(2*pi*f2/fs*a3)+cos(2*pi*f3/fs*a3);
X_R100=fft(x_R100,N);
subplot(2,2,3);
plot(w,abs(X_R100));
axis([0 1 0 50]);
xlabel('f/fs');
title('Rectangular Window, L=100');

x_H100=(cos(2*pi*f1/fs*a3)+cos(2*pi*f2/fs*a3)+cos(2*pi*f3/fs*a3)).*(0.54-0.46*cos(2*pi*a3/(L_100-1)));
X_H100=fft(x_H100,N);
subplot(2,2,4);
plot(w,abs(X_H100));
axis([0 1 0 50]);
xlabel('f/fs');
title('Hamming Window, L=100');