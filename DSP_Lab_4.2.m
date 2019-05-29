%DSP_Lab_4.2.m
figure(1)
subplot(2,2,1);
%5th-order
[z,p,k] = buttap(5);          % Butterworth filter prototype
[num,den] = zp2tf(z,p,k);     % Convert to transfer function form
[H,w]=freqs(num,den);         % Frequency response of analog filter
Hf=abs(H);
Hf2=power(Hf,2);
f=w/(2*pi);
plot(f,Hf2);
grid on
ylabel('H|f|^2')
xlabel('f')
title('5-order');
%10th-order
subplot(2,2,2);
[z,p,k] = buttap(10);          % Butterworth filter prototype
[num,den] = zp2tf(z,p,k);     % Convert to transfer function form
[H,w]=freqs(num,den);         % Frequency response of analog filter
Hf=abs(H);
Hf2=power(Hf,2);
f=w/(2*pi);
plot(f,Hf2);
grid on
ylabel('H|f|^2')
xlabel('f')
title('10-order');
%20th-order
subplot(2,2,3);
[z,p,k] = buttap(20);          % Butterworth filter prototype
[num,den] = zp2tf(z,p,k);     % Convert to transfer function form
[H,w]=freqs(num,den);         % Frequency response of analog filter
Hf=abs(H);
Hf2=power(Hf,2);
f=w/(2*pi);
plot(f,Hf2);
grid on
ylabel('H|f|^2')
xlabel('f')
title('20-order');
%2th-order
subplot(2,2,4);
[z,p,k] = buttap(2);          % Butterworth filter prototype
[num,den] = zp2tf(z,p,k);     % Convert to transfer function form
[H,w]=freqs(num,den);         % Frequency response of analog filter
Hf=abs(H);
Hf2=power(Hf,2);
f=w/(2*pi);
plot(f,Hf2);
grid on
ylabel('H|f|^2')
xlabel('f')
title('2-order');

