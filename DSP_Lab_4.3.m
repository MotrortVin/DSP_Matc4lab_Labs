%DSP_Lab_4.3.m
clc;
close all;
clear all;

%Parameters
Ap_ = 0.1;     % Min Passband Ripple (db)
Aa_ = 41;       % Min Stopband Attenuation (db)
Wpa = 0.4*pi;      % Lower passband Freq (rad/s)
Wsa = 0.45*pi;      % Lower stopband Freq (rad/s)
Wsb = 0.65*pi;      % Higher stopband Freq (rad/s)
Wpb = 0.7*pi;      % Higher stopband Freq (rad/s)
Ws  = 1;     % Sampling Freq (rad/s)

Df = min(Wsa-Wpa, Wpb-Wsb)/2/pi;
W1 = (Wpa+Wsa)/2;
W2 = (Wsb+Wpb)/2;

% 2. Choose Delta

delta_a = 10^(-Aa_/20);
c = 10 ^ ( Ap_ / 20);
delta_p = (c-1)/(c+1);
delta = min(delta_a, delta_p);

% 3. Get Aa from delta
Aa = -20*log10(delta);            % Actual stopband attenuation

% 4. Calculate alpha
if (Aa <= 21)
    alpha = 0;
elseif ((21 < Aa) && (Aa <= 50))
    alpha = 0.5842*(Aa-21)^0.4 + 0.07886*(Aa-21);
else
    alpha = 0.1102*(Aa - 8.7);
end

% 5. Calculate D and N

if (Aa <= 21)
    D = 0.9222;
else
    D = (Aa - 7.95)/14.36;
end

N = ceil ( Ws * D / Df +1);

if (mod(N,2) == 0)
    N = N + 1;
end

% 6. Form Kaiser window

M = (N-1)/2;
n = 0:1:N-1;
wn = kaiser(N,alpha)';
alpha_1=(N-1)/2;
m=n-alpha_1+eps;
a1=sin(pi*m)./(pi*m);
a2=sin(W1*m)./(pi*m);
a3=sin(W2*m)./(pi*m);

dn = a1 + a2 -a3;
h = wn.*dn;
figure(1)
stem(n,h);
title('Impulse Response h(n)');
xlabel('n');
ylabel('h(n)');
[x,k] = freqz(h,1,512);
figure(2)
plot(k/pi,abs(x));
title('kaiser带阻滤波器的幅频特性');
xlabel('ω/π');
ylabel('H|e^j^w|');