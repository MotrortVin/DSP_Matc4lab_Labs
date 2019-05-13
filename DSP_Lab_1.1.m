clc,clear
x = [1 1 1 1 3 3 3 3 1 1 1 2 2 2 2 1 1 1 1];
h = [1 0 1 -2];
fprintf('The result is:');
y = conv(x, h)
