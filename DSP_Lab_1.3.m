clc,clear
x = [1 1 1 1 3 3 3 3 1 1 1 2 2 2 2 1 1 1 1];
y = zeros(1,22);
w_0 = [x 0 0 0];
w_1 = zeros(1,22);
w_2 = zeros(1,22);
w_3 = zeros(1,22);
for n= 1:22
    if (n-1)<1
        w_1(n) = 0;
    else
        w_1(n) = w_0(n-1);
    end
    if (n-2) <1
        w_2(n) = 0;
    else
        w_2(n) = w_1(n-1);
    end
    if (n-3) <1
        w_3(n) = 0;
    else
        w_3(n) = w_2(n-1);
    end
    y(n) = w_0(n) + w_2(n) - 2*w_3(n);
end
j = 0 : 21;
stem(j, y)