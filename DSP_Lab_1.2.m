clc,clear
L = 5 ;
x = [1 1 1 1 3; 3 3 3 1 1; 1 2 2 2 2; 1 1 1 1 0];
y = zeros(1,24);
for i =0 : 3
    t = 0 : 4;
    h = dirac(t) + dirac(t-2);
    h = 1 * sign(h)-2 * sign(dirac(t-3));
    mid = [zeros(1,i*L) conv(x(i+1,:),h) zeros(1,15-i*L)];
    y = y + mid;
end
j = 0 : 23;
fprintf('The Output:');
stem(j, y)
