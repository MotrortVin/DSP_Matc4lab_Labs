%DSP_Lab_4.1.m
num=[1,0,7,-1,8,0]; %The last four digits are 0,7,1 and 8.
mot=[1,0.22,0.037,0.142,-0.107,-0.013];
figure(1);  
subplot(2,1,1); 
impz(num,mot,50); 
subplot(2,1,2); 
stepz(num,mot,50);
