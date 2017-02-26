clc
clear all
x = 1:10;
y = [0 0 0 -0.5 1 0.55 0 0 0 0];
xx = linspace(1,10,1000);
f1 = interp1(x,y,xx,'linear');
f2 = interp1(x,y,xx,'spline');
c1 = polyfit(x,y,3);
f_1 = polyval(c1,xx);
c2 = polyfit(x,y,5);
f_2 = polyval(c2,xx);
c3 = polyfit(x,y,7);
f_3 = polyval(c3,xx);
hold on;
plot(x,y,'*');
plot(xx,f1);
plot(xx,f2);
plot(xx,f_1);
plot(xx,f_2);
plot(xx,f_3);
legend('original','linear','spline','polynomial 3','polynomial 5','polynomial 7');
title('results with different interpolation methods')
xlabel('x');
ylabel('y');