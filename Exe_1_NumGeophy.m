%% Exerciese Numerical Modelling Module 1
%  Jiahui Kang
%  Shiyi Li
%% ------ Task 1 Cancellation ------
clc; clear all;
a = 1:10;
n1 = 10.^a;
n2 = 10.^a-10.^(-a);
n3 = n1-n2;

rel.error = abs((n3./(10.^(-a))-1)*100);

n1_sin = single(10.^a);
n2_sin = single(10.^a-10.^(-a));
n3_sin = single(n1_sin-n2_sin);
rel_sin.error = abs((n3_sin./(10.^(-a))-1)*100);

figure(1)
semilogy(a,rel.error);
hold on
semilogy(a,rel_sin.error);
legend('Double precision','Single precision','location','east');
title('Resulting Curve of Double & Single Precision');
xlabel('a');
ylabel('rel.error');
hold off
%% ------ Task 2 Numerical differentiation ------
clc; clear all;

% Find Analytical first derivative
syms x
f(x) = 3*cos(10*x)+5*sin(15*x);
df = diff(f(x),x);

% Find Central finite difference approximation

h = 0.01;
x_aprx = -pi:h:pi;
l = length(x_aprx);


df_function = @(x) 75*cos(15*x) - 30*sin(10*x);
F = f(x_aprx);
df_true = df_function(x_aprx);

df_cfd = zeros(1,l);
df_cfd(1) = df_true(1); % Boundary condition is not given.
df_cfd(l) = df_true(l); % Boundary condition is not given.
for i = 2:l-1
    df_cfd(i) = (f(x_aprx(i+1))-f(x_aprx(i-1)))/(2*h);
end

figure(2)
plot(x_aprx,df_cfd);
hold on;

% Determin coefficients for the fourth-order difference

h = 0.01;
M = [1 1 1 1 1;
    -2 -1 0 1 2;
    2 1/2 0 1/2 2;
    -4/3 -1/6 0 1/6 4/3
    2/3 1/24 0 1/24 2/3];
u = [0; 1/h; 0; 0; 0];
A = M\u;

% Compute the first derivative
df_fdm = zeros(1,l);

df_fdm(1) = df_true(1);
df_fdm(2) = df_true(2);
df_fdm(l-1) = df_true(l-1);
df_fdm(l) = df_true(l);% Unknown boundary conditions

for i = 3:l-2
    df_fdm(i) = dot(A,F(i-2:i+2));
end

plot(x_aprx,df_fdm,'*');

% Check the accuracy of the solutions

error_cfd = df_true - df_cfd;
error_fdm = df_true - df_fdm;

figure(3)
plot(x_aprx, error_cfd);
hold on
plot(x_aprx, error_fdm);
hold off

%% ------ Task 3 Interpolation ------

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

figure(4)
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


