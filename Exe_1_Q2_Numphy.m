%% ------ Task 2 Numerical differentiation ------
clc; clear all;

% Find Analytical first derivative
syms x
f(x) = 3*cos(10*x)+5*sin(15*x);
df = diff(f(x),x);
df_function = matlabFunction(df);

% Determin coefficients for the fourth-order difference
M = [1 1 1 1 1;
     -2 -1 0 1 2;
     2 1/2 0 1/2 2;
     -4/3 -1/6 0 1/6 4/3
     2/3 1/24 0 1/24 2/3];
u = [0; 1; 0; 0; 0];
A = M\u;

h = 0.01;
%     j = 1;
%     error_cfd = zeros(5);
%     error_fdm = zeros(5);
    
% for h = [0.1 0.01 0.001 0.001 0.0001]
x_aprx = -pi:h:pi;
l = length(x_aprx);
    
% Find analytical first derivative
F = f(x_aprx);
df_true = df_function(x_aprx);

% Find Central finite difference approximation
df_cfd = zeros(1,l);
df_cfd(1) = df_true(1); % Boundary condition is not given.
df_cfd(l) = df_true(l); % Boundary condition is not given.
    for i = 2:l-1
        df_cfd(i) = (f(x_aprx(i+1))-f(x_aprx(i-1)))/(2*h);
    end

figure(1)
set (gcf,'Position',[400,100,1280,840], 'color','w');
plot(x_aprx,df_true,'b');
hold on;
plot(x_aprx,df_cfd,'r--','linewidth',2);

% Compute the first derivative
df_fdm = zeros(1,l);

df_fdm(1) = df_true(1);
df_fdm(2) = df_true(2);
df_fdm(l-1) = df_true(l-1);
df_fdm(l) = df_true(l);     % Unknown boundary conditions

for i = 3:l-2
    df_fdm(i) = dot(A/h,F(i-2:i+2));
end

plot(x_aprx,df_fdm,'*');
legend('Analytical Solution','Central Finite Difference','Fourth-Order Finite Difference');
str_ttl = ['First Derivative Solutions(\Delta x = ', num2str(h),')'];
title(str_ttl,'fontsize',14);
xlabel('x','fontsize',16);
ylabel('\partial f/\partial x','fontsize',16);

% Check the accuracy of the solutions

error_max_cfd = max(df_true - df_cfd);
error_max_fdm = max(df_true - df_fdm);

dim = [.15 .3 .6 .6];
str = {['error\_max\_cfd = ', num2str(error_max_cfd)],...
       ['error\_max\_fdm = ', num2str(error_max_fdm)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

%     j = j+1;
% end
% 
% 
% figure(3)
% plot(x_aprx, error_order_cfd);
% hold on
% plot(x_aprx, error_order_fdm);
% hold off