% solution of the ODE y'' + 2y' + y = x^2 ;
% with the boundary conditions y(a) = y1 and y(b) = y2 ;

%%
clear;

clc;
a = 0;
b = 1;
y1 = 0.2;
y2 = 0.8;
n = 50; % number of points 



h = (b-a)/(n-1) ;

y(1) = y1;
y(n) = y2;

y_new(1) = y1;
y_new(n) = y2;

iterations = 0;

error_now = 1;
error_req = 10e-6;

%%
while error_now > error_req
    for i=2:(n-1)
        x = a + (i-1)*h;
        y_new(i) = ((h*x)^2 +  y(i-1)*(h-1) - y(i+1)*(h+1) )/(h*h - 2) ;
    end

    iterations = iterations + 1;
    error_now = 0;
    for i=2:(n-1)
        error_now = error_now + abs(y_new(i) - y(i)) ;
    end
    
    if(rem(iterations, 1000) == 0)
        iterations
        error_now
    end


    y = y_new ;
end


%% plotting my solution

x_arr = a:h:b ;
plot(x_arr, y, LineWidth=2) ;
title('My Solution');
xlim([a, b]);
ylim([y1, y2]) ;
disp(y) ;

%% plotting the solution using in buillt function to compare

syms ys(x)
ode = diff(ys, x, 2) + 2*diff(ys, x) + ys == x^2 ;
cond1 = ys(a) == y1;
cond2 = ys(b) == y2;

ysol(x) = dsolve(ode, [cond1, cond2]) ;

figure;
fplot(x, ysol, linewidth=2);
title('Matlab Solution')
xlim([a, b]);
ylim([y1, y2]) ;

