% Solution of the ODE d^2u/dx^2 = k*u ;
% with the boundary condition u(a) = u1 and u(b) = u2 using FDM

%%
clear;
clc;

a = 0;
b = 1;
u1 = 0;
u2 = 1;
k = 5;
n = 50;

h = (b-a)/(n-1) ;
x_arr = a:h:b;


u(1) = u1;
u(n) = u2;

u_new(1) = u1;
u_new(n) = u2;

error_now = 1;
error_req = 1e-6;

iterations = 0;

%%
while error_now > error_req
    for i=2:(n-1)
        u_new(i) = (u(i+1) + u(i-1))/(h*h*k + 2) ;
    end

    error_now = 0;
    for i=2:(n-1)
        error_now = error_now + abs(u(i) - u_new(i)) ;
    end

    iterations = iterations + 1;

    if(rem(iterations, 1000) == 0)
        iterations
        error_now
    end
    u = u_new ;
end

%%

subplot(2, 1, 1);
plot(x_arr, u, linewidth=2);
title('My Solution');
xlim([0, 1]);
ylim([-0.5, 1.5]) ;


%%

syms us(xs)

ode = diff(us, xs, 2) == k*us;
cond1 = us(a) == u1;
cond2 = us(b) == u2;

usol(xs) = dsolve(ode, [cond1, cond2]) ;
subplot(2, 1, 2) ;
title('Matlab Solution')
fplot(xs, usol, 'r', linewidth=2) ;
xlim([0, 1]);
ylim([-0.5, 1.5]) ;

