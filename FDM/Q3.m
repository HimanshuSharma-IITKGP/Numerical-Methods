% solution of d^2u/dx^2 = k (constant) using FDM
% with boundary conditions u(a) = u0 and u'(b) = p0 using FDM

%%
clear;
close all;
clc;

a = 0;
b = 1;
u0 = 1;
p0 = 1;
k = 5;
n = 30;

x_domain = b - a;
h = x_domain/(n-1) ;
x = a:h:b ; 


u(1) = u0;
u(n) = 0;

u_new(1) = u0;
u_new(n) = 0;

error_now = 1 ;
error_req = 1e-7 ;

iterations = 0 ;

%%
while error_now > error_req
    for i=2:n
        if(i ~= n)
            u_new(i) = 0.5*(u(i-1) + u(i+1)) - 0.5*k*h*h  ;
        else
            u_new(i) = u(i-1) + p0*h - 0.5*k*h*h ;
        end
    end

    error_now = 0;
    for i=2:n
        error_now = error_now + abs(u(i) - u_new(i)) ;
    end
    iterations = iterations + 1;
    u = u_new ;
end

%%
plot(x, u, linewidth=2 );
xlim([0, 1]);
ylim([-0.6, 1]) ;
title('My Solution')

%% plotting the exact solution

syms us(xs)
dus(xs) = diff(us, xs) ;

ode = diff(us, xs, 2) == k;

cond1 = us(a) == u0;
cond2 = dus(b) == p0;

usol(xs) = dsolve(ode, [cond1, cond2]) ;
figure ;
fplot(xs, usol, linewidth=2) ;
xlim([0, 1]) ;
ylim([-0.6, 1]);
title('Matlab Solution')

