% solution of the equation (dc/dt) = D*(d^2c/dr^2) + (2*D/r)*(dc/dr) , 0<r<1
% dc/dx = 0 at x = 0 and dc/dx = k at x = 1;

clear;
clc;
close all ;

a = 0;
b = 1;
c1 = 0;
c2 = 1;
n = 51;
k = 1;

h = (b-a)/(n-1) ;
dt = 0.0002;
r_arr = a:h:b ;

L_c = 1e-6;
D_c = 1e-12;
T_c = (L_c^2)/D_c ;
D = 1;
total_iterations = 0.1/dt ;


C = zeros(n, 1);
C_transient = zeros(total_iterations, n) ;



for iter=1:total_iterations
    for i = 1:n
        r = a + (i-1)*h ;
        if i == 1
            C(i) = C(i) + D*dt*2*(C(i+1) - C(i))/(h^2) ;
        elseif i == n
            C(i) = C(i) + D*dt*2*( h*k - C(i) + C(i-1) )/(h^2) + 2*D*dt*k/r ;
        else
            C(i) = C(i) + dt*D*(C(i+1) - 2*C(i) + C(i-1))/(h^2) + 2*dt*D*(C(i+1) - C(i-1))/(2*r*h) ;
        end

         C_transient(iter, :) = C ;
    end
    
%     plot(r_arr, C, linewidth=2);
%     ylim([0, 3.5])
%     pause(0.2);

end

hold on
plot(r_arr, C, linewidth=2);




