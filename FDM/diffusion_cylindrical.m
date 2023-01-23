%%
% solution of the equation (dc/dt) = D*(d^2c/dr^2) + (D/r)*(dc/dr) , 0<r<1
% dc/dr = 0 at r = 0 and dc/dr = k at r = 1;

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

C = zeros(n, 1);
C_new = zeros(n, 1);

total_iterations = 0.1/dt ;

C_transient = zeros(total_iterations, n);




for iter=1:total_iterations
    for i = 1:n
        r = a + (i-1)*h ;
        if i == 1
            C_new(i) = C(i) + D*dt*2*(C(i+1) - C(i))/(h^2) ;
        elseif i == n
            C_new(i) = C(i) + D*dt*2*( h*k - C(i) + C(i-1) )/(h^2) + D*dt*k/r ;
        else
            C_new(i) = C(i) + dt*D*(C(i+1) - 2*C(i) + C(i-1))/(h^2) + dt*D*(C(i+1) - C(i-1))/(2*r*h) ;
        end
    end


    C = C_new ;
    C_transient(iter, 1:n) = C;
    plot(r_arr, C, linewidth=2);
    ylim([0, 3.5])
%     pause(0.2);

end

hold on;
plot(r_arr, C, linewidth=2);


%% solution for displacement u
% solution of the Differential Equation: (d^2u/dr^2) + (2/r)*(du/dr) - (2/r^2)*u = (omega*(1+nu)/(3*(1-nu)))*(dc/dr)
% boundary conditions: sigma_r (at r = R) = 0 and sigma_r (at r = 0) is finite (here the equation is in dimensional form)

nu = 0.33;
U = zeros(n, 1) ;
U(1) = 0;
U(n) = 0;

c_charactersitic = 10e27;
omega_characteristic = 10e-27;
E0 = 10e9;

beta = omega_characteristic*c_charactersitic ;

E = 1;
omega = 1;
sigma_r = zeros(n, 1);
sigma_r(n) = 0;

for iter = 1:total_iterations
    for i = 2:(n-1)

        A = beta*(omega*(1+nu)/(3*(1-nu)))*(C_transient(iter, i+1) - C_transient(iter, i-1))/(2*h);
        r = a + (i-1)*h;

        U(i) = (U(i+1)*(1+h/r) + U(i-1)*(1-h/r)-A*h*h)/(2 + 2*h*h/(r*r)) ;
        epsilon_r = (U(i+1) - U(i-1))/(2*h);
        epsilon_theta = U(i)/r ;

        sigma_r(i) = E*((1-nu)*epsilon_r + 2*nu*epsilon_theta - beta*omega*C_transient(iter, i)*(1+nu)/3)/(1-nu-2*nu*nu)  ;
    end
end

figure ;
plot(r_arr, U, LineWidth=2) ;
figure ;
plot(r_arr, sigma_r, linewidth=2);













