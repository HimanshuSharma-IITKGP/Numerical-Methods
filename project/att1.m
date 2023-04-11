%%
clear;
clc;


D_0 = 16;
c_ref = 53.98e27;
E_0 = 90.13e9;
E_b_outer = 1.5*E_0;
R_0 = 10e-6;
j_outer = 10e34;
omega = 0.014e-27 ;
nu = 0.28;
xi = 4;

h = 0.05;
r_tilde = 0:h:1 ;
dt = 0.001;
t_final = 1;
n_points = (1/h) + 1; 
n_iterations = t_final/dt ;

%%
c_tilde = zeros(n_iterations, n_points) ;
c_tilde(1, :) = 0; % initial condition




for iter = 1:n_iterations-1
    for i = 1:n_points
        if(i == 1) % boundary condition on left
            c1 = c_tilde(iter, i+1);
            c_tilde(iter+1, i) =  c_tilde(iter, i) + dt*( (c_tilde(iter, i+1) - c1)/(2*h*h) + (c_tilde(iter, i+1) - 2*c_tilde(iter, i) + c1)/(h^2) ) ;
        elseif(i == n_points) % boundary condition on right
            c1 = c_tilde(iter, i-1) + 2*h*xi*(1-c_tilde(iter, i)) ;
            c_tilde(iter+1, i) = c_tilde(iter, i) + dt*( (c1 - c_tilde(iter, i-1))/(2*h*r_tilde(i)) + (c1 - 2*c_tilde(iter, i) + c_tilde(iter, i-1))/(h^2) ) ;
        else
            c_tilde(iter+1, i) = c_tilde(iter, i) + dt*( (c_tilde(iter, i+1) - c_tilde(iter, i-1))/(2*h*r_tilde(i)) + (c_tilde(iter, i+1) - 2*c_tilde(iter, i) + c_tilde(iter, i-1))/(h^2) ) ;
        end

        if(isnan(c_tilde(iter+1, i)))
            disp(iter);
            disp(i) ;
            break
        end

%         pause(0.1);
%         plot(r_tilde, c_tilde(iter+1, :), linewidth=2)

        
    end
end

plot(r_tilde, c_tilde(n_iterations, :), linewidth=2);
title('Final Concentration', 'Color','blue');


%%

alpha = omega*c_ref*(1+nu)/(3*(1-nu)) ;

u_tilde = zeros(n_iterations, n_points) ;
u_tilde(:, 1) = 0;
u_tilde(:, n_points) = 1 ;



for iter = 1:n_iterations-1
    for i = 2:n_points-1
         l = 1/r_tilde(i)^2 + 2/h^2 ;
         u_tilde(iter+1, i) = ( (u_tilde(iter, i+1) +  u_tilde(iter, i-1))/(h^2) +  (u_tilde(iter, i+1)- u_tilde(iter, i-1))/(2*h*r_tilde(i)) - alpha*(c_tilde(iter, i+1)-c_tilde(iter, i-1))/(2*h) )/l ;
    
%          plot(r_tilde, u_tilde(iter, :), linewidth=2);
%          pause(0.01);
    end
end

figure
plot(r_tilde, u_tilde(n_iterations, :), linewidth=2) ;
title('Displacement', 'Color', 'blue')

%%
sigma_rr_tilde = zeros(1, n_points) ;
sigma_tt_tilde = zeros(1, n_points) ;
epsilon_rr = u_tilde(n_iterations, n_points)/r_tilde(n_points);
epsilon_tt = epsilon_rr ;

for i = 1:n_points
    a = c_tilde(n_iterations, i)*omega*c_ref/3 ;
    sigma_rr_tilde(i) = (nu*epsilon_tt + (1-nu)*epsilon_rr - a*(1+nu))/((1+nu)*(1-2*nu)) ;
    sigma_tt_tilde(i) = (nu*epsilon_rr + (1-nu)*epsilon_tt - a*(1+nu))/((1+nu)*(1-2*nu)) ;
end

figure
plot(r_tilde, sigma_tt_tilde, linewidth=2) ;
title('Hoop Stress', 'Color', 'Blue')
figure
plot(r_tilde, sigma_tt_tilde, linewidth=2) ;
title('Radial Stress', 'Color', 'Blue');

