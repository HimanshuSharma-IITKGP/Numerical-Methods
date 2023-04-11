%%
clear;
clc;
close all;

D_0 = 16;
c_ref = 53.98e27;
E_0 = 90.13e9;
E_b_outer = 1.5*E_0;
R_0 = 10e-6;
j_outer = 10e34;
T = 300;
R_g = 8.314; 
omega = 0.014e-27 ;
nu = 0.28;
xi = 4;
eta_E = -0.1464 ;
chi_max = 4.4 ;

h = 0.05;
r_tilde = 0:h:1 ;
dt = 0.001;
t_final = 1;
n_points = (1/h) + 1; 
n_iterations = t_final/dt ;


c_tilde = zeros(n_iterations, n_points);
u_tilde = zeros(n_iterations, n_points);
u_tilde(:, n_points) = 1;
sigma_rr_tilde = zeros(n_iterations, n_points);
sigma_tt_tilde = zeros(n_iterations, n_points);


sigma_h_tilde = zeros(n_iterations, n_points);




%%
alpha = omega*c_ref*(1+nu)/(3*(1-nu)) ;
b = omega*E_0/(R_g*T) ;


for iter = 1:n_iterations-1
    for i = 1:n_points
        l = 1/r_tilde(i)^2 + 2/h^2 ;
        a = omega*c_ref*c_tilde(iter, i)/3 ;
        beta = b*(1 + eta_E*chi_max*c_tilde(iter, i)) ;

        if (i==1)
            c1 = c_tilde(iter, i+1) ;
            c_tilde(iter+1, i) = c_tilde(iter, i) + dt*  ( (c_tilde(iter, i+1) - 2*c_tilde(iter, i) + c1 )/(h^2) ...
                                                         + (c_tilde(iter, i+1) - c1)/(2*h*h)  ...
                                                         - beta*c_tilde(iter, i)*(sigma_h_tilde(iter, i+1) - sigma_h_tilde(iter, i))/(h*h)  ...
                                                         - beta*c_tilde(iter, i)*(sigma_h_tilde(iter, i+2) - 2*sigma_h_tilde(iter, i+1) + sigma_h_tilde(iter, i))/(h^2) ) ;


            u_tilde(iter+1, i) = 0;

            
            epsilon_tt =  (u_tilde(iter, i+1)-u_tilde(iter, i))/(h)  ;
            epsilon_rr =  epsilon_tt ;
            
            sigma_h_tilde(iter+1, i) = (epsilon_tt + epsilon_rr - 2*a*(1+nu))/(3*(1-2*nu)) ;
            sigma_rr_tilde(iter+1, i) = (nu*epsilon_tt + (1-nu)*epsilon_rr - a*(1+nu))/((1+nu)*(1-2*nu)) ;
            sigma_tt_tilde(iter+1, i) = (nu*epsilon_rr + (1-nu)*epsilon_tt - a*(1+nu))/((1+nu)*(1-2*nu)) ;

        elseif(i==n_points)
            c2 = c_tilde(iter, i-1) + 2*h*xi*(1-c_tilde(iter, i)) ;
            c_tilde(iter+1, i) = c_tilde(iter, i) + dt*  ( (c2 - 2*c_tilde(iter, i) + c_tilde(iter, i-1) )/(h^2) ...
                                                         + (c2 - c_tilde(iter, i-1))/(2*h*r_tilde(i))  ...
                                                         - beta*c_tilde(iter, i)*(sigma_h_tilde(iter, i) - sigma_h_tilde(iter, i-1))/(h*r_tilde(i))  ...
                                                         - beta*c_tilde(iter, i)*(sigma_h_tilde(iter, i) - 2*sigma_h_tilde(iter, i-1) + sigma_h_tilde(iter, i-2))/(h^2) ) ;
 
            u_tilde(iter+1, i) = 1;
    
            epsilon_rr =  u_tilde(iter, i)/r_tilde(i) ;
            epsilon_tt =  (u_tilde(iter, i)-u_tilde(iter, i-1))/(h)  ;

            
            sigma_h_tilde(iter+1, i) = (epsilon_tt + epsilon_rr - 2*a*(1+nu))/(3*(1-2*nu)) ;
            sigma_rr_tilde(iter+1, i) = (nu*epsilon_tt + (1-nu)*epsilon_rr - a*(1+nu))/((1+nu)*(1-2*nu)) ;
            sigma_tt_tilde(iter+1, i) = (nu*epsilon_rr + (1-nu)*epsilon_tt - a*(1+nu))/((1+nu)*(1-2*nu)) ;
        else
            c_tilde(iter+1, i) = c_tilde(iter, i) + dt*  ( (c_tilde(iter, i+1) - 2*c_tilde(iter, i) + c_tilde(iter, i-1) )/(h^2) ...
                                                         + (c_tilde(iter, i+1) - c_tilde(iter, i-1))/(2*h*r_tilde(i))  ...
                                                         - beta*c_tilde(iter, i)*(sigma_h_tilde(iter, i+1) - sigma_h_tilde(iter, i-1))/(2*h*r_tilde(i))  ...
                                                         - beta*c_tilde(iter, i)*(sigma_h_tilde(iter, i+1) - 2*sigma_h_tilde(iter, i) + sigma_h_tilde(iter, i-1))/(h^2) ) ;

           
            u_tilde(iter+1, i) = ( (u_tilde(iter, i+1) +  u_tilde(iter, i-1))/(h^2) +  (u_tilde(iter, i+1)- u_tilde(iter, i-1))/(2*h*r_tilde(i)) - alpha*(c_tilde(iter, i+1)-c_tilde(iter, i-1))/(2*h) )/l ;
    
            epsilon_rr =  u_tilde(iter, i)/r_tilde(i) ;
            epsilon_tt =  (u_tilde(iter, i+1)-u_tilde(iter, i-1))/(2*h)  ;

            
            sigma_h_tilde(iter+1, i) = (epsilon_tt + epsilon_rr - 2*a*(1+nu))/(3*(1-2*nu)) ;
            sigma_rr_tilde(iter+1, i) = (nu*epsilon_tt + (1-nu)*epsilon_rr - a*(1+nu))/((1+nu)*(1-2*nu)) ;
            sigma_tt_tilde(iter+1, i) = (nu*epsilon_rr + (1-nu)*epsilon_tt - a*(1+nu))/((1+nu)*(1-2*nu)) ;

        end
    end
end


plot(r_tilde, c_tilde(n_iterations, :), linewidth=2);
figure
plot(r_tilde, u_tilde(n_iterations, :), linewidth=2);
figure
plot(r_tilde, sigma_h_tilde(n_iterations, :), linewidth=2);
title('Hydrostatic stress', 'Color', 'Blue');
figure
plot(r_tilde, sigma_tt_tilde(n_iterations, :), linewidth=2) ;
title('Hoop Stress', 'Color', 'Blue')
figure
plot(r_tilde, sigma_tt_tilde(n_iterations, :), linewidth=2) ;
title('Radial Stress', 'Color', 'Blue');


