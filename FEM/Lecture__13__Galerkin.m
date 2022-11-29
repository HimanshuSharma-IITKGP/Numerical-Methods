% solution of the ODE d^2u/dx^2 + u - x^2 = 0 with boundary conditions
% u(0) = 0 and u'(1) = 1 using Galerkin's weighted residual method

clear; clc;
syms x;
fx = x^2 ;

n = 2;
phi_0 = x;

phi_array = [(2*x - x^2), (x^2 - (2*x^3)/3)];
B = createB(n, phi_array, x);
f = createF(n, phi_array, x, phi_0, fx);

c = pinv(B)*f;
hold all;
axis on;
sol = createSolution(c, phi_array, phi_0, n);
fplot(x, sol, linewidth=2);
uSol = matlabFunction(sol);
exact_sol = ((2*cos(1-x) - sin(x))/cos(1)) + x^2 - 2;
fplot(x, exact_sol);
xlim([0, 1]);
ylim([0, 1.5]);


function [U] = createSolution(c, phi_array,phi_0 ,n)
    U = 0;
    for i=1:n
        U = U+c(i)*phi_array(i);
    end
    U = U + phi_0;
    disp(U);
    U = simplify(U);

    disp(U);
end

function [B] = createB(n, phi_array, x)
    B = zeros(n, n);
    for i=1:n
        for j=1:n
            B(i, j) = bilinearFunctional(phi_array, i, j, x);
        end
    end
end

function [f] = createF(n, phi_array, x, phi_0, fx)
    
    f = zeros(n, 1) ;

    for i=1:n
        f(i) = linearFunctional(phi_array, phi_0, i, x, fx);
    end
end


% function [phi_array] = create_phi(n, x)
% 
%     phi1 = x*(x-1)^2;
% 
% 
%     phi_array(1) = phi1;
%     for i=2:n
%          phi_array(i) = (x^(i-1))*phi1;
%     end
%    
% end


function [B_ij] = bilinearFunctional(phi_array, i, j, x)
    phi_i = phi_array(i);
    phi_j = phi_array(j);


    B_ij = int(linearOperator(phi_j, x)*phi_i, 0, 1);
%     disp(B_ij);
end

function [f_i] = linearFunctional(phi_array, phi_0, i, x, fx)
    phi_i = phi_array(i);

    f_i = int(((fx - linearOperator(phi_0, x))*phi_i), 0, 1);
%     disp(f_i);
end

function [D] = linearOperator(y, x)
    D = diff(y, x, 2) + y ;
end