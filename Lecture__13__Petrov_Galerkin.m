% Solution of the ODE d^2u/dx^2 + u - x^2 = 0 with the boundary condition
% u(0) = 0 and u'(1) = 1 using Petrov-Galerkin method

clear; clc;

syms x;
n = 2;
phi_0 = x;
phi_array = [(2*x - x^2), (x^2 - (2*x^3)/3)];
psi_array = create_psi(x, n);

B = createB(x, n, phi_array, psi_array) ;
f = createF(x, n, psi_array, phi_0);

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


function [B] = createB(x, n, phi_array, psi_array)
    B = zeros(n, n);

    for i=1:n
        for j =1:n
            B(i, j) = bilinearFunctional(phi_array, psi_array, i, j, x);
        end
    end
end

function [f] = createF(x, n, psi_array, phi_0)
    f = zeros(n, 1);

    for i=1:n
        f(i) = linearFunctional(psi_array, phi_0, i, x);
    end

end

% function [phi_array] = create_phi(x, n)
%     
%     phi_1 = x*(x-1)^2;
% 
%     phi_array(1) = phi_0;
% 
%     for i=2:n
%         phi_array(i) = x^(i-1)*phi_1;
%     end
% end


function [psi_array] = create_psi(x, n)
    
    psi_1 = x;

    psi_array(1) = psi_1;

    for i=2:n
        psi_array(i) = x*(i-1)*x;
    end
end

function [f_i] = linearFunctional(psi_array, phi_0, i, x)
    psi_i = psi_array(i);

    f_i = int((x^2-linearOperator(phi_0, x))*psi_i, 0, 1) ;
    
end

function [B_ij] = bilinearFunctional(phi_array, psi_array, i, j, x);
    psi_i = psi_array(i);
    phi_j = phi_array(j);

    B_ij = int(psi_i*linearOperator(phi_j, x), 0, 1);
end


function  [D] = linearOperator(y, x)
    D = diff(y, x, 2) + y ;
end

