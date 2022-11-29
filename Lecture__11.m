% solution of the equation d/dx[(1+x)du/dx] = 0 using Rayleigh-Ritz method
% with boundary condition u(0) = 0 and u(1) = 1;

clc; clear;
syms x;
phi_0 = x;

n = 2;

phi_array = create_phi(n, x);
B = createB(n, phi_array, x);
f = createF(n, phi_array, x, phi_0);

c = pinv(B)*f;
axis on;
sol = createSolution(c, phi_array, phi_0, n);
fplot(x, sol, linewidth=2);
uSol = matlabFunction(sol);

% uExact = log(1+x)/log(2)



function [U] = createSolution(c, phi_array,phi_0 ,n)
    U = 0;
    for i=1:n
        U = U+c(i)*phi_array(i);
    end
    U = U + phi_0;
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

function [f] = createF(n, phi_array, x, phi_0)
    
    f = zeros(n, 1) ;

    for i=1:n
        f(i) = linearFunctional(phi_array, phi_0, i, x);
    end
end


function [phi_array] = create_phi(n, x)

    
    phi1 = x*(1-x);


    phi_array(1) = phi1;
    for i=2:n
         phi_array(i) = (x^(i-1))*phi1;
    end
   
%     disp(phi_array)
end

function [B_ij] = bilinearFunctional(phi_array, i, j, x)
    phi_i = phi_array(i);
    phi_j = phi_array(j);

    dDxi = diff(phi_i, x);
    dDxj = diff(phi_j, x);

    B_ij = int(dDxi*(1+x)*dDxj, 0, 1);
%     disp(B_ij);
end


function [f_i] = linearFunctional(phi_array,phi_0, i, x)
    phi_i = phi_array(i);

    dDxi = diff(phi_i, x);
    dDx0 = diff(phi_0, x);

    f_i = -int(dDxi*(1+x)*dDx0, 0, 1);
%     disp(f_i);
end