% solution for d^2u/dx^2 = k(constant) for dirichlet boundary conditions
% u(x1) = u1 and u(x2) = u2

clear; clc;
syms x;

n = 3;
x1 = 0;
x2 = 1;
u1 = 0;
u2 = 0;
k = 5;

phi_0 = (u2-u1)*(x-x1)/(x2-x1) + u1 ; % phi_0(x1) = u1 and phi_0(x2) = u2



phi_array = create_phi(n, x, x1, x2);
B = createB(n, phi_array, x, x1, x2);
f = createF(n, phi_array, x1, x2, k);
    
c = pinv(B)*f;
hold on;
axis on;
sol = createSolution(c, phi_array, phi_0, n);
fplot(x, sol, linewidth=2);
uSol = matlabFunction(sol);



function [U] = createSolution(c, phi_array,phi_0 ,n)
    U = 0;
    for i=1:n
        U = U+c(i)*phi_array(i);
    end
    U = U + phi_0;
    U = simplify(U);

    disp(U);
end



function [B] = createB(n, phi_array, x, x1, x2)
    B = zeros(n, n);
    for i=1:n
        for j=1:n
            B(i, j) = bilinearFunctional(phi_array, i, j, x, x1, x2);
        end
    end
end


function [f] = createF(n, phi_array, x1, x2, k)
    
    f = zeros(n, 1) ;

    for i=1:n
        f(i) = linearFunctional(phi_array, i, x1, x2, k);
    end
end

function [phi_array] = create_phi(n, x, x1, x2)
    
    phi1 = (x2-x)*(x1-x); % phi_i(x1) = phi_i(x2) = 0

    phi_array(1) = phi1;
    for i=2:n
         phi_array(i) = (x^(i-1))*phi1;
    end

end

function [B_ij] = bilinearFunctional(phi_array, i, j, x, x1, x2)
    phi_i = phi_array(i);
    phi_j = phi_array(j);

    dDxi = diff(phi_i, x);
    dDxj = diff(phi_j, x);

    B_ij = int(dDxi*dDxj, x1, x2);
end


function [f_i] = linearFunctional(phi_array, i,x1, x2, k)
    f_i = -k*int(phi_array(i), x1, x2);
end