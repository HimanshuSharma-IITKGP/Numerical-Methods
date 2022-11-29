% solution of d^2u/dx^2 = k(constant) for neumann boundary conditions
% u'(x1) = p1 and u(x2) = u2 using Galerkin's Method

clear; clc;
syms x;

n = 1;
k = 5;
fx = k;
x1 = 0;
x2 = 1;
p1 = 0;
u2 = 1;

phi_0 = p1*x + (u2 - p1*x2) ;


phi_array = create_phi(n, x, x1, x2);

B = createB(n, phi_array, x, x1, x2);
f = createF(n, phi_array, x, phi_0, fx, x1, x2);

c = pinv(B)*f;
hold all;
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
    disp(U);
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

function [f] = createF(n, phi_array, x, phi_0, fx, x1, x2)
    
    f = zeros(n, 1) ;

    for i=1:n
        f(i) = linearFunctional(phi_array, phi_0, i, fx, x, x1, x2);
    end
end


function [phi_array] = create_phi(n, x, x1, x2)

    phi1 = x^2 - 2*x1*x + 2*x1*x2 - x2^2 ;


    phi_array(1) = phi1;
    for i=2:n
         phi_array(i) = ((x-x1)^(i-1))*phi1;
    end
   
end

function [B_ij] = bilinearFunctional(phi_array, i, j, x, x1, x2)
    phi_i = phi_array(i);
    phi_j = phi_array(j);


    B_ij = int(linearOperator(phi_j, x)*phi_i, x1, x2);
%     disp(B_ij);
end


function [f_i] = linearFunctional(phi_array, phi_0, i,fx, x, x1, x2)
    phi_i = phi_array(i);

    f_i = int(((fx - linearOperator(phi_0, x))*phi_i), x1, x2);
end

function [D] = linearOperator(y, x)
    D = diff(y, x, 2) ;
end

