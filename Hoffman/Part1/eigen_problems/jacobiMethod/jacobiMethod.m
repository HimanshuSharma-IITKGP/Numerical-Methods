

% A is a real symmetric matrix

function [A, V] = jacobiMethod(A, V)
   
    
    if(norm(diag(diag(A)) - A)<(1e-4))
        A = A;
        V = V;
        return ;
    end



    [p, q] = getMaxIndex(A);

    J = jacobiMatrix(A, p, q);
    V = V*J ;
    A = (J')*A*J ;

    jacobiMethod(A, V);
    
end