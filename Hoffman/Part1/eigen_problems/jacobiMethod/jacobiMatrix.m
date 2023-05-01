

function [J] = jacobiMatrix(A, p, q)
    
    
    if(A(q, q) == A(p, p))
        c = 1/sqrt(2);
        s = 1/sqrt(2);
    elseif(A(p ,q) == 0)
        s = 0;
        c = 1;
    else 
        phi = (A(q, q) - A(p, p))/(2*A(p ,q)) ;
        t = sign(phi)/(abs(phi) + sqrt(phi^2 + 1));
        c = 1/sqrt(1 + t^2);
        s = c*t ;
    end

    J = eye(size(A));
    J(p, p) = c;
    J(p, q) = s;
    J(q, p) = -s;
    J(q, q) = c;
end