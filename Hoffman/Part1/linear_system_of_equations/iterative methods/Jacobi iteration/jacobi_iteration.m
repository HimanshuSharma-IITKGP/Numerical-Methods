


function [X_new, iter] = jacobi_iteration(A, b, X_0, conv)
    
    X = X_0;
    X_new = X ;
    disp(conv)
    err = 1;
    iter = 0;
    while(err>conv)

        for i = 1:size(X_0, 1)
            S = 0;
            for j = 1:size(X_0, 1)
                S = S + A(i, j)*X(j) ;
            end
            R = b(i) - S ;
            X_new(i) = X(i) + R/(A(i, i));
        end
        
        err = norm(X_new - X);
        X = X_new ;
        iter = iter + 1;
    end
    
end