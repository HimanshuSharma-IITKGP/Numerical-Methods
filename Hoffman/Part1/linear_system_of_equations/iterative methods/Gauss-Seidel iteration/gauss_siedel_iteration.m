

function [X_new, iter] = gauss_siedel_iteration(A, b, X_0, conv)
    
    X = X_0;
    X_new = X;
    err = 1;
    iter = 0;
    while(err>conv)
        for i = 1:size(X_0, 1)
            S1 = 0;
            for j = 1:(i-1)
                S1 = S1 + A(i, j)*X_new(j);
            end
    
            S2 = 0;
            for j = i:size(X_0, 1)
                S2 = S2 + A(i, j)*X(j);
            end

            R = b(i) - S1 - S2;
    
            X_new(i) = X(i) + R/(A(i, i));
        end
    
        err = norm(X_new - X);
        X = X_new ;
        iter = iter + 1 ;
    end
end