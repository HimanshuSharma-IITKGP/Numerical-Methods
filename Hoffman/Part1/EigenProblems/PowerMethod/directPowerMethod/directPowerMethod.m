A = [8 -2 -2; -2 4 -2; -2 -2 13];
[vec, eig, iter] = directPowerMethodFunction(A, [1 1 1]', 1e-6)

function [v_new, lamda, iter] = directPowerMethodFunction(A, v0, conv)
   
    err = 1;
    v_new = zeros(size(v0));

    lamda = 0 ;

    iter = 0;
    while(err>conv)
        y = A*v0;
        lamda = max(abs(y));
        if(lamda==0)
            v_new = -1;
            break ;
        end
        v_new = y/lamda ;
        
        err = norm(v_new-v0)/norm(v0);
        
        v0 = v_new ;
        iter = iter + 1;
    end
end