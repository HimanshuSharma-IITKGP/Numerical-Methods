function [x_new, iter] = NewtonRaphsonMethod(fx, x0, conv)
    syms x ;
    dfdx = matlabFunction(diff(fx(x), x), 'vars', {x});

    err = 100;
    iter = 0 ;
    
    while(err>conv)
        if(dfdx(x0)==0)
            disp('failed to converge')
            return 
        end

        x_new = x0 - fx(x0)/dfdx(x0);
        iter = iter + 1;
        err = abs((x_new - x0)/x0) ;

        x0 = x_new ;
    end
end