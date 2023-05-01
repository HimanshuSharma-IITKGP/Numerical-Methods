function [x_new, iter] = NewtonRaphsonMethodVisual(fx, x0, conv)
    syms x ;
    dfdx = matlabFunction(diff(fx(x), x), 'vars', {x});

    err = 100;
    iter = 0 ;
    
    figure;
    xlim([x0-2, x0+2]);
    fplot(fx,([x0-2, x0+2]));
    hold on;
    plot([x0-2, x0+2], [0 0]);


    while(err>conv)
        if(dfdx(x0)==0)
            disp('failed to converge')
            return 
        end

        x_new = x0 - fx(x0)/dfdx(x0);
        iter = iter + 1;
        err = abs((x_new - x0)/x0) ;

        hold on;
        plot([x0, x0], [0, fx(x0)], '-');
        plot([x0 x_new], [fx(x0), 0], '--');

        x0 = x_new ;
    end
end