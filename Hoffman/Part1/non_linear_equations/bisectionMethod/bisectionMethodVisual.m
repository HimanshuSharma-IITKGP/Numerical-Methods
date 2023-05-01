function [c, iter] = bisectionMethodVisual(fn, a, b, conv, max_iter)

    figure;
    xlim([a-0.5, b+0.5]);
    fplot(fn,[a-0.5, b+0.5]);
    hold on;
    plot([a-0.5, b+0.5], [0 0]);
    if(fn(a) == 0)
        c = a;
        return 
    end 

    if(fn(b) == 0)
        c = b;
        return 
    end 

    if(fn(a)*fn(b)>0)
        disp("enter some other range") ;
        return ;
    end
    
    hold on;
    plot([a a], [0, fn(a)]);
    plot([b b], [0, fn(b)]);


    err = 100;
    iter = 1; 
    while(err>conv)

        if(iter == max_iter)
            disp("Failed to converge") ;
            return ;
        end
        c = (a+b)/2 ;

        hold on;
        plot([c c], [0 fn(c)]);


        if(fn(c)*fn(a)<0)
            b = c;
        elseif(fn(c)*fn(b)<0)
            a = c ;
        else
            break ;
        end

        iter = iter + 1;
        err = abs(a-b)/abs(c) ;

    end
end