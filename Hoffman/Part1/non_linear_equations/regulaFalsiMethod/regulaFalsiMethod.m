function [w, iter] = regulaFalsiMethod(fn, a, b, conv)

    if(fn(a) == 0)
        w = a;
        return 
    end 

    if(fn(b) == 0)
        w = b;
        return 
    end 

    if(fn(a)*fn(b)>0)
        disp("enter some other range") ;
        return ;
    end
   

    err = 100;

    iter = 0;
    while(err>conv)
        % disp(iter)
        w = (fn(b)*a-fn(a)*b)/(fn(b) - fn(a)) ;

        
        if(fn(w)*fn(a)<0)
            b = w;
        elseif(fn(w)*fn(b)<0)
            a = w ;
        else
            break ;
        end
        
        iter = iter + 1;
        err = abs(fn(w))/abs(a); % residual

    end
end