function [c] = bisectionMethod(fn, a, b, conv)

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
   

    err = 100;

    while(err>conv)

        c = (a+b)/2 ;
        
        if(fn(c)*fn(a)<0)
            b = c;
        elseif(fn(c)*fn(b)<0)
            a = c ;
        else
            break ;
        end

        err = abs(a-b)/abs(c) ;

    end
end