

function [p2] = secantMethod(fn, p0, p1, conv)
    
    if(fn(p0) == 0)
        p2 = a;
        return 
    end 

    if(fn(p1) == 0)
        p2 = b;
        return 
    end 

    if(fn(p0)*fn(p1)>0)
        disp("enter some other range") ;
        return ;
    end
   


    err = 100;
    while(err>conv)
        q0 = fn(p0);
        q1 = fn(p1);

        p2 = p0 - q0*(p1-p0)/(q1-q0) ; 

        err = abs(fn(p2)/p0) ;

        if(fn(p2)*fn(p1)<0)
            p0 = p2;
        elseif(fn(p2)*fn(p0)<0)
            p1 = p2;
        end
    end
end