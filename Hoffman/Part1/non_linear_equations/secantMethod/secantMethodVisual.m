

function [p2, iter] = secantMethodVisual(fn, p0, p1, conv)
    
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
   

    figure;
    xlim([p0-0.5, p1+0.5]);
    fplot(fn,[p0-0.5, p1+0.5]);
    hold on;
    plot([p0-0.5, p1+0.5], [0 0]);


    err = 100;
    iter = 0;
    while(err>conv)
        q0 = fn(p0);
        q1 = fn(p1);
        
        hold on;
        plot([p0 p0],[0 q0], '-')
        plot([p1 p1],[0 q1], '-')
        plot([p0, p1], [q0, q1], '--');

        p2 = p0 - q0*(p1-p0)/(q1-q0) ; 



        err = abs(fn(p2)/p0) ;

        if(fn(p2)*fn(p1)<0)
            p0 = p2;
        elseif(fn(p2)*fn(p0)<0)
            p1 = p2;
        end

        iter = iter + 1;
    end
end