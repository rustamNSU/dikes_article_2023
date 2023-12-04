% Newton's method
% ********************************************************************
%  * Find root of equation func(x) = 0 using Newton's method 
%  * == INPUT == 
%  * func       -- function to solve
%  * x0         -- initial approximation of root
%  * tol        -- Newton's method tolerance
%  * itt_max    -- max number of iterations
%  * del_x      -- step for derivative calculation
%  * == OUTPUT ==
%  * @return x_Newton -- root of the equation func(x) = 0
%  ********************************************************************/

function x_Newton = NewtonMethod(func, x_0, del_x, itt_max, tol)

    x_old = x_0;
    x_Newton = x_old;    
    
    itt_Newton = 0;
    Res_Newton = 1;

    while ( itt_Newton < itt_max ) && ( Res_Newton > tol )
%         func(x_old + del_x) - func(x_old)
        derivative = ( func(x_old + del_x) - func(x_old) ) / del_x;
        x_Newton = x_old - func(x_old) / derivative;

        %w_init < 0 or complex
        ix_Newton0 = find((x_Newton<=0)|(abs(imag(x_Newton))>0));
        x_Newton(ix_Newton0) = 1.0e-8;   
        
        if ( x_Newton <= x_old )
            x_Newton = x_old;
        end
        
%         if ( x_Newton < 0 )
%            x_Newton = 1.0e-8;
%         end

        Res_Newton = abs((x_Newton - x_old)/x_old);

        x_old = x_Newton;
%         Res_Newton = norm(func(x_old));

        
        itt_Newton = itt_Newton + 1;

    end
    
    if ( itt_Newton >= itt_max )
        error('Newton_method didnt converge!')
    end
    
end