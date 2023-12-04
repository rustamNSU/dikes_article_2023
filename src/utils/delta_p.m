function result = delta_p(Kh, Ch, p)
    %% Validate input
    Kh = min([max([0.0 Kh]) 1.0]);
    Ch = max([0.0 Ch]);
    largeLeakoff = 1000.0;
    
    betam = 3.147345190264944;
    betamt = 2.533559440826570;
    b0 = 0.991179982305673;
    
    if Ch < largeLeakoff
        delta = betam^3/3*solution_f(Kh, b0*Ch, betam^3/3)*(1+b0*Ch);
    else
        delta = betam^3/3*expanded_solution_f(Kh, b0*Ch, betam^3/3)*(1+b0*Ch);
    end
    
    if Ch < largeLeakoff
        result = delta * (1 - p +...
            p*solution_f(Kh, b0*Ch, betam^3/3)*(betam^3+betamt^4*Ch));
    else
        result = delta * (1 - p +...
            p*expanded_solution_f(Kh, b0*Ch, betam^3/3)*(betam^3+betamt^4*Ch));
    end
end


%% Solution of approximate diff. equation, see [1] (A 1), [2] (3.11)
function f = solution_f(Kh, Ch, C1)
    f = (1-Kh^3-1.5*Ch*(1-Kh^2)+3*Ch^2*(1-Kh)-6*Ch^3*atanh((1-Kh)/(2*Ch+1+Kh)))/(3*C1);
end %solution_f


%% Large leakoff solution_f
function f = expanded_solution_f(Kh, Ch, C1)
    f = ((1-Kh^4)/(4*Ch)-(1-Kh^5)/(5*Ch^2)+(1-Kh^6)/(6*Ch^3))/C1;
end %expanded_solution_f