function w = calculate_asymptotic_width(Ep, Kp, Cp, mup, s, v)
    %% Define initial solution like vertex
    w0 = max(calculate_vertex_width(Ep, Kp, Cp, mup,s, v));
   
    
    %% Algorithm parameters
    ASYMPTOTIC_TOLERANCE = 1e-8;
    ASYMPTOTIC_MAX_ITERATIONS = 200;
    DELTA_X = 1e-8;
    f = @(w)implicit_asymptotic(Ep, Kp, Cp, mup, w, s, v);
    w = NewtonMethod(f, w0, DELTA_X, ASYMPTOTIC_MAX_ITERATIONS, ASYMPTOTIC_TOLERANCE); 
end

