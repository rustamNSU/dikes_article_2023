function sNew = compute_survey_distance(Ep, Kp, Cp, mup, w, sOld, dt)
    %% Define initial solution like vertex
    s0 = compute_vertex_distance(Ep, Kp, Cp, mup, w, sOld, dt);
    if s0 <= sOld
        sNew = sOld;
        return
    end
    
    
    %% Algorithm parameters
    ASYMPTOTIC_TOLERANCE = 1e-8;
    ASYMPTOTIC_MAX_ITERATIONS = 200;
    DELTA_X = 1e-8;
    f = @(s)asymptotic_s(Ep, Kp, Cp, mup, w, s, sOld, dt);
    sNew = NewtonMethod(f, s0, DELTA_X, ASYMPTOTIC_MAX_ITERATIONS, ASYMPTOTIC_TOLERANCE); 
end


%%
function residual = asymptotic_s(Ep, Kp, Cp, mup, w, s0, sOld, dt)
    MIN_VELOCITY = 1.e-10;
    V = max([(s0 - sOld) / dt MIN_VELOCITY]);
    residual = implicit_asymptotic(Ep, Kp, Cp, mup, w, s0, V);
end %asymptotic_s
