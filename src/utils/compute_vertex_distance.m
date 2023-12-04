% Compute distance from asymptotics, Dontsov
function s = compute_vertex_distance(Ep, Kp, Cp, mup, w, sOld, dt)
    betam = 3.147345190264944;
    betamt = 2.533559440826570;
    
    sk = (Ep * w / Kp)^2;
    sm = nthroot(sOld^3 + (3.0 * Ep * dt * w^3)/(mup * betam^3), 3);
    smt = realmax;
    if (Cp > 1e-18)
        smt = nthroot(sOld^6 + (3.0 * Ep^2 * dt * w^8)/(2.0 * betamt^8 * mup^2 * Cp^2), 6);
    end
    s = min([sk sm smt]);
end

