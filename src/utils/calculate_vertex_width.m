% calculate vertex solutions
% for different fracture tip growth regimes
function [wk, wm, wmt] = calculate_vertex_width(Ep, Kp, Cp, mup,s, v)
    betam = 3.147345190264944;
    betamt = 2.533559440826570;
    
    wk = Kp / Ep * sqrt(s);
    wmt = betamt * nthroot(v * s^5 * (2 * mup * Cp / Ep)^2, 8);
    wm = betam * nthroot(s^2 * mup * v / Ep, 3);
end