function residual = implicit_asymptotic(Ep, Kp, Cp, mup, w, s0, V)
    s_hat = s0^2 * V * mup / (Ep * w^3); 
    K_hat = Kp * sqrt(s0) / (Ep * w);
    C_hat = 2 * sqrt(s0) * Cp / (w * sqrt(V));
    g_delta = fcn_g_del(K_hat, C_hat);
    
    residual = s_hat - g_delta;
end

