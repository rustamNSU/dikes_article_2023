%% generate elasticity matrix
function C = generate_elasticity_matrix(reservoir, mesh, settings)
    C = zeros(mesh.n, mesh.n);
    if settings.FRACTURE_TYPE == "KGD"
        Ep = reservoir.Ep / (4.0*pi);
        for m = 1:mesh.n
            for k = 1:mesh.n
                C(m,k) = -Ep *...
                    (1.0/(mesh.xc(m) - mesh.rightBoundary(k)) -...
                    1.0/(mesh.xc(m) - mesh.leftBoundary(k)));
            end
        end
    else
        h = reservoir.h;
        Ep = reservoir.Ep * 2 / (pi * pi * h);
        cmk = zeros(1, mesh.n);
        for k = 1:mesh.n
            cmk(k) = -Ep *...
                (pkn_g(2 * (mesh.xc(1) - mesh.rightBoundary(k)) / h) -...
                pkn_g(2 * (mesh.xc(1) - mesh.leftBoundary(k)) / h));
        end
        for m = 1:mesh.n
            for k = 1:mesh.n
                % C(m,k) = -Ep *...
                %     (pkn_g(2 * (mesh.xc(m) - mesh.rightBoundary(k)) / h) -...
                %     pkn_g(2 * (mesh.xc(m) - mesh.leftBoundary(k)) / h));
                C(m,k) = cmk(abs(k-m) + 1);
            end
        end
    end
end


function g = pkn_g(s)
    a = ellipticE(1.0 / (1.0 + s.^2));
    g = sqrt(1.0 + s.^2) / s * a;
end