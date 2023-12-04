%% generate elasticity matrix
function C = generate_elasticity_matrix(Ep, mesh)
    C = zeros(mesh.n, mesh.n);
    for m = 1:mesh.n
        for k = 1:mesh.n
%             if abs(k-m) > 4
%                 continue;
%             end
            C(m,k) = -Ep/(4.0*pi) *...
                (1.0/(mesh.xc(m) - mesh.rightBoundary(k)) -...
                 1.0/(mesh.xc(m) - mesh.leftBoundary(k)));
        end
    end
end