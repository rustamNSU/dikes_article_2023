function [Acc, Act, Atc, Att] = generate_compressibility_mobility_matrix(mesh, frac, settings)
    %% Fill global mobility matrix
    A = zeros(mesh.n, mesh.n);
    channel = frac.channel;
    tip     = frac.tip;
    fracElements = sort([channel, tip]);
    MIN_CONDUCTIVE_WIDTH = 1e-4;
    meanViscosity = mean(frac.mu(fracElements), "all");
    for i = 1:length(fracElements)
        element = fracElements(i);
        wl3 = max(frac.width(element-1), MIN_CONDUCTIVE_WIDTH)^3;
        wc3 = max(frac.width(element), MIN_CONDUCTIVE_WIDTH)^3;
        wr3 = max(frac.width(element+1), MIN_CONDUCTIVE_WIDTH)^3;
        
        rhol = frac.rho(element-1);
        rhoc = frac.rho(element);
        rhor = frac.rho(element+1);
        
        mul = frac.mu(element-1);
        muc = frac.mu(element);
        mur = frac.mu(element+1);
        if meanViscosity > 1e12
            mul = +Inf;
            muc = +Inf;
            mur = +Inf;
        end
        
        dxl = mesh.dx(element-1);
        dxc = mesh.dx(element);
        dxr = mesh.dx(element+1);

        coef_frac = 12;
        if settings.FRACTURE_TYPE == "PKN"
            coef_frac = pi^2;
        end
    
        
        %% Left tip element (boundary)
        if i==1
            mobilityR = (rhor*wr3*dxr/mur + rhoc*wc3*dxc/muc) / (dxr+dxc) / (coef_frac * dxc^2);
            A(element, element) = -mobilityR;
            A(element, element+1) = mobilityR;
            continue
        end
        
        
        %% Right tip element (boundary)
        if i==length(fracElements)
            mobilityL = (rhol*wl3*dxl/mul + rhoc*wc3*dxc/muc) / (dxl+dxc) / (coef_frac * dxc^2);
            A(element, element) = -mobilityL;
            A(element, element-1) = mobilityL;
            continue
        end
        
        
        %% Inner elements
        mobilityL = (rhol*wl3*dxl/mul + rhoc*wc3*dxc/muc) / (dxl+dxc) / (coef_frac * dxc^2);
        mobilityR = (rhor*wr3*dxr/mur + rhoc*wc3*dxc/muc) / (dxr+dxc) / (coef_frac * dxc^2);
        A(element, element-1) = mobilityL;
        A(element, element) = -(mobilityL + mobilityR);
        A(element, element+1) = mobilityR;
    end %Fill global mobility matrix
    
    Acc = A(channel, channel);
    Act = A(channel, tip);
    Atc = A(tip, channel);
    Att = A(tip, tip);
end

