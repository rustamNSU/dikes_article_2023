function [Acc, Act, Atc, Att] = generate_mobility_matrix(...
        mesh, width, channel, tip, mup)
    %% Fill global mobility matrix
    A = zeros(mesh.n, mesh.n);
    fracElements = sort([channel, tip]);
    for i = 1:length(fracElements)
        element = fracElements(i);
        wl3 = width(element-1)^3;
        wc3 = width(element)^3;
        wr3 = width(element+1)^3;
        dxl = mesh.dx(element-1);
        dxc = mesh.dx(element);
        dxr = mesh.dx(element+1);
        
        mobilityL = (wl3*dxl + wc3*dxc) / (dxl+dxc) / (mup * dxc^2);
        mobilityR = (wr3*dxr + wc3*dxc) / (dxr+dxc) / (mup * dxc^2);
        
        %% Left tip element (boundary)
        if i==1
            A(element, element) = -mobilityR;
            A(element, element+1) = mobilityR;
            continue
        end
        
        %% Right tip element (boundary)
        if i==length(fracElements)
            A(element, element) = -mobilityL;
            A(element, element-1) = mobilityL;
            continue
        end
        
        %% Inner elements
        A(element, element-1) = mobilityL;
        A(element, element) = -(mobilityL + mobilityR);
        A(element, element+1) = mobilityR;
    end %Fill global mobility matrix
    
    Acc = A(channel, channel);
    Act = A(channel, tip);
    Atc = A(tip, channel);
    Att = A(tip, tip);
end

