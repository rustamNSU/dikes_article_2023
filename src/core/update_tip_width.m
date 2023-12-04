function update_tip_width(startTime, endTime,...
        mesh, reservoir, settings, newFrac, oldFrac) 
    Ep = reservoir.Ep;
    Kp = reservoir.Kp;
    Cp = reservoir.Cp;
    dt = endTime - startTime;
    
    sOld = oldFrac.distanceToSurvey;
    sNew = newFrac.distanceToSurvey;
    w = newFrac.width(newFrac.survey);
    mu = newFrac.mu(newFrac.survey);
    V = (sNew - sOld) / dt;
    newFrac.front_velocity = V;
    
    %% Find distance to front for all tip elements
    leftTip = newFrac.leftTip;
    rightTip = newFrac.rightTip;
    leftDistance = mesh.rightBoundary(leftTip) - newFrac.front(1);
    rightDistance = newFrac.front(2) - mesh.leftBoundary(rightTip);
    
    
    %% Calculate zero moments    
    leftZeroMoments = arrayfun(@(s) calculate_moments(Ep, Kp, Cp, 12 * mu(1), s, V(1), w(1)),...
        leftDistance);
    rightZeroMoments = arrayfun(@(s) calculate_moments(Ep, Kp, Cp, 12 * mu(2), s, V(2), w(2)),...
        rightDistance);
    
    %% Average tip width by volume
    leftAverageWidth = (leftZeroMoments - [0 leftZeroMoments(1:end-1)]) ./...
        newFrac.mesh.dx(leftTip);
    rightAverageWidth = (rightZeroMoments - [rightZeroMoments(2:end) 0]) ./...
        newFrac.mesh.dx(rightTip);
    
    %% Filter: cut width if tip width greater than survey width
%     leftAverageWidth(leftAverageWidth > w(1)) = max([0.2*w(1) settings.MIN_WIDTH]);
%     rightAverageWidth(rightAverageWidth > w(2)) = max([0.2*w(2) settings.MIN_WIDTH]);
    leftAverageWidth(leftAverageWidth > w(1)) = w(1);
    rightAverageWidth(rightAverageWidth > w(2)) = w(2);
    
    leftAverageWidth(leftAverageWidth < settings.MIN_WIDTH) = settings.MIN_WIDTH;
    rightAverageWidth(rightAverageWidth < settings.MIN_WIDTH) = settings.MIN_WIDTH;
    
    %% Update width in newFrac tip elements
    newFrac.width(leftTip) = leftAverageWidth;
    newFrac.width(rightTip) = rightAverageWidth;

    newFrac.update_tip_parameters();
end


%%
function V0 = calculate_moments(Ep, Kp, Cp, mup, s, v, w)
    ZERO_VELOCITY = 1e-12;
    if v < ZERO_VELOCITY
        Kp = w * Ep / sqrt(s);
        wk = Kp * sqrt(s) / Ep;
        V0 = 2 * wk * s / 3;
        return
    end
    
    wa = calculate_asymptotic_width(Ep, Kp, Cp, mup, s, v);
    Kh = Kp * sqrt(s) / (Ep * wa);
    Ch = 2 * sqrt(s) * Cp / (sqrt(v) * wa);
    p0 = 0.377;
    delta_p0 = delta_p(Kh, Ch, p0);
    
    V0 = 2 * wa * s / (3 + delta_p0);
end %calculate_moments

