function [sError, nTipElements] = update_survey_distance(startTime, endTime,...
        mesh, reservoir, settings, newFrac, oldFrac)
    sOld = oldFrac.distanceToSurvey;
    sNew = newFrac.distanceToSurvey;
    sPrev = sNew;
    w = newFrac.width(newFrac.survey);
    if settings.FRACTURE_TYPE == "PKN"
        w = 4 / pi .* w;
    end
    mu = newFrac.mu(newFrac.survey);
    
    Ep = reservoir.Ep;
    Kp = reservoir.Kp;
    Cp = reservoir.Cp;
    dt = endTime - startTime;
    
    %% calculate sNew
    for i = 1:length(sNew)
        sNew(i) = compute_survey_distance(Ep, Kp, Cp, 12 * mu(i), w(i), sOld(i), dt);
    end
    sNew = (1-settings.FRONT_RELAXATION)*sPrev + settings.FRONT_RELAXATION*sNew;
    
    if settings.FIXED_BOTTOM_TIP
        sNew(1) = sPrev(1);
    end
    newFrac.set_survey_distance(sNew);
    
    %%
    sError = norm(sNew - sPrev, Inf) / norm(sNew, Inf);
    nTipElements = max([length(newFrac.leftTip), length(newFrac.rightTip)]);
end

