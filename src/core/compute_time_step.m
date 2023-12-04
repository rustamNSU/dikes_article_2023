function compute_time_step(startTime, endTime, mesh, reservoir, schedule, settings, state, newFrac, oldFrac, level) 
    tab = "";
    if level ~= 0
        tab = "  ";
    end
    fprintf(tab + "========================================\n");
    fprintf(tab + 'Time step: %5.1f s -> %5.1f s\n', startTime, endTime);
    fprintf(tab + '-- dt = %5.1f s\n', endTime - startTime);
    %% Update survey and tip elements
    newFrac.update_survey_elements(); 
    oldFrac.deepCopy(newFrac);
    newFrac.time = endTime;
    
    %% Initialize time step algorithm parameters
    nIter = 0;
    MError = 1.0;
    sError = 1.0;
    nTipElements = 1;
    statusLubrication.isConvergence = false;
    dt = endTime - startTime;
    if dt < 1e-3
        a=10;
    end
    if dt < 1e-4
        error('Stop calculation because very small timestep');
    end
    dMInjected = schedule.injected_mass(startTime, endTime);
    
    %% Compute newFrac
    while (MError > settings.MASS_TOLERANCE ||...
                sError > settings.FRONT_TOLERANCE ||...
                nIter < 2) &&...
           nIter < settings.MAX_NONLINEAR_ITERATION

        meanViscosity = mean(newFrac.mu(newFrac.channel), "all");
        %% Front update
        if meanViscosity <= 1e10
            [sError, nTipElements] = update_survey_distance(startTime, endTime, mesh, reservoir, settings, newFrac, oldFrac);
            if nTipElements > settings.MAX_TIP_ELEMENTS
                fprintf(2, 'Limit maximum of new tip elements\n');
                break
            end
            
            % if (nTipElements > 1 && endTime > 101)
            %     a = 10;
            % end
            %% Tip update
            update_tip_width(startTime, endTime, mesh, reservoir, settings, newFrac, oldFrac);
        else
            sError = 0.0;
        end
        
        %% Lubrication
        statusLubrication = solve_compressibility_lubrication(startTime, endTime, mesh, reservoir, schedule, settings, state, newFrac, oldFrac);
        if not(statusLubrication.isConvergence)
            fprintf(2, 'Lubrication doesn''t converge\n');
            fprintf(2, 'Lubrication iterations = %d\n', statusLubrication.iterations);
            break
        end
        
        %% Check convergence
        dM = newFrac.get_fracture_fluid_mass() - oldFrac.get_fracture_fluid_mass();
        MError = abs(dM - dMInjected) / max([newFrac.get_fracture_fluid_mass(), dMInjected]);
        nIter = nIter + 1;
    end
    
    
    %% Display status
    fprintf(tab + '-- iteration = %d\n', nIter);
    fprintf(tab + '-- Merr = %5.10f\n', MError);
    fprintf(tab + '-- sError = %5.10f\n', sError);
    fprintf(tab + '-- M = %5.10f\n', newFrac.get_fracture_fluid_mass());
    
    
    %% Adaptive time step
    isDivideTimeStep = nIter >= settings.MAX_NONLINEAR_ITERATION ||...
        not(statusLubrication.isConvergence) ||...
        nTipElements > settings.MAX_TIP_ELEMENTS;
    if isDivideTimeStep
        fprintf(2, 'Time step: %5.6f -> %5.6f didn''t converge!\n', startTime, endTime);
        middleTime = 0.5 * (startTime + endTime);
        newFrac.deepCopy(oldFrac);
        compute_time_step(startTime, middleTime, mesh, reservoir, schedule, settings, state, newFrac, oldFrac, level+1);
        compute_time_step(middleTime, endTime, mesh, reservoir, schedule, settings, state, newFrac, oldFrac, level+1);
        return
    end %adaptive time step

    fprintf(tab + 'Total mass = %0.5f\n', newFrac.get_fracture_fluid_mass());
    fprintf(tab + "========================================\n");
end %compute_time_step

