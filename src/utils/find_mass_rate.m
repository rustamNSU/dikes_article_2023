function find_mass_rate(newFrac, oldFrac, mesh, mQ, dt)
    newFrac.mass_rate(1:end) = 0;
    newFrac.flux_velocity(1:end) = 0;
    mNew = newFrac.get_mass();
    mOld = oldFrac.get_mass();
    
    %% q_i = 0 in fracture tips
    internalInds = newFrac.get_fracture_elements();
    internalInds(end) = [];
    for ind = internalInds
        newFrac.mass_rate(ind+1) = (mOld(ind) - mNew(ind) + mesh.dx(ind) * mQ(ind)) / dt + newFrac.mass_rate(ind);
        newFrac.flux_velocity(ind+1) = newFrac.mass_rate(ind+1) / (0.5 * (newFrac.rho(ind) + newFrac.rho(ind+1)));
%         if newFrac.mass_rate(ind+1) > 0
%             newFrac.flux_velocity(ind+1) = newFrac.mass_rate(ind+1) / newFrac.rho(ind);
%         else
%             newFrac.flux_velocity(ind+1) = newFrac.mass_rate(ind+1) / newFrac.rho(ind+1);
%         end
    end
end