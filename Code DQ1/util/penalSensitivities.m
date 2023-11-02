function dc = penalSensitivities(xProj,ce)
    global penal_type Emin E0 penal

    if strcmp(penal_type, 'SIMP')
        dc = -penal*(E0-Emin)*xProj.^(penal-1).*ce;
    elseif strcmp(penal_type, 'RAMP')
        dc = -(E0-Emin)*(1+penal)*((1+penal*(1-xProj)).^-2).*ce;
    elseif strcmp(penal_type, 'POL')
        dc = -(E0-Emin)*(1+8*xProj.^3).*ce/3;
    end
end