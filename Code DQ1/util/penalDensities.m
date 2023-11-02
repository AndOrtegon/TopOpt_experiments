function xPen = penalDensities(xProj)
    global penal_type Emin E0 penal

    if strcmp(penal_type, 'SIMP') 
        xPen = Emin + (E0-Emin)*xProj.^penal;
    elseif strcmp(penal_type, 'RAMP') 
        xPen = Emin + (E0-Emin)*xProj./(1+penal*(1-xProj));
    elseif strcmp(penal_type, 'POL') 
        xPen = Emin + (E0-Emin)*(xProj+2*xProj.^4)/3;
    end
end