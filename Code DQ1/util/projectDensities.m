function xProj = projectDensities(x)
    global Hproj Heta Hbeta

    if Hproj
        xProj = (tanh(Hbeta*Heta) + tanh(Hbeta*(x - Heta)))/(tanh(Hbeta*Heta) + tanh(Hbeta*(1 - Heta)));
    else
        xProj = x;
    end
end