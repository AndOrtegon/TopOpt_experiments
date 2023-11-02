function [c, vol, dc, dv] = analyze(xPhys,xProj)
global ft KE F H Hs edofMat freedofs iK jK Hproj Hbeta Heta nelx nely E0 Emin problem penal penal_type method U

if strcmp(problem,'Elasticity')
    % Penalize
    if strcmp(penal_type, 'SIMP')
        xPen = Emin + (E0-Emin)*xProj.^penal;
    elseif strcmp(penal_type, 'RAMP')
        xPen = Emin + (E0-Emin)*xProj./(1+penal*(1-xProj));
    elseif strcmp(penal_type, 'POL')
        xPen = Emin + (E0-Emin)*(xProj+2*xProj.^4)/3;
    end
    
    % Assemble matrix
    sK = reshape(KE(:)*xPen(:)',64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;

    U = zeros(2*(nely+1)*(nelx+1),1);
    
    % Solve by one of the methods
    if strcmp(method,'direct')
        U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    elseif strcmp(method,'iterative')
        tol = 1e-5 ;
        maxit = 1e4 ;

        ME.identifier = [];
        try
            L = ichol(K(freedofs,freedofs));
        catch ME
        end

        if (strcmp(ME.identifier,'MATLAB:ichol:Breakdown'))
            msg = ['ichol encountered nonpositive pivot, using no preconditioner.'];

            % you might consider tring different preconditioners (e.g. LU) in
            % the case ichol breaks down. We will default to no preconditioner:
            U(freedofs) = pcg(K(freedofs,freedofs),...
                F(freedofs),...
                tol,maxit, ...
                [] ... % preconditioner
                );
        else
            msg = [];
            U(freedofs) = pcg(K(freedofs,freedofs),...
                F(freedofs), ...
                tol, maxit, ...
                L,L.' ... % preconditioner
                );
        end
        disp(msg)
    end

    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    c  = sum(sum(xPen.*ce));

    % Find sensitivities
    if strcmp(penal_type, 'SIMP')
        dc = -penal*(E0-Emin)*xProj.^(penal-1).*ce;
    elseif strcmp(penal_type, 'RAMP')
        dc = -(E0-Emin)*(1+penal)*((1+penal*(1-xProj)).^-2).*ce;
    elseif strcmp(penal_type, 'POL')
        dc = -(E0-Emin)*(1+8*xProj.^3).*ce/3;
    end

    dv = ones(nely,nelx);
    %% HEAVISIDE PROJECTION / MODIFICATION OF SENSITIVITIES
    if Hproj
        dc(:) = (Hbeta*sech(Hbeta*(xPhys(:) - Heta)).^2)/(tanh(Hbeta*Heta) + tanh(Hbeta*(1 - Heta))).*dc(:);
        dv(:) = (Hbeta*sech(Hbeta*(xPhys(:) - Heta)).^2)/(tanh(Hbeta*Heta) + tanh(Hbeta*(1 - Heta))).*dv(:);
    end

elseif strcmp(problem,'Heat')
    % Penalize
    if strcmp(penal_type, 'SIMP')
        xPen = Emin + (E0-Emin)*xProj.^penal;
    elseif strcmp(penal_type, 'RAMP')
        xPen = Emin + (E0-Emin)*xProj./(1+penal*(1-xProj));
    elseif strcmp(penal_type, 'POL')
        xPen = Emin + (E0-Emin)*(xProj+2*xProj.^4)/3;
    end
    
    % Assemble matrix
    sK = reshape(KE(:)*xPen(:)',16*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;

    U = zeros((nely+1)*(nelx+1),1);
    
    % Solve by one of the methods
    if strcmp(method,'direct')
        U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    elseif strcmp(method,'iterative')
        tol = 1e-5 ;
        maxit = 1e4 ;

        ME.identifier = [];
        try
            L = ichol(K(freedofs,freedofs));
        catch ME
        end

        if (strcmp(ME.identifier,'MATLAB:ichol:Breakdown'))
            msg = ['ichol encountered nonpositive pivot, using no preconditioner.'];

            % you might consider tring different preconditioners (e.g. LU) in
            % the case ichol breaks down. We will default to no preconditioner:
            U(freedofs) = pcg(K(freedofs,freedofs),...
                F(freedofs),...
                tol,maxit, ...
                [] ... % preconditioner
                );
        else
            msg = [];
            U(freedofs) = pcg(K(freedofs,freedofs),...
                F(freedofs), ...
                tol, maxit, ...
                L,L.' ... % preconditioner
                );
        end
        disp(msg)
    end

    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    c = sum(sum(xPen.*ce));

    % Find sensitivities
    if strcmp(penal_type, 'SIMP')
        dc = -(E0-Emin)*penal*xProj.^(penal-1).*ce;
    elseif strcmp(penal_type, 'RAMP')
        dc = -(E0-Emin)*(1+penal)*((1+penal*(1-xProj)).^-2).*ce;
    elseif strcmp(penal_type, 'POL')
        dc = -(E0-Emin)*(1+8*xProj.^3).*ce/3;
    end

    dv = ones(nely,nelx);
end

%% FILTERING/MODIFICATION OF SENSITIVITIES
if ft == 1
    dc(:) = H*(xPhys(:).*dc(:))./Hs./max(1e-3,xPhys(:));
elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
end

vol = sum(xProj(:)) ;
end