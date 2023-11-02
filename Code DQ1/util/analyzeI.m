function [c, vol, dc, dv] = analyzeI(x)
global ft Ke F H Hs edofMat freedofs iK jK nelx nely E0 Emin problem penal penal_type He

xPen = penalDensities(x) ;
if strcmp(problem,'Compliance')
elseif strcmp(problem,'Heat')
    xPenRes = kron(reshape(xPen,4,nely*nelx),ones(16,1)) ;
    sK = reshape(Ke(:).*xPenRes,4*16*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;

    U = zeros((nely+1)*(nelx+1),1);
    U(freedofs) = K(freedofs,freedofs)\F(freedofs) ;

    % Compliance
    cb = zeros(4,nely,nelx) ;
    cb(1,:,:) = reshape(sum((U(edofMat)*Ke(:,:,1)).*U(edofMat),2),nely,nelx);
    cb(2,:,:) = reshape(sum((U(edofMat)*Ke(:,:,2)).*U(edofMat),2),nely,nelx);
    cb(3,:,:) = reshape(sum((U(edofMat)*Ke(:,:,3)).*U(edofMat),2),nely,nelx);
    cb(4,:,:) = reshape(sum((U(edofMat)*Ke(:,:,4)).*U(edofMat),2),nely,nelx);

    c = sum(sum(sum(xPen.*cb)));

    xProj = x ;
    % Sensitivities
    if strcmp(penal_type, 'SIMP')
        dc = -(E0-Emin)*penal*xProj.^(penal-1).*cb;
    elseif strcmp(penal_type, 'RAMP')
        dc = -(E0-Emin)*(1+penal)*((1+penal*(1-xProj)).^-2).*cb;
    end

    dv = ones(4,nely,nelx);
end

%% FILTERING/MODIFICATION OF SENSITIVITIES
if ft == 2
    dc = He*reshape(ord2pos(dc),[],1);
    dv = He*reshape(ord2pos(dv),[],1);

    dc = pos2ord(reshape(dc,2*nely,2*nelx)) ;
    dv = pos2ord(reshape(dv,2*nely,2*nelx)) ;

    dc = dc(:) ;
    dv = dv(:) ;
end

vol = sum(x(:)) ;
end
