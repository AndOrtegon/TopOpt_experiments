function x = disaggregated(step_p) 

addpath(genpath('util'))

%% Declaration of global variable
global alldofs nelx nely volfrac penal He rmin ft E0 Emin Ke F H Hs edofMat freedofs iK jK ml penal_type Hproj Hbeta Heta problem epsilon lambda movp
global method plotting printing

movp = false ;
method = 'direct' ;

plotting = true ;
printing = true ;

% % Heat problem
problem = 'Heat' ;

nelx    = 300;
nely    = 300;
volfrac = 0.4;
penal   = 3;
rmin    = 2;
ft      = 2;
cont    = true ;
penal_type = 'SIMP';

% Move limit
ml      = 0.2;   %1.0 means essentially no move limits

% Heaviside projection
Hproj   = false;
Hbeta   = 8;
Heta    = 0.5;

%% Material properties
E0 = 1 ;
Emin = 1e-3 ;

Ke = zeros(4,4,4) ;

Ke(:,:,1) = [ 3 -1 -1 -1 ;
             -1  2  0 -1 ;
             -1  0  1  0 ;
             -1 -1  0  2 ]/12 ;

Ke(:,:,2) = [ 2 -1 -1  0 ;
             -1  3 -1 -1 ;
             -1 -1  2  0 ;
              0 -1  0  1 ]/12 ;

Ke(:,:,3) = [ 1  0 -1  0 ;
              0  2 -1 -1 ;
             -1 -1  3 -1 ;
              0 -1 -1  2 ]/12 ;

Ke(:,:,4) = [ 2  0 -1 -1 ;
              0  1  0 -1 ;
             -1  0  2 -1 ;
             -1 -1 -1  3 ]/12 ;

nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,4)+repmat([0 nely+[1 0] -1],nelx*nely,1);

iK      = reshape(kron(edofMat,ones(16,1))',[],1) ;
jK      = reshape(kron(edofMat,ones(4,4))' ,[],1);


F         = sparse((nely+1)*(nelx+1),1);
fixeddofs = nely/2+1-(nely/20):2:nely/2+1+(nely/20);
alldofs   = 1:(nely+1)*(nelx+1);
freedofs  = setdiff(alldofs,fixeddofs);
F(:,1)    = 1/(length(freedofs));

%% Filter
posX = kron(ones(2*nely,nelx),[1/3,2/3]) + kron(0:nelx-1,ones(2*nely,2)) ;
posY = kron(ones(nely,2*nelx),[1/3;2/3]) + kron((0:nely-1)',ones(2,2*nelx)) ;

iH = ones(nelx*nely*6,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k  = 0;
for i1 = 1:2*nelx % Iterar en todos los x numeración 2-indice
    for j1 = 1:2*nely % Iterar en todos los y numeración 2-indice
        e1 = (i1-1)*(2*nely)+j1; % Posición numeración 1-indice
        for i2 = max(1,i1-5):min(2*nelx,i1+5) % Iterar entre vecinos en x
            for j2 = max(1,j1-5):min(2*nely,j1+5) % Iterar entre vecinos en y
                if sqrt((posX(i1,j1)-posX(i2,j2))^2+(posY(i1,j1)-posY(i2,j2))^2) < rmin % Calcular distancia y determinar si menor a 0.7
                    e2 = (i2-1)*(2*nely)+j2; % Posición en numeración 1-indice de elemento vecino
                    k = k+1; % Contar elemento vecino
                    iH(k) = e1; % Elemento central en filas
                    jH(k) = e2; % Elemento vecino en columnas
                    sH(k) = max(0,rmin-sqrt((posX(i1,j1)-posX(i2,j2))^2+(posY(i1,j1)-posY(i2,j2))^2)); % Valor de la relación
                end
            end
        end
    end
end
H = sparse(iH,jH,sH,4*nely*nelx,4*nely*nelx);
Hs = sum(H,2);
He = spdiags(1./Hs,0,4*nely*nelx,4*nely*nelx)*H ;

%%
x     = volfrac*ones(4,nely,nelx) ;
xDis  = filterDensities(ord2pos(x)) ;
xPhys = pos2ord(xDis) ;

%% Optimization
hist=0;

% step_p = 0.001;
target_p = penal;

if strcmp(penal_type, 'SIMP')
    init_p = 1.0;
elseif strcmp(penal_type, 'RAMP')
    init_p = 0.0;
end

% INITIALIZE MMA PARAMETERS
change = 1;
outeriter = 0;
maxoutit  = 1300;

%     hist = zeros(maxoutit-1,nely*nelx) ;

kkttol  = 1e-4;
obj_tol = 1e-5;

TOLF = 1e-5;
% TOLF = 1e-7;
ndv = size(x(:),1);
ncons = 1;

a0 = 1;
ai = zeros(ncons,1);
ci = 1000*ones(ncons,1);
d = 0;

low = ones(ndv,1);
upp = ones(ndv,1);

xold1=x(:);
xold2=x(:);
levels = linspace(0,1,20);
% START ITERATION
%%%% The iterations start:
kktnorm = kkttol+10;
loop = 0;
old_obj = 0;
rel_obj_change = 10*obj_tol;
tic
while loop < maxoutit
    %% CONTINUATION
    if cont
        if loop == 0
            penal = init_p;
        else
            penal = min(penal + step_p, target_p);
        end
    end

    %% Optimization
    [c, vol, dc, dv] = analyzeI(xPhys) ;

    xval = x(:);
    xmin = max(0,xval-ml);
    xmax = min(1,xval+ml);

    f0val = c ; % To change the functional add a term to this and the gradient line.
    df0dx = dc(:);
    fval  = vol/(volfrac*nely*nelx*4) - 1;
    dfdx  = dv(:)'/(volfrac*nely*nelx*4);

    if printing 
        fprintf('Optm    It:%4i  Obj:%7.3f  Cons: %4.3f  KKT-norm.:%7.3f  RChange: %4.3f\n', ...
            loop,c,fval,kktnorm,rel_obj_change);
    end
    
    if 0 < loop
        rel_obj_change = abs((f0val - old_obj)/old_obj);

        if(abs(f0val-old_obj)<=TOLF*(1+f0val) && penal >= target_p)
            disp('TOLF requirement satisfy')
            break
        end
    end

    % The residual vector of the KKT conditions is calculated:
    if 0 < loop
        [~,kktnorm,~] = ...
            kktcheck(ncons,ndv,xmma,ymma,zmma,lam,xsi,eta,muu,zet,s, ...
            xmin,xmax,df0dx,fval,dfdx,a0,ai,ci,d);
    end

    %%%% The MMA subproblem is solved at the point PARAM_VALUE
    [xmma,ymma,zmma,lam,xsi,eta,muu,zet,s,low,upp] = ...
        mmasub(ncons,ndv,loop,xval,xmin,xmax,xold1, ...
        xold2, f0val,df0dx,fval,dfdx,low,upp,a0,ai,ci,d);

    xnew = reshape(xmma,4,nely,nelx);

    % Filtering and projection
    xDis  = filterDensities(ord2pos(xnew)) ;
    xPhys = pos2ord(xDis) ;
%     xPhys = xnew ;

    % Update MMA
    xold2 = xold1(:);
    xold1 = x(:);
    x     = xnew;

    old_obj = f0val;
    f_star = c;

    % PLOT DENSITIES
    if plotting
        figure(1)
        title('Optimization')
        subplot(1,1,1)
        colormap(gray);
        fv = contourf(posX,posY,1-xDis,levels,'EdgeColor','none') ; caxis([0 1]); axis equal; axis off; drawnow;
    end
    loop = loop + 1;
end
toc

if ~printing
    fprintf('Optm    It:%4i  Obj:%7.3f  Cons: %4.3f  KKT-norm.:%7.3f  RChange: %4.3f\n', ...
        loop,c,fval,kktnorm,rel_obj_change);
end

figure(1) 
title('Optimization') 
subplot(1,1,1)
xG = xGraph(xPhys,10) ;
colormap(gray); imagesc(1-xG); caxis([0 1]); axis equal; axis off; drawnow;
