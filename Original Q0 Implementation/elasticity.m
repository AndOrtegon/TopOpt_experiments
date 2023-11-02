clear;
clc;

addpath(genpath('util'))
%% Declaration of global variable
global alldofs nelx nely volfrac penal He rmin ft E0 Emin KE F U H Hs edofMat freedofs iK jK ml penal_type Hproj Hbeta Heta problem
global method plotting printing TOLF
method = 'direct' ; 

%% Optimization input
% Problem 'Elasticity' or 'Heat'
problem = 'Elasticity' ;

nelx    = 120;          % Number of elements in x
nely    = 40;           % Number of elements in y
volfrac = 0.4;          % Maximum fraction of material 
penal   = 3;            % Value of the penalization
rmin    = 2;            % Radius of the filter
ft      = 2;            % Type or filter
step    = 0.02 ;        % Continuation step, for no continuation set to (penal-1) in SIMP and penal in RAMP
penal_type = 'SIMP';    % Type of penalization

TOLF = 1e-7 ;           % Optimization tolerance
plotting= true ;        % Plot design during optimization
printing = true ;       % Print values in the console during optimization


% Move limit
ml      = 0.2;   %1.0 means essentially no move limits

% Heaviside projection
Hproj   = false;
Hbeta   = 8;
Heta    = 0.5;

%% MATERIAL PROPERTIES
E0   = 1;
Emin = 1e-9;
nu   = 0.3;

% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];

KE      = 1/(1-nu^2)*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11])/24;
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK      = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK      = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

%%% Hay tres problemas posibles para el caso del compliance,
%%% comentar/descomentar para ver

% Half-MMB Beam
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);

% 2d cantilever
%     F = sparse(2*(nely+1)*(nelx+1),1,-1,2*(nely+1)*(nelx+1),1);
%     U = zeros(2*(nely+1)*(nelx+1),1);
%     fixeddofs = 1:1:2*(nely+1) ;

% Define loads and supports: Mitchell Bridge
%     F = sparse(2*(nely+1)*(nelx/2+1),1,-1,2*(nely+1)*(nelx+1),1);
%     U = zeros(2*(nely+1)*(nelx+1),1);
%     fixeddofs = union([2*(nely+1),2*(nely+1)-1],[2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1)-1]);


alldofs  = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);

%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k  = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
He = spdiags(1./Hs,0,nely*nelx,nely*nelx)*H ;
%% Initialize Design
% Random initial conditions
%rng(1);
% x = rand(nely,nelx);

% Uniform initial conditions
x = volfrac*ones(nely,nelx);

%% Optimization phase
% First optimization
[x_star,f_star] = optimization_phase(x,step);
