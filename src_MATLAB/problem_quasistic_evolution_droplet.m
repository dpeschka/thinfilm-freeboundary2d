%
% Mixed FE approach for strong contact line dissipation limit
% of thin-film free boundary problem
%
clear all

addpath('fem/')
addpath('dof/')
addpath('io/')

nfe      =     2;   % degree of fe-space velocity

params.vol0     =   1.0;
params.mu       =   1.0;
params.sigma    =  1.00;
params.g1       =  0.0;
params.g2       =  3.25;
params.delta    =  0.02;

% read mesh
%[x,y,npoint,nelement,e2p,~,~] = readtriamesh('mesh/disc_finer');
[x,y,npoint,nelement,e2p,~,~] = readtriamesh('mesh/disc');
[e2pb,~] = extract_e2p_boundary(e2p);
x=x/3;
y=y/3;

% construct quadrature rules & FE spaces
qr1D = fem_quadrature_rule( 1 , 3);
qr2D = fem_quadrature_rule( 2 , 5);

FEb = fem_functions_1d(nfe, qr1D);
FE  = fem_functions_2d(nfe, qr2D);

[dof ,ndof ,~ ] = generate_dofs(e2p , FE ,  true);
[dofb,ndofb,dm] = generate_dofs(e2pb, FEb,  true);

x = x(1:ndof);
y = y(1:ndof);

if nfe==1
    idp = e2pb(:,1:2);
else
    idp = e2pb;
end
idp = unique(idp(:));

FEobj.FE       = FE;
FEobj.FEb      = FEb;
FEobj.dof      = dof;
FEobj.ndof     = ndof;
FEobj.dofb     = dofb;
FEobj.ndofb    = ndofb;
FEobj.dm       = dm;
FEobj.idp      = idp;
FEobj.npoint   = npoint;
FEobj.nelement = nelement;
FEobj.e2p      = e2p;
FEobj.e2pb     = e2pb;

x0 = x;
y0 = y;

%% iterate
t      = 0;
Nmax   = 500;
tau    = 0.01;

SCHEME = 'RICH2';
%%
for nt=1:Nmax
    fprintf('%i, %i of %i\n',0,nt,Nmax)
    t = t + tau;
    if SCHEME == 'SEMI1'
        [x,y,h] = evolve(x,y,tau,params,FEobj);
    end
    if SCHEME == 'RICH2'
        [x1,y1,h1] = evolve(x,y,tau,params,FEobj);
        [x2,y2,h2] = nevolve(x,y,tau/2,2,params,FEobj);
        x = 2*x2-x1;
        y = 2*y2-y1;
        h = 2*h2-h1;
    end
    if SCHEME == 'RICH3'
        [x1,y1,h1] = evolve(x,y,tau,params,FEobj);
        [x2,y2,h2] = nevolve(x,y,tau/2,2,params,FEobj);
        [x4,y4,h4] = nevolve(x,y,tau/4,4,params,FEobj);
        a1 = +8/3;
        a2 = -6/3;
        a3 = +1/3;
        x = a1*x4+a2*x2+a3*x1;
        y = a1*y4+a2*y2+a3*y1;
        h = a1*h4+a2*h2+a3*h1;
    end
    
    % x = xn;
    % y = yn;
    hh = geth(x,y,params,FEobj);
    trisurf(FEobj.e2p(:,1:3),x,y,hh)
    hold on
    plot(x,y,'r.')
    hold off
    %drawnow
    view(2)
    axis equal
    drawnow
    
end