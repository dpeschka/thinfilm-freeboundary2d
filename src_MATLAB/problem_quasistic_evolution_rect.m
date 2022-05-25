%
% Mixed FE approach for an incompressible Stokes problem
%
clear all
addpath('fem/')
addpath('dof/')
addpath('io/')

nfe      =     2;   % degree of fe-space velocity

params.vol0     =   1.25;
params.mu       =   1.0;
params.sigma    =   1.0;
params.g1       =   0.0;
params.g2       =   0.0;
params.delta    =   0.1;

% read mesh
[x,y,npoint,nelement,e2p,~,~] = readtriamesh('mesh/domain_rectangle');
x  = 2*x;
x  = x - 0.5;

% construct quadrature rules & FE spaces
qr1D = fem_quadrature_rule( 1 , 3);
qr2D = fem_quadrature_rule( 2 , 5);

FEb = fem_functions_1d(nfe, qr1D);
FE  = fem_functions_2d(nfe, qr2D);

[e2pb,~] = extract_e2p_boundary(e2p); 

ymin = min(y(:));
ymax = max(y(:));
yyy  = mean(y(e2pb(:,1:2)),2);
xxx  = mean(x(e2pb(:,1:2)),2);
selb = (abs(yyy-ymin)<1d-6)|(abs(yyy-ymax)<1d-6);
sell = (abs(x(1:npoint)+0.5)<1d-6);
selr = (abs(x(1:npoint)-0.5)<1d-6);
                                          

[dof  ,ndof  ,~  ] = generate_dofs(e2p , FE ,  true);
[dofb ,ndofb ,dm ] = generate_dofs(e2pb         , FEb,  true );

[dofb1,ndofb1,dm1] = generate_dofs(e2pb( selb,:), FEb,  false);
[dofb2,ndofb2,dm2] = generate_dofs(e2pb(~selb,:), FEb,  false);
[dofb3,ndofb3,dm3] = generate_dofs(e2pb(~selb,:), FEb,  true );

x = x( 1:ndof);
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

FEobj.selb     = selb;

FEobj.dofb1    = dofb1;
FEobj.ndofb1   = ndofb1;
FEobj.dm1      = dm1;

FEobj.dofb2    = dofb2;
FEobj.ndofb2   = ndofb2;
FEobj.dm2      = dm2;

FEobj.dofb3    = dofb3;
FEobj.ndofb3   = ndofb3;
FEobj.dm3      = dm3;

FEobj.idp      = idp;
FEobj.npoint   = max(max(e2p(:,1:3)));
FEobj.nelement = size(e2p,1);
FEobj.e2p      = e2p;
FEobj.e2pb     = e2pb;

%%
x0 = x;
y0 = y;
x  = x + 0.05*x.*cos(2/5*pi*y);

% iterate
t      = 0;
Nmax   = 1200;
tau    = 0.01;
[h] = geth_rect(x,y,params,FEobj);

trisurf(e2p(:,1:3),x,y,h)

[temp,iis]=sort(y(sell));

%%
for nt=1:Nmax
  fprintf('%i of %i\n',nt,Nmax)
  t = t + tau;
  x0 = x;
  y0 = y;
  
  % RICH2 discretization
  [x1,y1,h1] = evolve_rect(x,y,tau,params,FEobj);
  [x2,y2,h2] = nevolve_rect(x,y,tau/2,2,params,FEobj);
  x = 2*x2-x1;
  y = 2*y2-y1;
  h = 2*h2-h1;
  
  velx = (x-x0)/tau;
  vely = (y-y0)/tau;

  trisurf(FEobj.e2p(:,1:3),x,y,h)
  
  xlim([-1 1])
  ylim([0 5])
  zlim([0 2])
  
  xxb = x(sell);
  yyb = y(sell);
  xxb = xxb(iis);
  yyb = yyb(iis);
  
  hold on
  plot(xxb,yyb,'r-','LineWidth',2)
  hold off
  
  xdiff = abs(2*interp1(yyb,xxb,2.5));
  
  axis equal
  shading interp
  view(2)
  title(sprintf('%f t=%f',xdiff,t))
  
  drawnow
  tau = 0.005*xdiff/max(abs(velx)); 
end