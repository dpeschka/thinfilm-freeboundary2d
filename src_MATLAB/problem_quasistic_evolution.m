%
% Mixed FE approach for strong contact line dissipation limit
% of thin-film free boundary problem
%
clear all

addpath('fem/')
addpath('dof/')
addpath('io/')

% different parameter settings
% w=1.0 a g1= 0 g2=3.25    delta=0.01/2
% w=0.8 b g1= 0 g2=3.25/w
% w=0.6 c g1= 0 g2=3.25/w
% w=1.0 d g1= 5 ...
% w=1.0 e g1=10
% w=1.0 f g1=20
% w=1.0 g g1=40
% w=1.0 h g1=80
% w=1.4 i g1=80
% w=1.0 j g1= 0 delta=0.01*2
%       k
w        =  1.0;
mu       =  1.0;
nfe      =    2;   % degree of fe-space velocity
vol0     =  1.0;
tau      = 0.01/(2*w);

sigma    = 1.00;
g1       = 0.00;
g2       = 3.25*w;
delta    = 0.01*2;

sigma    = 1.00;
g1       = 0.00;
g2       = 3.25*w;
delta    = 0.01*2;

params.sigma = sigma;
params.g1    = g1;
params.g2    = g2;
params.delta = delta;
params.vol0  = vol0;
params.mu    = mu;

%% init

% read mesh
[x,y,npoint,nelement,e2p,~,~] = readtriamesh('mesh/disc');
[e2pb,~] = extract_e2p_boundary(e2p);
x=x/3;
y=y/3;

%%
% construct quadrature rules & FE spaces
qr1D = fem_quadrature_rule( 1 , 3);
qr2D = fem_quadrature_rule( 2 , 5);

FEb = fem_functions_1d(nfe, qr1D);
FE  = fem_functions_2d(nfe, qr2D);

[dof ,ndof ,~ ] = generate_dofs(e2p , FE ,  true);
[dofb,ndofb,dm] = generate_dofs(e2pb, FEb,  true);

x = x(1:ndof);
y = y(1:ndof);

ee = 1+0*x;
zz =   0*x;
if nfe==1
    idp = e2pb(:,1:2);
else
    idp = e2pb;
end
idp = unique(idp(:));

%% iterate
t = 0;

for i=1:1400
    i
    t = t + tau;
    
    %% FE part
    % setup FE dofs/matrices
    [edet ,dFinv ]=vec_transformation_2d(e2p ,x,y,FE );
    [edetb,dFinvb]=vec_transformation_2d(e2pb,x,y,FEb);
    
    aa  = vec_localstiff  (edet,dFinv,FE);
    mm  = vec_localmass   (edet      ,FE);
    ccx = vec_localconv   (edet,dFinv,FE,ee,zz,dof);
    ccy = vec_localconv   (edet,dFinv,FE,zz,ee,dof);
    
    [ii ,jj ] = distribute_dofs (dof ,FE );
    [iib,jjb] = distribute_dofs (dofb,FEb);
    
    % build standard matrices
    A  = sparse(ii(:),jj(:),aa(:) ,ndof,ndof);
    M  = sparse(ii(:),jj(:),mm(:) ,ndof,ndof);
    Cx = sparse(ii(:),jj(:),ccx(:),ndof,ndof);
    Cy = sparse(ii(:),jj(:),ccy(:),ndof,ndof);
    
    % solve stationary droplet shape on a given domain with constraint
    % that droplet volume is given by vol0
    v   = M*ones(ndof,1);
    L   = [A+g1*M v;v' 0];
    rhs = [g2*M*x;vol0];
    nbd = length(idp);
    ii  = 1:nbd;jj = idp;
    B   = sparse(ii,jj,ones(nbd,1),nbd,ndof+1);
    u   = [L B';B sparse(nbd,nbd)]\[rhs;zeros(nbd,1)];
    h   = u(1:ndof);
    
    %% ALE part
    % build codim-1 problem with mean curvature
    % and velocity from contact line model
    hx = M\(Cx*h);
    hy = M\(Cy*h);
    
    aab     = vec_localstiff (edetb,dFinvb,FEb);
    mmb     = vec_localmass  (edetb       ,FEb);
    
    ix1 = find(dm>0);
    ix2 = dm(ix1);
    
    xt = x(ix1);
    yt = y(ix1);
    
    [ux,uy]=flowrule(hx,hy,sigma);
    
    uxt = ux(ix1);
    uyt = uy(ix1);
    
    Ab = sparse(iib(:),jjb(:),aab(:),ndofb,ndofb);
    Mb = sparse(iib(:),jjb(:),mmb(:),ndofb,ndofb);
    
    % mean curvature velocity
    uxn = (Mb+tau*delta*Ab)\(Mb*uxt - delta*Ab*xt);
    uyn = (Mb+tau*delta*Ab)\(Mb*uyt - delta*Ab*yt);
    
    % solve ALE problem, which optimizes the extension of uxn,uyn
    % but uses the normal component of uxn,uyn
    
    % compute normal
    nx = -hx./sqrt(hx.^2+hy.^2);
    ny = -hy./sqrt(hx.^2+hy.^2);
    
    % build extension problem
    Z = sparse(ndof,ndof);
    A = [A Z;Z A];
    
    % build normal constraint
    iib = [1:nbd 1:nbd];
    jjb = [idp;ndof+idp];
    aab = [nx(idp);ny(idp)];
    B = sparse(iib(:),jjb(:),aab(:),nbd,2*ndof);
    
    uxmod = ux;
    uymod = uy;
    uxmod(ix1) = uxn;
    uymod(ix1) = uyn;
    
    % compute ALE extension
    ue = [A B';B sparse(nbd,nbd)]\[zeros(2*ndof,1);B*[uxmod;uymod]];
    
    ux_curv = ue(1:ndof);
    uy_curv = ue(ndof+(1:ndof));
    
    % update mesh
    x = x + tau*ux_curv;
    y = y + tau*uy_curv;
    
    % output to paraview
    uuu{1}=h;
    uuu{2}=ux_curv;
    uuu{3}=uy_curv;
    vnames{1}='h';
    vnames{2}='ux';
    vnames{3}='uy';
    x0 = mean(x);
    %outvtk_triP2_scalar(sprintf('sliding_droplet%4.4i.vtk',i),x,y,e2p(:,1:6),uuu,vnames);
    
    % output to screen
    x0 = 0;
    trisurf(e2p(:,1:3),x(1:npoint)-x0,y(1:npoint),0*x(1:npoint),h(1:npoint))
    shading interp
    hold on
    quiver(x-x0,y,ux,uy,1)
    quiver(x-x0,y,ux_curv,uy_curv,1)
    hold off
    view(2)
    axis equal tight
    grid on
    drawnow    
end