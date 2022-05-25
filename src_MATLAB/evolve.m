function [xn,yn,h] = evolve(x,y,tau,params,FEobj)

ee = 1+0*x;
zz =   0*x;

sigma = params.sigma ;
g1    = params.g1    ;
g2    = params.g2    ;
delta = params.delta ;
vol0  = params.vol0  ;
mu    = params.mu    ;

FE   = FEobj.FE  ;
FEb  = FEobj.FEb ;
dof  = FEobj.dof ;
ndof = FEobj.ndof;

dofb     = FEobj.dofb    ;
ndofb    = FEobj.ndofb   ;
dm       = FEobj.dm      ;
idp      = FEobj.idp     ;
npoint   = FEobj.npoint  ;
nelement = FEobj.nelement;
e2p      = FEobj.e2p     ;
e2pb     = FEobj.e2pb    ;

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

% build rhs, constraints, and solve
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

% build normal constrain
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
xn = x + tau*ux_curv;
yn = y + tau*uy_curv;