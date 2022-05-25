function [xn,yn,h] = evolve_rect(x,y,tau,params,FEobj)

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

selb   = FEobj.selb      ;
dofb1  = FEobj.dofb1     ;
ndofb1 = FEobj.ndofb1    ;
dm1    = FEobj.dm1       ;

dofb2  = FEobj.dofb2     ;
ndofb2 = FEobj.ndofb2    ;
dm2    = FEobj.dm2       ;

dofb3  = FEobj.dofb3     ;
ndofb3 = FEobj.ndofb3    ;
dm3    = FEobj.dm3       ;

%% FE part
% setup FE dofs/matrices
[edet ,dFinv ]=vec_transformation_2d(e2p          ,x,y,FE );
[edetb,dFinvb]=vec_transformation_2d(e2pb(~selb,:),x,y,FEb);

aa  = vec_localstiff  (edet,dFinv,FE);
mm  = vec_localmass   (edet      ,FE);
ccx = vec_localconv   (edet,dFinv,FE,ee,zz,dof);
ccy = vec_localconv   (edet,dFinv,FE,zz,ee,dof);

[ii   , jj  ] = distribute_dofs (dof ,FE );
[iib3 , jjb3] = distribute_dofs (dofb3,FEb);

% build standard matrices
A  = sparse(ii(:),jj(:),aa(:) ,ndof,ndof);
M  = sparse(ii(:),jj(:),mm(:) ,ndof,ndof);
Cx = sparse(ii(:),jj(:),ccx(:),ndof,ndof);
Cy = sparse(ii(:),jj(:),ccy(:),ndof,ndof);

[h] = geth_rect(x,y,params,FEobj);


%% ALE part
% build codim-1 problem with mean curvature
% and velocity from contact line model
hx = M\(Cx*h);
hy = M\(Cy*h);

aab     = vec_localstiff (edetb,dFinvb,FEb);
mmb     = vec_localmass  (edetb       ,FEb);

ix1 = find(dm3>0);
% ix2 = dm3(ix1);

xt = x(ix1);
yt = y(ix1);

[ux,uy]=flowrule(hx,hy,sigma);

uxt = ux(ix1);
uyt = uy(ix1);

Ab = sparse(iib3(:),jjb3(:),aab(:),ndofb3,ndofb3);
Mb = sparse(iib3(:),jjb3(:),mmb(:),ndofb3,ndofb3);

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
idp2 = unique(dofb2(:));nbd2=length(idp2);
idp1 = unique(dofb1(:));nbd1=length(idp1);
iib = [1:nbd2 1:nbd2 nbd2+(1:nbd1)];

jjb = [idp2;ndof+idp2;ndof+idp1];
aab = [nx(idp2);ny(idp2);ones(nbd1,1)];
nbd = nbd1+nbd2;
B   = sparse(iib(:),jjb(:),aab(:),nbd,2*ndof);

uxmod = ux;
uymod = uy;
uxmod(idp2) = uxn;
uymod(idp2) = uyn;
uymod(idp1) = 0;

% compute ALE extension
rhse = [zeros(2*ndof,1);B*[uxmod;uymod]];
Ae   = [A B';B sparse(nbd,nbd)];
ue = Ae\rhse;

ux_curv = ue(1:ndof);
uy_curv = ue(ndof+(1:ndof));

% update mesh
xn = x + tau*ux_curv;
yn = y + tau*uy_curv;