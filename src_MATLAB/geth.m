function [h] = geth(x,y,params,FEobj)

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

