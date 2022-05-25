function [ux,uy]=flowrule(hx,hy,sigma)

h2 = hx.^2+hy.^2;

% vel = (1./sqrt(h2)).*(h2/2-sigma);
vel = (h2/2-sigma);

ux = -vel.*hx./sqrt(h2);
uy = -vel.*hy./sqrt(h2);