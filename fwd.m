function  V=fwd(x0,z0,pp)
% Gravity field evaluation from 2D basin models
% constant density contrast assumed
%
% INPUT
%(x0,z0,xx,zz,Dens)
% pp=[1xc 2DUMMY 3Lx 4DUMMY 5zt  6zb 7DUMMY 8DUMMY 9DUMMY 10DUMMY 11dens]
%
% OUTPUT
% v = [gz gzx gzz gx]
%
np=length(pp(:,1));
nx=length(x0);
V=zeros(nx,4);
for k=1:np
    xx=[pp(k,1)-pp(k,3)/2 pp(k,1)+pp(k,3)/2];
    zz=[pp(k,5)  pp(k,6)];
    Dens=pp(k,11);
    v=mod2dg(x0,z0,xx,zz,Dens);
    V=V+v;
end
%
% |--.------------------------
% |
% |
% | xx(1)  xx(2)
% |    _____.......zz(1) 
% |    |   |
% |    |   |
% |    |___|.......zz(2)
% |
% |
%
%  Gravity anomaly from right-sided 2D prisms
%  input
%	xx - vertex x-postions (km)
%	zz - vertex z-positions (km)
%   Density contrast- g/cm3
%	obs- xx e zz sao vetores com dois elementos
%            (x0,z0) sao vetores com a n coordenadas das estacoes
%  output - gravity anomaly (mGal)
%           gravity derivatives (mGal/km)
%
function v=mod2dg(x0,z0,xx,zz,Dens)
hgg=0.0001;			%gravimeter height above the ground (10 cm);
ud=ones(size(x0));
x1=x0-xx(1)+eps;x1=x1.*ud*100000;
x2=x0-xx(2)+eps;x2=x2.*ud*100000;
z1=z0-zz(1)+eps;z1=z1.*ud*100000;
z2=z0-zz(2)+eps;z2=z2.*ud*100000;
x=x2;z=z2;
 W=log(x.*x+z.*z);
  d1=x.*W+2*z.*atan(x./z);
  D1=z.*W+2*x.*atan(z./x);
 dx1=W;
 dz1=atan(x./z);
x=x2;z=z1;
 W=log(x.*x+z.*z);
  d2=x.*W+2*z.*atan(x./z);
  D2=z.*W+2*x.*atan(z./x);
 dx2=W;
 dz2=atan(x./z);
x=x1;z=z2;
 W=log(x.*x+z.*z);
  d3=x.*W+2*z.*atan(x./z);
  D3=z.*W+2*x.*atan(z./x);
 dx3=W;
 dz3=atan(x./z);
x=x1;z=z1;
 W=log(x.*x+z.*z);
  d4=x.*W+2*z.*atan(x./z);
  D4=z.*W+2*x.*atan(z./x);
 dx4=W;
 dz4=atan(x./z);
 gama=0.667e-7;
 d =  -gama*Dens*(d1-d2-d3+d4)*1000; 		    % z-component (mGal)
 D =  -gama*Dens*(D1-D2-D3+D4)*1000;			% x-component (mGal)
 dx=  -gama*Dens*(dx1-dx2-dx3+dx4)*100000000;	% x-derivative of the gz component (mGal/km)    
 dz=-2*gama*Dens*(dz1-dz2-dz3+dz4)*100000000; 	% z-derivative of the gz component (mGal/km)
 v=[d dx dz D];
 return
