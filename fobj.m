function Qp = fobj(PP)
global x0 z0 d0 Sint pframe DX MU
pp=pframe;
pp(:,6)=PP;
wp=0;
V2d=fwd(x0,z0,pp);
dc=V2d(:,1);
res=d0-dc;
vp=sum(abs(diff(PP)));
vd=norm(res);
Qp=vd+MU*vp;
return

%
%

