clear
clc
global x0 z0 d0 Sint pframe DX MU
%
% 1) Data inversion parameters
nP=43;          % number of prisms
dns=-0.25;      %[gcm-3]
NoiseLevel=0.4; %[mGal]
zMax=6;         %[km] Maximum depth for the basin 
MU=0.1;        % regularization parameter
xbi=10;
xbf=50;         % left and wright limits of the basin
%
% 2) Data input
load dc.dat
x0=dc(:,1);
z0=dc(:,2);
nx=length(x0);
d0=dc(:,5);
d0=d0+[randn(nx,1)-0.5]*2*NoiseLevel;
dd=mean(diff(x0));
%
% 3) Mass excess and cross-section evaluation
Mint=2*1.1924*1e7*dd*sum(d0);
Sint=Mint/(dns*1000);
%
% 4) Basin model parametrization
DX=(xbf-xbi)/nP;
vv=[xbi:DX:xbf-DX]';
%          pp=[1xc 2NaN 3Lx 4NaN 5zt  6zb 7NaN 8NaN 9NaN 10NaN 11dens]
pframe=repmat([0.0 NaN   DX  NaN 0.0  0.0  NaN  NaN  NaN   NaN   dns],nP,1);
pframe(1:nP,1)=vv(1:nP)+DX/2;
%
% 5) Model constraints and inversion
wz=1e-6*Sint/(xbf-xbi);
P0=ones(nP,1)*wz;
lb=zeros(nP,1);
ub=ones(nP,1)*zMax;
Aeq=ones(1,nP)*DX; %[segment width in km]
beq=Sint*1e-6;     %[cross-section in km2]
tic
[X,FVAL]=fmincon(@fobj,P0,[],[],Aeq,beq,lb,ub,[]);
toc
%
% 6) Model response evaluation
ps=pframe;
ps(:,6)=X;
V2d=fwd(x0,z0,ps);
dc=V2d(:,1);
%
% 7) Figures
figure
subplot(211)
plot(x0,d0,'ok','MarkerFaceColor',0.8*[1 1 1])
ylabel('Gravity Data (mGal)')
hold on;plot(x0,dc,'-r',x0,d0-dc,'--k','LineWidth',2);hold off
ww=axis;;ww(3)=-40;ww(4)=5;axis(ww);
subplot(212)
[xv,zv]=pp_CrossSectionTRUE;
fill(xv,zv,0.8*[1 1 1],'EdgeColor',0.8*[1 1 1],'LineWidth',1)
[xm,zm]=pp_CrossSection(nP,ps,DX);
hold on;plot(xm,zm,'-r','LineWidth',2);hold off
uu=axis;uu(1:2)=ww(1:2);uu(4)=5;axis(uu);axis ij
xlabel('Distance (km)')
ylabel('Depth (km)')
%
[Msol,Ssol]=pp_kgm(nP,ps);
load pp.mod -ascii
load np.mod -ascii
[Mmod,Smod]=pp_kgm(np,pp);
text(45,2,['S_{data}=' num2str(Ssol*1e-6,'%-5.1f') 'km^2'])
text(45,3,['S_{model}=' num2str(Smod*1e-6,'%-5.1f') 'km^2'])
dS=sqrt(abs(Sint-Smod))/1000;
hold on;fill([0 dS dS 0],[0 0 dS dS],'k');hold off