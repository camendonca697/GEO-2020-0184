% synthetic data evaluation from basin models in figures 2, 3 and 4
%
% type and enter:
% prog1_SyntheticData
% prog2_DataInversion
%
clear
%
% 1) Data acquisition parameters
dd=1.1; % station interval along the profile (km) 
xi=0;
xf=60;  % profile starting at xi and ending at xf
x0=[xi:dd:xf]';
nx=length(x0);
zo=0.0000; % gravimeter height above the ground (km)
z0=zo*ones(nx,1);
%
% 2) Basin model
dns=-0.25; % density gcm-3
xbi=10;
xbf=50; % x-position with left and right endings of the basin
nMOD=input(' Your choice 2=Figure 2;  3=Figure 3; 4=Figure 4; 5= DIY ');
if nMOD==2
    zb=[0.1 0.2 0.4 0.4 0.6 0.8 1.0 1.2 2.2 3.2 4.2 4.5 4.5 3.5 0.5 0.45 0.40 0.25 0.25 0.20 0.10]';
elseif nMOD==3
    zb=[ 1.2 2.2 3.2 4.2 4.5 4.5 3.5 0.5 0.45 0.1 0.2 0.4 0.4 0.6 0.8 1.0 0.40 0.25 0.25 0.20 0.10]';
elseif nMOD==4
    zb=[ 0.20 0.10 1.2 2.2 3.2 4.2  0.45 0.1 0.2 0.4 0.4 0.6 4.5 4.5 3.5 0.5 0.8 1.0 0.40 0.25 0.25]';
else
    disp('Try your own model...');return
    zb=[0.2 0.5 1.5 2.5 3.5 4.0 3.5 1.0 0.5 0.5]
end
% zb= depth of the bottom of each piece-wise segment representing the interface 
np=length(zb);      % number of segments representing the interface
px=(xbf-xbi)/np;    % width of the segment 
%pp=[1xc        2NaN 3Lx 4NaN 5zt  6zb   7NaN 8NaN 9NaN 10NaN 11dens]
pp= [xbi+px/2    NaN  px  NaN 0.0  zb(1)  NaN  NaN  NaN   NaN   dns];
for k=2:np
    wx=xbi+(k-1)*px+px/2;
    wp=[wx   NaN  px NaN 0.0  zb(k) NaN  NaN  NaN  NaN dns];
    pp=[pp;wp];
end
%
% 3) Forward model evaluation
V2d=fwd(x0,z0,pp);
dc=V2d(:,1);
gInc=atand(V2d(:,1)./V2d(:,4));
% Gaussian noise 0.4 mGal
rnd=[randn(nx,1)-0.5]*2*0.4;
d0=dc+rnd;
%
% 4) Save synthetic data and model
vv=[x0 z0 d0 rnd dc gInc];
save dc.dat vv -ascii
save pp.mod pp -ascii
save np.mod np -ascii
save px.mod px -ascii
%
% 5) Mass evaluation from data and model
Mint=2*1.1924*1e7*dd*sum(d0); % Mass evaluation from the data set
Sint=Mint/(dns*1000);         % Cross-section from Mint (assuming known density dns)
[Mmod, Smod]=pp_kgm(np,pp);   % Mass (Mmod) and cross-section area (Smod) evaluation directly from the model
%
% 6) Gravity profile and basin picture
figure
subplot(211)
plot(x0,d0,'ok',x0,dc,'-k')
ylabel('Gravity Data (mGal)')
subplot(212)
[vx,vz]=pp_CrossSection(np,pp,px);
plot(vx,vz,'-r')
hold on;
fill(vx,vz,'y');
hold off
axis ij;
axis([xi xf 0 5])
text(1,2,['S_{data }=' num2str(Sint*1e-6,'%-5.1f') 'km^2'])
text(1,3,['S_{model}=' num2str(Smod*1e-6,'%-5.1f') 'km^2'])
text(39,2,['M_{data }=' num2str(Mint*1e-9,'%-5.2f') 'x10^{9} kg/m'])
text(39,3,['M_{model}=' num2str(Mmod*1e-9,'%-5.2f') 'x10^{9} kg/m'])
dS=sqrt(abs(Sint-Smod))/1000;
hold on;
fill([0 dS dS 0],[0 0 dS dS],'k');
hold off
xlabel('Distance (km)')
ylabel('Depth (km)')
