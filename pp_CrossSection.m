function [x,z] = pp_CrossSection(np,pp,px)
x1=pp(1,1)-px/2;
x=[x1;x1];
z1=0.0;z2=pp(1,6);
z=[z1;z2];
for k=2:np
    x1=pp(k,1)-px/2;
    x=[x;x1;x1];
    z1=pp(k-1,6);z2=pp(k,6);
    z=[z;z1;z2];
end
x1=pp(np,1)+px/2;
x=[x;x1;x1];
z=[z;z2;0.0];
end

