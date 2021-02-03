function [M,S] = pp_kgm(np,pp)
%
% Mass and cross section evaluation from basin model in pp
% pp=[1xc 2NaN 3Lx 4NaN 5zt  6zb 7NaN 8NaN 9NaN 10NaN 11dens]
M=0; %[kg]
S=0; %[m2]
for k=1:np
    ap=pp(k,3)*[pp(k,6)-pp(k,5)];
    ap=ap*1e6;
    M=M+ap*pp(k,11)*1000;
    S=S+ap;
end


