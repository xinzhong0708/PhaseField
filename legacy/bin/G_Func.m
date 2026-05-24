function [g_norm,g_tot] = G_Func(g0,td,c,T,P,add,nN,g0_N,E_sc,v)
%Load parameter
P             =  P/1e8;
mtpl          =  td.mtpl;
zt            =  td.zt;
alp           =  td.alp;
w             =  td.w;
%Real c and penalty c
c_real        =  c(1:length(g0));
c_pena        =  c(length(g0)+1:end);
%Calculate free energy
%Mechanical mixing
g_mech        =  c_real*g0';
%If melt
z             =  c_real*zt;
if strcmp(td.phase_name,'Melt(H18)')
    [mtpl,z,zt]   =  temkin_H18('Melt(H18)',td,c_real,z);
end
%Ideal mixing energy
g_id          =  T*8.3144*sum(mtpl.*(z.*softlog(z+double(z==0)) - c_real*(zt.*softlog(zt+double(zt==0)))),2); % Configurational Entropy
%Nonideal mixing energy
alp           = [1 T P]*alp;
W             =  w(:,:,1)   +   w(:,:,2)*T +   w(:,:,3)*P;
g_nid         =  0;
for i = 1:size(W,1)
    for j = 1:size(W,2)
        g_nid =  g_nid + c_real(:,i).*c_real(:,j).*alp(i).*alp(j)./(c_real*alp')./(alp(i)+alp(j)).*W(i,j);
    end
end

%Add mechanical+ideal+nonideal parts
g_tot         =  g_mech + g_id + g_nid;

%Add dependent endmember
if add
    deltaG    =  8.314*T*5;
    % g_tot     =  g_tot + c_pena*(g0_N'+deltaG) + 8.314*T*c_pena*softlog(c_pena');
    g_tot     =  g_tot + sum( g0_N.*c_pena + 2e7.*c_pena.^2 );
    n         = [td.n_em(:,1:end-1);nN];
else
    n         =  td.n_em(:,1:end-1);
end

%Normalize G to 1 cation
g_tot         =  g_tot/E_sc/v;
cation        =  c*sum(n,2);
g_tot         =  (g_tot);
g_norm        =  (g_tot/cation);
end

function y = softlog(x)
eps=1e-4;
y = log(x+sqrt(x.^2+eps^2))-log(2);
end


