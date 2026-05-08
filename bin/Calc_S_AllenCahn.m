function STATE = Calc_S_AllenCahn(STATE,PARAM,MODEL)
%CALC_S_ALLENCAHN Allen-Cahn source term using structured variables.
phi = STATE.phi;
omg = STATE.omg;
LL  = PARAM.L;
Np  = size(phi,3);
s   = repmat({zeros(size(phi,1),size(phi,2))},1,Np);
for alp = 1:Np
    s{alp} = -LL .* MODEL.dgdphi(phi(:,:,alp));
    for beta = 1:Np
        s{alp} = s{alp} - LL .* MODEL.dpdphi(alp,beta,phi) .* omg(:,:,beta);
    end
end
STATE.S_AC = s;
end