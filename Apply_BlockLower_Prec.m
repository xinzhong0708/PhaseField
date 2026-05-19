function z = Apply_BlockLower_Prec(r,Nphi,Lphi,Uphi,Lmu,Umu,A21,A11,A22)
%APPLY_BLOCKLOWER_PREC
%
% Apply approximate inverse of block lower preconditioner:
%
%   P = [Aphi   0
%        Amuphi Amu]
%
% so:
%   zphi = Aphi^{-1} rphi
%   zmu  = Amu^{-1} (rmu - Amuphi*zphi)

n = numel(r);

rphi = r(1:Nphi);
rmu  = r(Nphi+1:n);

% ------------------------------------------------------------
% Solve phi block
% ------------------------------------------------------------
if Nphi > 0

    if ~isempty(Lphi)
        zphi = Uphi \ (Lphi \ rphi);
    else
        d = diag(A11);
        zphi = rphi ./ max(abs(d),eps);
    end

else

    zphi = zeros(0,1);

end

% ------------------------------------------------------------
% Solve mu block with lower coupling correction
% ------------------------------------------------------------
rmu_eff = rmu - A21*zphi;

if ~isempty(Lmu)
    zmu = Umu \ (Lmu \ rmu_eff);
else
    d = diag(A22);
    zmu = rmu_eff ./ max(abs(d),eps);
end

z = [zphi; zmu];

end