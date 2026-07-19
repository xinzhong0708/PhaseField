function PARAM = Compute_M(STATE,PARAM,MODEL,PHYS)
%COMPUTE_M Build heterogeneous diagonal elemental mobility from phase mobilities.
%
%
% Phase order follows MODEL.phs_name / phase index.
% Element order follows STATE.E / PARAM.Ne.
%
% At interfaces:
%   M_e(x) = sum_phase p_phase(x) * M_phase,e
%
% Output:
%   PARAM.M{ie,ie} = diagonal elemental mobility
%   PARAM.M{ie,je} = zero for ie ~= je

[ny,nx,~] = size(STATE.p);

if isfield(PARAM,'Ne')
    Ne = PARAM.Ne;
else
    Ne = numel(STATE.E);
end

if isfield(MODEL,'phs_name')
    Nphase = numel(MODEL.phs_name);
else
    Nphase = max(MODEL.phase_index);
end

if ~isfield(PHYS,'M_phs')
    error('Compute_M: PHYS.M_phs is missing.')
end

Mraw = PHYS.M_phs;

% ------------------------------------------------------------
% Interpret PHYS.M_phs
% ------------------------------------------------------------
if isscalar(Mraw)

    M_phase_elem = Mraw*ones(Nphase,Ne);

elseif isvector(Mraw)

    Mraw = Mraw(:);

    if numel(Mraw) ~= Nphase
        error('Compute_M: vector PHYS.M_phs must have length Nphase.')
    end

    M_phase_elem = repmat(Mraw,1,Ne);     % Nphase x Ne

else

    [n1,n2] = size(Mraw);

    if n1 == Nphase && n2 == Ne

        % rows = phase, columns = element
        M_phase_elem = Mraw;

    elseif n1 == Ne && n2 == Nphase

        % rows = element, columns = phase
        M_phase_elem = Mraw.';

    else

        error('Compute_M: PHYS.M_phs must be scalar, Nphase x Ne, or Ne x Nphase.')

    end
end

% ------------------------------------------------------------
% Collapse grain p to phase p
% ------------------------------------------------------------
p_phase = zeros(ny,nx,Nphase);

for iph = 1:Nphase

    grains = find(MODEL.phase_index == iph);

    if ~isempty(grains)
        p_phase(:,:,iph) = sum(STATE.p(:,:,grains),3);
    end
end

psum = sum(p_phase,3);

% ------------------------------------------------------------
% Allocate full elemental mobility matrix
% ------------------------------------------------------------
PARAM.M = cell(Ne,Ne);

for ie = 1:Ne
    for je = 1:Ne
        PARAM.M{ie,je} = zeros(ny,nx);
    end
end

PARAM.M_diag = cell(1,Ne);

% ------------------------------------------------------------
% Build diagonal mobility field for each element
% ------------------------------------------------------------
for ie = 1:Ne

    Mgrid = zeros(ny,nx);

    for iph = 1:Nphase
        Mgrid = Mgrid + M_phase_elem(iph,ie).*p_phase(:,:,iph);
    end

    % Safe normalization, usually psum is already 1
    mask = psum > eps;

    Mtmp = Mgrid;
    Mtmp(mask)  = Mgrid(mask)./psum(mask);
    Mtmp(~mask) = mean(M_phase_elem(:,ie));

    PARAM.M{ie,ie}  = Mtmp;
    PARAM.M_diag{ie} = Mtmp;

end

% Store phase-element mobility table for diagnostics
PARAM.M_phase_elem = M_phase_elem;

end