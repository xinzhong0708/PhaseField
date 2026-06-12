function PARAM = Compute_M(STATE,PARAM,MODEL,PHYS)
%COMPUTE_M Build heterogeneous elemental mobility from phase mobilities.
%
% PHYS.M_phs can be:
%
%   1) vector, length = Nphase:
%        same phase mobility for all elements
%
%   2) matrix, Nphase x Ne:
%        rows    = phases
%        columns = elements
%
%   3) matrix, Ne x Nphase:
%        rows    = elements
%        columns = phases
%
% Phase order follows MODEL.phs_name / phase index.
% Element order follows STATE.E / PARAM.Ne.
%
% At interfaces:
%   M_e(x) = sum_phase p_phase(x) * M_phase,e

[ny,nx,Ngrain] = size(STATE.p);

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
if isvector(Mraw)

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

        error('Compute_M: PHYS.M_phs must be Nphase x Ne or Ne x Nphase.')

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
% Build one mobility field for each element
% ------------------------------------------------------------
PARAM.M = cell(1,Ne);

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

    PARAM.M{ie} = Mtmp;

end

end