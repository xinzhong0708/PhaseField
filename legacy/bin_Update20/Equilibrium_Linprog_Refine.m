function OUT = Equilibrium_Linprog_Refine(T,P,phs_name,Nsys,options)
%EQUILIBRIUM_LINPROG_REFINE Cation-basis linprog equilibrium with refinement.
%
% Usage:
%
%   OUT = Equilibrium_Linprog_Refine(T,P,phs_name,Nsys)
%   OUT = Equilibrium_Linprog_Refine(T,P,phs_name,Nsys,options)
%
% Inputs:
%   T          temperature, K
%   P          pressure, Pa
%   phs_name   phase names, e.g. {'Orthopyroxene','Garnet','Kyanite','Quartz'}
%   Nsys       cation bulk composition:
%                [Fe Mg Ca Al] or [Fe Mg Ca Al Si]
%
% Output OUT:
%   OUT.phase_prop      phase proportions on cation basis
%   OUT.phase_comp      endmember proportions for each phase
%   OUT.alph            pseudocompound amount on cation basis
%   OUT.g               pseudocompound Gibbs energy on cation basis
%   OUT.Npc             pseudocompound composition on cation basis
%   OUT.pc_id           phase id of each pseudocompound
%   OUT.p_ref           refined endmember grid
%   OUT.mass_residual   cation mass-balance residual
%
% Notes:
%   The minimization basis is cation-normalized:
%       g   = g_raw / cation_sum
%       Npc = Npc_cation / cation_sum
%
%   O is not used as an independent constraint.
%   This avoids the formula-unit vs cation-basis inconsistency.

if nargin < 5 || isempty(options)
    options = struct();
end

% ------------------------------------------------------------
% Defaults
% ------------------------------------------------------------
solmod      = 'solution_models_PFM';
Cname       = {'Fe','Mg','Ca','Al','Si','O'};

dz          = 1/4;
nref        = 8;
active_tol  = 1e-7;
active_rel  = 1e-5;

z_window    = 0.085;
dz_fact     = 1.5;
ref_fact    = 1.25;

plot_ref    = 0;
verbose     = 1;

if isfield(options,'solmod'),     solmod     = options.solmod;     end
if isfield(options,'Cname'),      Cname      = options.Cname;      end
if isfield(options,'dz'),         dz         = options.dz;         end
if isfield(options,'nref'),       nref       = options.nref;       end
if isfield(options,'active_tol'), active_tol = options.active_tol; end
if isfield(options,'active_rel'), active_rel = options.active_rel; end
if isfield(options,'z_window'),   z_window   = options.z_window;   end
if isfield(options,'dz_fact'),    dz_fact    = options.dz_fact;    end
if isfield(options,'ref_fact'),   ref_fact   = options.ref_fact;   end
if isfield(options,'plot_ref'),   plot_ref   = options.plot_ref;   end
if isfield(options,'verbose'),    verbose    = options.verbose;    end

nph = numel(phs_name);

if isscalar(z_window)
    z_window = z_window*ones(1,nph);
end

if isscalar(dz_fact)
    dz_fact = dz_fact*ones(1,nph);
end

if numel(z_window) ~= nph
    error('options.z_window must be scalar or length(phs_name).')
end

if numel(dz_fact) ~= nph
    error('options.dz_fact must be scalar or length(phs_name).')
end

% ------------------------------------------------------------
% Bulk composition
% ------------------------------------------------------------
Nsys = Nsys(:);

if numel(Nsys) == 4
    Nsys = [Nsys; 1 - sum(Nsys)];
elseif numel(Nsys) ~= 5
    error('Nsys must be [Fe Mg Ca Al] or [Fe Mg Ca Al Si].')
end

if any(Nsys < -1e-12)
    error('Nsys contains negative component.')
end

Nsys(abs(Nsys) < 1e-14) = 0;
Nsys = Nsys./sum(Nsys);

% ------------------------------------------------------------
% Initialize thermodynamic models
% ------------------------------------------------------------
td = init_thermo(phs_name,Cname,solmod);

for iph = 1:nph
    if isfield(td(iph),'dz')
        td(iph).dz(:) = dz;
    end
end

p_ref = props_generate(td);

% ------------------------------------------------------------
% Refinement loop
% ------------------------------------------------------------
g_best     = [];
Npc_best   = [];
pc_id_best = [];
alph_best  = [];
gmin_best  = inf;

for iref = 1:nref

    % Evaluate pseudocompounds
    [g_raw,Npc_raw,pc_id] = Eval_Pseudocompounds_Local(T,P,phs_name,td,p_ref);

    % Convert to cation basis
    [g_cat,Npc_cat,cat_sum] = Normalize_Cation_Basis_Local(g_raw,Npc_raw);

    % Linprog
    LB = zeros(numel(g_cat),1);
    UB = inf(numel(g_cat),1);

    [alph,gmin,exitflag] = linprog(g_cat,[],[],Npc_cat,Nsys,LB,UB);

    if exitflag <= 0 || isempty(alph)
        error('linprog failed during refinement step %d.',iref)
    end

    alph = alph(:);
    alph(abs(alph) < 1e-14) = 0;

    % On cation basis sum(alph) should be 1, but normalize for roundoff.
    alph = alph./sum(alph);

    res = Npc_cat*alph - Nsys;

    if verbose == 1
        fprintf('ref %2d: gmin = %.12e, max mass residual = %.3e, active = %d\n', ...
            iref,gmin,max(abs(res)),nnz(alph > active_tol))
    end

    g_best     = g_cat;
    Npc_best   = Npc_cat;
    pc_id_best = pc_id(:);
    alph_best  = alph;
    gmin_best  = gmin;

    if iref == nref
        break
    end

    % Refine around active pseudocompounds
    p_old = p_ref;

    for iph = 1:nph

        idp = find(pc_id_best == iph);

        if isempty(idp)
            continue
        end

        a_phase = alph_best(idp);
        p_phase = p_ref{iph};

        if isempty(p_phase)
            continue
        end

        atol = max(active_tol,active_rel*max(a_phase));

        active_local = find(a_phase > atol);

        if isempty(active_local)
            [amax,imax] = max(a_phase);
            if amax > active_tol
                active_local = imax;
            end
        end

        if isempty(active_local)
            continue
        end

        zwin = z_window(iph)/(ref_fact^(iref-1));
        step = zwin/dz_fact(iph);

        z_new = [];

        for ia = active_local(:).'

            z0 = p_phase(ia,:);
            zloc = Local_Simplex_Refine_Local(z0,step,zwin);
            z_new = [z_new; zloc]; %#ok<AGROW>

        end

        p_ref{iph} = Unique_Simplex_Rows_Local([p_phase; z_new]);

    end

    if plot_ref == 1
        Plot_Refinement_Local(p_old,p_ref,phs_name,iref)
        drawnow
    end

end

% ------------------------------------------------------------
% Postprocess phase proportions and compositions
% ------------------------------------------------------------
phase_prop = zeros(1,nph);
phase_comp = cell(1,nph);

for iph = 1:nph

    idp = find(pc_id_best == iph);

    if isempty(idp)
        phase_prop(iph) = 0;
        phase_comp{iph} = [];
        continue
    end

    phase_prop(iph) = sum(alph_best(idp));

    if phase_prop(iph) > 0
        phase_comp{iph} = alph_best(idp).' * p_ref{iph} / phase_prop(iph);
    else
        phase_comp{iph} = zeros(1,size(p_ref{iph},2));
    end

end

% ------------------------------------------------------------
% Output
% ------------------------------------------------------------
OUT = struct();

OUT.T             = T;
OUT.P             = P;
OUT.phs_name      = phs_name;
OUT.Cname         = Cname(1:5);
OUT.Nsys          = Nsys(:).';

OUT.phase_prop    = phase_prop;
OUT.phase_comp    = phase_comp;

OUT.alph          = alph_best;
OUT.g             = g_best;
OUT.Npc           = Npc_best;
OUT.pc_id         = pc_id_best;
OUT.p_ref         = p_ref;
OUT.gmin          = gmin_best;
OUT.mass_residual = Npc_best*alph_best - Nsys;

OUT.cat_sum       = cat_sum;

% ------------------------------------------------------------
% Print result
% ------------------------------------------------------------
if verbose == 1

    fprintf('\nEquilibrium result at T = %.6g K, P = %.6g Pa\n',T,P)
    fprintf('max cation mass residual = %.6e\n',max(abs(OUT.mass_residual)))

    fprintf('\nPhase proportions:\n')
    for iph = 1:nph
        if phase_prop(iph) > 1e-8
            fprintf('  %-18s %.10f\n',phs_name{iph},phase_prop(iph))
        end
    end

    fprintf('\nPhase endmember compositions:\n')
    for iph = 1:nph
        if phase_prop(iph) > 1e-8
            fprintf('\n%s\n',phs_name{iph})
            disp(phase_comp{iph})
        end
    end

end

end


% =========================================================================
% Evaluate pseudocompounds
% =========================================================================
function [g,Npc,pc_id] = Eval_Pseudocompounds_Local(T,P,phs_name,td,p_ref)

[g0,v0]       = tl_g0(T,P,td);
[g,Npc,pc_id] = tl_gibbs_energy(T,P,phs_name,td,p_ref,g0,v0);

g     = g(:);
pc_id = pc_id(:);

end


% =========================================================================
% Convert raw formula-unit basis to cation-normalized basis
% =========================================================================
function [g_cat,Npc_cat,cat_sum] = Normalize_Cation_Basis_Local(g_raw,Npc_raw)

% Cation rows are Fe Mg Ca Al Si.
% O row is ignored as dependent.
Npc_cation = Npc_raw(1:end-1,:);

cat_sum = sum(Npc_cation,1);

if any(cat_sum <= 0)
    error('Normalize_Cation_Basis_Local: some pseudocompounds have zero cation sum.')
end

Npc_cat = Npc_cation ./ cat_sum;
g_cat   = g_raw(:) ./ cat_sum(:);

end


% =========================================================================
% Local simplex refinement around one active composition
% =========================================================================
function Z = Local_Simplex_Refine_Local(z0,step,zwin)

z0 = z0(:).';
z0(z0 < 1e-14) = 0;
z0 = z0./sum(z0);

n = numel(z0);

Z = z0;

if n <= 1
    return
end

if step <= 0 || zwin <= 0
    return
end

nmove = max(1,floor(zwin/step));

for imove = 1:nmove

    s = imove*step;

    for i = 1:n
        for j = 1:n

            if i == j
                continue
            end

            if z0(j) >= s

                z = z0;
                z(i) = z(i) + s;
                z(j) = z(j) - s;

                if all(z >= -1e-12) && abs(sum(z)-1) < 1e-10
                    z(z < 0) = 0;
                    z = z./sum(z);
                    Z = [Z; z]; %#ok<AGROW>
                end
            end
        end
    end
end

end


% =========================================================================
% Unique rows on simplex
% =========================================================================
function Z = Unique_Simplex_Rows_Local(Z)

if isempty(Z)
    return
end

Z(Z < 1e-14) = 0;

s = sum(Z,2);
good = s > 0;
Z = Z(good,:);
s = s(good);

Z = Z./s;

tol = 1e-10;
Zr = round(Z/tol)*tol;

[~,ia] = unique(Zr,'rows','stable');
Z = Z(ia,:);

end


% =========================================================================
% Optional plot
% =========================================================================
function Plot_Refinement_Local(p_old,p_new,phs_name,iref)

figure(100)
clf

nph = numel(phs_name);

for iph = 1:nph

    subplot(1,nph,iph)
    hold on

    A = p_old{iph};
    B = p_new{iph};

    if size(A,2) >= 2
        plot(A(:,1),A(:,2),'ko','MarkerSize',4)
        plot(B(:,1),B(:,2),'r.','MarkerSize',8)
        xlabel('z_1')
        ylabel('z_2')
    else
        plot(A(:,1),zeros(size(A,1),1),'ko','MarkerSize',4)
        plot(B(:,1),zeros(size(B,1),1),'r.','MarkerSize',8)
        xlabel('z_1')
    end

    title(sprintf('%s ref %d',phs_name{iph},iref))
    box on

end

end