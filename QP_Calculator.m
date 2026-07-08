function [c_new,mu_e,chi,DIAG] = QP_Calculator(pars,p,c_old,E,eta,opts)
%QP_CALCULATOR Local pseudocompound QP equilibrium prototype.
%
% This is a prototype replacement for LE_Calculator in interface / spinodal
% risk nodes.
%
% Main idea:
%
%   Each active phase alpha is represented by two local pseudocompounds:
%
%       c_alpha_minus = c_old - dc*v
%       c_alpha_plus  = c_old + dc*v
%
%   Unknowns are their amounts x >= 0:
%
%       sum_q x_alpha_q = p_alpha
%
%   The local objective is:
%
%       min  sum_q x_q*g_q
%          + eta/2 * |A*x - E|^2
%          + prox/2 * sum_q x_q*|c_q - c_old|^2
%
%   where A columns are elemental compositions e_q.
%
% Notes:
%   - This avoids using the raw spinodal Hessian as Newton tangent.
%   - It can be used even for a pure cpx node next to quartz.
%   - It returns an effective single composition per phase by averaging
%     the local pseudocompounds.
%
% Expected input shape:
%   p      : Nnode x Nphase or Nphase x Nnode
%   E      : Ne x Nnode or Nnode x Ne
%   c_old  : cell{iph}, each nc x Nnode or Nnode x nc
%
% Output:
%   c_new  : same orientation as c_old
%   mu_e   : Ne x Nnode
%   chi    : cell{Ne,Ne}, each 1 x Nnode
%   DIAG   : diagnostics

if nargin < 6 || isempty(opts)
    opts = struct();
end

if ~isstruct(opts)
    opts = struct();
end

% ------------------------------------------------------------
% Defaults
% ------------------------------------------------------------
p_cut          = 1e-12;
split_dc       = 0.01;
prox_weight    = 1e2;
qp_reg         = 1e-14;
chi_floor      = 1e-8;
relax_c        = 1.0;
dc_eff_cap     = inf;
use_quadprog   = 1;
slack_weight   = 1e8;
use_neg_hessian_dir = 1;

if isfield(opts,'p_cut'),        p_cut        = opts.p_cut;        end
if isfield(opts,'split_dc'),     split_dc     = opts.split_dc;     end
if isfield(opts,'prox_weight'),  prox_weight  = opts.prox_weight;  end
if isfield(opts,'qp_reg'),       qp_reg       = opts.qp_reg;       end
if isfield(opts,'chi_floor'),    chi_floor    = opts.chi_floor;    end
if isfield(opts,'relax_c'),      relax_c      = opts.relax_c;      end
if isfield(opts,'dc_eff_cap'),   dc_eff_cap   = opts.dc_eff_cap;   end
if isfield(opts,'use_quadprog'), use_quadprog = opts.use_quadprog; end
if isfield(opts,'slack_weight'), slack_weight = opts.slack_weight; end
if isfield(opts,'use_neg_hessian_dir')
    use_neg_hessian_dir = opts.use_neg_hessian_dir;
end

Nphase = numel(pars);

% ------------------------------------------------------------
% Orient p as Nnode x Nphase
% ------------------------------------------------------------
if size(p,2) == Nphase
    pmat = p;
elseif size(p,1) == Nphase
    pmat = p.';
else
    error('QP_Calculator: p must be Nnode x Nphase or Nphase x Nnode.')
end

Nnode = size(pmat,1);

% ------------------------------------------------------------
% Orient E as Ne x Nnode
% ------------------------------------------------------------
if size(E,2) == Nnode
    Emat = E;
elseif size(E,1) == Nnode
    Emat = E.';
else
    error('QP_Calculator: E must be Ne x Nnode or Nnode x Ne.')
end

Ne = size(Emat,1);

% ------------------------------------------------------------
% Eta vector
% ------------------------------------------------------------
if isscalar(eta)
    eta_vec = eta*ones(1,Nnode);
else
    eta_vec = eta(:).';
    if numel(eta_vec) ~= Nnode
        eta_vec = eta_vec(1)*ones(1,Nnode);
    end
end

eta_vec = max(eta_vec,0);

% ------------------------------------------------------------
% Orient c_old as nc x Nnode for each phase
% ------------------------------------------------------------
c_cols     = cell(1,Nphase);
c_was_t    = false(1,Nphase);
c_new_cols = cell(1,Nphase);

for iph = 1:Nphase
    [c_cols{iph},c_was_t(iph)] = C_To_Cols_Local(c_old{iph},Nnode);
    c_new_cols{iph} = c_cols{iph};
end

% ------------------------------------------------------------
% Allocate output
% ------------------------------------------------------------
mu_e = zeros(Ne,Nnode);

chi = cell(Ne,Ne);
for ie = 1:Ne
    for je = 1:Ne
        chi{ie,je} = zeros(1,Nnode);
    end
end

DIAG = struct();
DIAG.exitflag   = zeros(1,Nnode);
DIAG.fail       = false(1,Nnode);
DIAG.nvar       = zeros(1,Nnode);
DIAG.mass_err   = zeros(Ne,Nnode);
DIAG.solver     = strings(1,Nnode);

% ------------------------------------------------------------
% Main node loop
% ------------------------------------------------------------
for in = 1:Nnode

    pn = pmat(in,:).';
    active = find(pn > p_cut);

    if isempty(active)

        DIAG.fail(in) = true;
        continue

    end

    En      = Emat(:,in);
    eta_n   = eta_vec(in);

    cand_phase = [];
    cand_c     = {};
    cand_g     = [];
    cand_e     = [];
    cand_d2    = [];

    % --------------------------------------------------------
    % Build two local pseudocompounds for each active phase
    % --------------------------------------------------------
    for ia = 1:numel(active)

        iph = active(ia);

        c0 = c_cols{iph}(:,in);

        v = Split_Direction_Local(pars{iph},c0,opts,use_neg_hessian_dir);
        v = v(:);

        if numel(v) ~= numel(c0)
            v = zeros(size(c0));
            v(1) = 1;
        end

        if max(abs(v)) > 0
            v = v./max(abs(v));
        else
            v(:) = 0;
            v(1) = 1;
        end

        cA = c0 - split_dc*v;
        cB = c0 + split_dc*v;

        [okA,gA,eA] = Try_Eval_Local(pars{iph},cA,Ne);
        [okB,gB,eB] = Try_Eval_Local(pars{iph},cB,Ne);

        % Backtrack if a candidate is invalid
        dc_try = split_dc;

        while (~okA || ~okB) && dc_try > split_dc*1e-4

            dc_try = 0.5*dc_try;

            cA = c0 - dc_try*v;
            cB = c0 + dc_try*v;

            [okA,gA,eA] = Try_Eval_Local(pars{iph},cA,Ne);
            [okB,gB,eB] = Try_Eval_Local(pars{iph},cB,Ne);

        end

        % If still invalid, fall back to old composition twice
        if ~okA || ~okB
            cA = c0;
            cB = c0;

            [ok0,g0,e0] = Try_Eval_Local(pars{iph},c0,Ne);

            if ~ok0
                DIAG.fail(in) = true;
                continue
            end

            gA = g0; gB = g0;
            eA = e0; eB = e0;
        end

        cand_phase = [cand_phase; iph; iph]; %#ok<AGROW>

        cand_c{end+1} = cA; %#ok<AGROW>
        cand_c{end+1} = cB; %#ok<AGROW>

        cand_g = [cand_g; gA; gB]; %#ok<AGROW>
        cand_e = [cand_e eA eB]; %#ok<AGROW>

        cand_d2 = [cand_d2; sum((cA-c0).^2); sum((cB-c0).^2)]; %#ok<AGROW>

    end

    if DIAG.fail(in)
        continue
    end

    nq = numel(cand_g);

    A = cand_e;
    g = cand_g(:);

    % Phase-sum constraints:
    %   sum_q in phase alpha x_q = p_alpha
    S = zeros(numel(active),nq);

    for ia = 1:numel(active)
        iph = active(ia);
        S(ia,cand_phase == iph) = 1;
    end

    bS = pn(active);

    % --------------------------------------------------------
    % Solve local QP
    % --------------------------------------------------------
    if use_quadprog == 1 && exist('quadprog','file') == 2

        H = eta_n*(A.'*A);
        H = 0.5*(H + H.') + qp_reg*eye(nq);

        f = g - eta_n*(A.'*En) + 0.5*prox_weight*cand_d2;

        lb = zeros(nq,1);
        ub = inf(nq,1);

        try
            qopt = optimoptions('quadprog','Display','off');
            [x,~,exitflag] = quadprog(H,f,[],[],S,bS,lb,ub,[],qopt);
            solver_name = "quadprog";
        catch
            [x,exitflag] = Linprog_Slack_Local(A,g,cand_d2,S,bS,En,prox_weight,slack_weight);
            solver_name = "linprog_slack";
        end

        if isempty(x) || exitflag <= 0
            [x,exitflag] = Linprog_Slack_Local(A,g,cand_d2,S,bS,En,prox_weight,slack_weight);
            solver_name = "linprog_slack";
        end

    else

        [x,exitflag] = Linprog_Slack_Local(A,g,cand_d2,S,bS,En,prox_weight,slack_weight);
        solver_name = "linprog_slack";

    end

    if isempty(x) || exitflag <= 0

        DIAG.fail(in)     = true;
        DIAG.exitflag(in) = exitflag;
        DIAG.solver(in)   = solver_name;
        continue

    end

    x(abs(x) < 1e-14) = 0;

    % Small renormalization per phase
    for ia = 1:numel(active)
        iph = active(ia);
        idq = find(cand_phase == iph);
        sx  = sum(x(idq));

        if sx > 0
            x(idq) = x(idq) * pn(iph)/sx;
        end
    end

    % --------------------------------------------------------
    % Effective composition per phase
    % --------------------------------------------------------
    for ia = 1:numel(active)

        iph = active(ia);
        idq = find(cand_phase == iph);

        piph = pn(iph);

        if piph <= p_cut
            continue
        end

        c_eff = zeros(size(c_cols{iph}(:,in)));

        for iq = idq(:).'
            c_eff = c_eff + x(iq)*cand_c{iq};
        end

        c_eff = c_eff./piph;

        c0 = c_cols{iph}(:,in);

        % Optional relaxation/cap for safety
        dc_eff = c_eff - c0;

        if isfinite(dc_eff_cap)
            md = max(abs(dc_eff));
            if md > dc_eff_cap
                dc_eff = dc_eff * dc_eff_cap/md;
            end
        end

        c_eff = c0 + relax_c*dc_eff;

        c_new_cols{iph}(:,in) = c_eff;

    end

    % --------------------------------------------------------
    % Penalty chemical potential
    % sign convention:
    %   mu = eta*(E_target - E_model)
    % --------------------------------------------------------
    Efit = A*x;
    mu_e(:,in) = eta_n*(En - Efit);

    DIAG.mass_err(:,in) = Efit - En;
    DIAG.exitflag(in)   = exitflag;
    DIAG.nvar(in)       = nq;
    DIAG.solver(in)     = solver_name;

end

% ------------------------------------------------------------
% Return c in original orientation
% ------------------------------------------------------------
c_new = c_old;

for iph = 1:Nphase
    c_new{iph} = Cols_To_C_Local(c_new_cols{iph},c_was_t(iph));
end

% ------------------------------------------------------------
% Stable chi prototype
% ------------------------------------------------------------
for ie = 1:Ne
    for je = 1:Ne
        if ie == je
            chi{ie,je}(:) = chi_floor;
        else
            chi{ie,je}(:) = 0;
        end
    end
end

DIAG.max_mass_err = max(abs(DIAG.mass_err),[],1);

end


% =========================================================================
% Orient c to nc x Nnode
% =========================================================================
function [C,was_t] = C_To_Cols_Local(Cin,Nnode)

was_t = false;

if size(Cin,2) == Nnode

    C = Cin;

elseif size(Cin,1) == Nnode

    C = Cin.';
    was_t = true;

else

    error('QP_Calculator: cannot orient c array.')

end

end


% =========================================================================
% Return c to original orientation
% =========================================================================
function Cout = Cols_To_C_Local(C,was_t)

if was_t
    Cout = C.';
else
    Cout = C;
end

end


% =========================================================================
% Try evaluating PhaseThermo
% =========================================================================
function [ok,g,e,R] = Try_Eval_Local(par,c,Ne)

ok = false;
g  = inf;
e  = nan(Ne,1);
R  = [];

try

    R = PhaseThermo(par,c);

    g = Extract_G_Local(R);
    e = Extract_E_Local(R,Ne);

    g = g(1);
    e = e(:,1);

    if isfinite(g) && all(isfinite(e))
        ok = true;
    end

catch

    ok = false;

end

end


% =========================================================================
% Extract Gibbs/free energy
% =========================================================================
function g = Extract_G_Local(R)

if isfield(R,'g')
    g = R.g;
elseif isfield(R,'G')
    g = R.G;
elseif isfield(R,'f')
    g = R.f;
elseif isfield(R,'F')
    g = R.F;
else
    error('QP_Calculator: PhaseThermo output has no g/G/f/F field.')
end

g = g(:);

end


% =========================================================================
% Extract elemental composition
% =========================================================================
function e = Extract_E_Local(R,Ne)

if isfield(R,'e')
    e = R.e;
elseif isfield(R,'E')
    e = R.E;
elseif isfield(R,'N')
    e = R.N;
elseif isfield(R,'comp')
    e = R.comp;
else
    error('QP_Calculator: PhaseThermo output has no e/E/N/comp field.')
end

if size(e,1) ~= Ne && size(e,2) == Ne
    e = e.';
end

if size(e,1) > Ne
    e = e(1:Ne,:);
end

if size(e,1) ~= Ne
    error('QP_Calculator: extracted e has wrong number of rows.')
end

end


% =========================================================================
% Split direction
% =========================================================================
function v = Split_Direction_Local(par,c0,opts,use_neg_hessian_dir)

nc = numel(c0);
v  = zeros(nc,1);
v(1) = 1;

% User-given split direction
if isfield(opts,'split_dir') && ~isempty(opts.split_dir)

    if iscell(opts.split_dir)
        % If caller passes cell, the caller should already select correct
        % phase outside.  Here use first nonempty direction as prototype.
        for i = 1:numel(opts.split_dir)
            if ~isempty(opts.split_dir{i})
                vv = opts.split_dir{i}(:);
                if numel(vv) == nc
                    v = vv;
                    return
                end
            end
        end
    else
        vv = opts.split_dir(:);
        if numel(vv) == nc
            v = vv;
            return
        end
    end

end

% User-given split component index
if isfield(opts,'split_component') && ~isempty(opts.split_component)

    ic = opts.split_component;

    if ic >= 1 && ic <= nc
        v = zeros(nc,1);
        v(ic) = 1;
        return
    end

end

% Try most negative Hessian direction
if use_neg_hessian_dir == 1

    try

        R = PhaseThermo(par,c0);

        H = [];

        if isfield(R,'H_c')
            H = R.H_c;
        elseif isfield(R,'Hc')
            H = R.Hc;
        elseif isfield(R,'H')
            H = R.H;
        elseif isfield(R,'h')
            H = R.h;
        end

        if ~isempty(H)

            H = 0.5*(H + H.');
            [V,D] = eig(H);
            lam = diag(D);

            [lmin,imin] = min(lam);

            if isfinite(lmin)
                v = real(V(:,imin));
                return
            end

        end

    catch

        % Keep default direction

    end
end

end


% =========================================================================
% Linprog fallback with L1 mass slack
% =========================================================================
function [x,exitflag] = Linprog_Slack_Local(A,g,cand_d2,S,bS,E,prox_weight,slack_weight)

Ne = size(A,1);
nq = size(A,2);

% Variables:
%   y = [x ; rplus ; rminus]
%
% Mass slack:
%   A*x + rplus - rminus = E

f = [g(:) + 0.5*prox_weight*cand_d2(:); ...
     slack_weight*ones(2*Ne,1)];

Aeq1 = [A, eye(Ne), -eye(Ne)];
beq1 = E;

Aeq2 = [S, zeros(size(S,1),2*Ne)];
beq2 = bS;

Aeq = [Aeq1; Aeq2];
beq = [beq1; beq2];

lb = zeros(nq + 2*Ne,1);
ub = inf(nq + 2*Ne,1);

try

    lopt = optimoptions('linprog','Display','off');
    [y,~,exitflag] = linprog(f,[],[],Aeq,beq,lb,ub,lopt);

catch

    [y,~,exitflag] = linprog(f,[],[],Aeq,beq,lb,ub);

end

if isempty(y)
    x = [];
else
    x = y(1:nq);
end

end