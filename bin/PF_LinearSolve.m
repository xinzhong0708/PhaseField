function [sol,DIAG] = PF_LinearSolve(A,R,NUM,x0)
%PF_LINEARSOLVE Solve sparse linear system with selectable method.
%
% Usage:
%   [sol,LDIAG] = PF_LinearSolve(A,R,NUM);
%   [sol,LDIAG] = PF_LinearSolve(A,R,NUM,x0);
%
% NUM controls:
%   NUM.linear_solver       = 'direct';       % direct, direct_colamd,
%                                             % direct_symamd, direct_dissect,
%                                             % gmres_ilu, bicgstab_ilu, tfqmr_ilu
%   NUM.linear_perm         = 'symamd';       % iterative: symamd, colamd, none, dissect
%   NUM.linear_tol          = 1e-6;
%   NUM.linear_maxit        = 100;
%   NUM.gmres_restart       = 50;
%   NUM.ilu_type            = 'ilutp';        % ilutp, crout, nofill
%   NUM.ilu_droptol         = 1e-3;
%   NUM.ilu_thresh          = 0.1;
%   NUM.ilu_milu            = 'row';          % only for crout
%   NUM.ilu_reuse           = 0;              % reuse ILU factors across calls
%   NUM.ilu_reuse_steps     = 5;              % rebuild after this many uses
%   NUM.ilu_rebuild         = 1;              % rebuild once if cached ILU fails
%   NUM.linear_cache_id     = 'default';      % separate caches for ACCH/CHLE/etc.
%   NUM.ilu_cache_check_pattern = 1;          % also require same sparsity pattern
%   NUM.clear_ilu_cache     = 0;              % clear all persistent ILU caches
%   NUM.clear_ilu_cache_id  = '';             % clear one named cache only
%   NUM.direct_fallback     = 1;
%
% Notes:
%   - The returned sol always satisfies A*sol ~= R in the original ordering.
%   - ILU reuse is approximate. If the matrix changes too much, Krylov may
%     slow down or fail. With direct_fallback = 1, failed iterative solves
%     are replaced by direct LU.

persistent ILU_CACHE

if nargin < 3 || isempty(NUM)
    NUM = struct();
end
if nargin < 4
    x0 = [];
end

if isfield(NUM,'clear_ilu_cache') && NUM.clear_ilu_cache == 1
    ILU_CACHE = [];
elseif isfield(NUM,'clear_ilu_cache_id') && ~isempty(NUM.clear_ilu_cache_id)
    ILU_CACHE = ClearCacheID_Local(ILU_CACHE,NUM.clear_ilu_cache_id);
end

R = R(:);
N = size(A,1);

if size(A,2) ~= N || numel(R) ~= N
    error('PF_LinearSolve: A must be square and R must have matching length.')
end

solver = 'direct';
if isfield(NUM,'linear_solver') && ~isempty(NUM.linear_solver)
    solver = NUM.linear_solver;
end
solver = lower(solver);

% Defaults
linear_tol          = 1e-6;
linear_maxit        = 100;
gmres_restart       = 50;
ilu_type            = 'ilutp';
ilu_droptol         = 1e-3;
ilu_thresh          = 0.1;
ilu_milu            = 'row';
ilu_reuse           = 0;
ilu_reuse_steps     = 5;
ilu_rebuild         = 1;
linear_cache_id     = 'default';
ilu_cache_check_pattern = 1;
direct_fallback     = 1;

if isfield(NUM,'linear_tol')
    linear_tol = NUM.linear_tol;
end
if isfield(NUM,'linear_maxit')
    linear_maxit = NUM.linear_maxit;
end
if isfield(NUM,'gmres_restart')
    gmres_restart = NUM.gmres_restart;
end
if isfield(NUM,'ilu_type') && ~isempty(NUM.ilu_type)
    ilu_type = lower(NUM.ilu_type);
end
if isfield(NUM,'ilu_droptol')
    ilu_droptol = NUM.ilu_droptol;
end
if isfield(NUM,'ilu_thresh')
    ilu_thresh = NUM.ilu_thresh;
end
if isfield(NUM,'ilu_milu') && ~isempty(NUM.ilu_milu)
    ilu_milu = NUM.ilu_milu;
end
if isfield(NUM,'ilu_reuse')
    ilu_reuse = NUM.ilu_reuse;
end
if isfield(NUM,'ilu_reuse_steps')
    ilu_reuse_steps = NUM.ilu_reuse_steps;
end
if isfield(NUM,'ilu_rebuild')
    ilu_rebuild = NUM.ilu_rebuild;
end
if isfield(NUM,'linear_cache_id') && ~isempty(NUM.linear_cache_id)
    linear_cache_id = char(NUM.linear_cache_id);
end
if isfield(NUM,'ilu_cache_check_pattern')
    ilu_cache_check_pattern = NUM.ilu_cache_check_pattern;
end
if isfield(NUM,'direct_fallback')
    direct_fallback = NUM.direct_fallback;
end

if ilu_reuse_steps < 1
    ilu_reuse_steps = 1;
end

DIAG                      = struct();
DIAG.linear_solver        = solver;
DIAG.flag                 = 0;
DIAG.relres               = 0;
DIAG.iter                 = [0 0];
DIAG.fallback_used        = 0;
DIAG.perm_method          = '';
DIAG.prec_time            = 0;
DIAG.solve_time           = 0;
DIAG.total_time           = 0;
DIAG.matrix_size          = N;
DIAG.nnz                  = nnz(A);
DIAG.ilu_type             = ilu_type;
DIAG.ilu_droptol          = ilu_droptol;
DIAG.ilu_thresh           = ilu_thresh;
DIAG.ilu_reuse            = ilu_reuse;
DIAG.ilu_reuse_steps      = ilu_reuse_steps;
DIAG.linear_cache_id      = linear_cache_id;
DIAG.ilu_cache_check_pattern = ilu_cache_check_pattern;
DIAG.ilu_cache_hit        = 0;
DIAG.ilu_cache_age        = -1;
DIAG.ilu_rebuilt          = 0;
DIAG.ilu_rebuild          = 0;
DIAG.ilu_nnz_L            = 0;
DIAG.ilu_nnz_U            = 0;
DIAG.error_message        = '';

t_all = tic;

% -------------------------------------------------------------------------
% Direct sparse solve
% -------------------------------------------------------------------------
if strcmpi(solver,'direct') || strcmpi(solver,'direct_colamd') || ...
        strcmpi(solver,'direct_symamd') || strcmpi(solver,'direct_dissect') || ...
        strcmpi(solver,'direct_none')

    if strcmpi(solver,'direct_symamd')
        perm_method = 'symamd';
    elseif strcmpi(solver,'direct_dissect')
        perm_method = 'dissect';
    elseif strcmpi(solver,'direct_none')
        perm_method = 'none';
    else
        perm_method = 'colamd';
    end

    t_solve = tic;
    sol = DirectSolve_Local(A,R,perm_method);
    DIAG.solve_time = toc(t_solve);

    DIAG.perm_method = perm_method;
    DIAG.relres      = norm(A*sol - R)/max(norm(R),eps);
    DIAG.total_time  = toc(t_all);
    return
end

% -------------------------------------------------------------------------
% Iterative solve with optional ILU preconditioner
% -------------------------------------------------------------------------
perm_method = 'symamd';
if isfield(NUM,'linear_perm') && ~isempty(NUM.linear_perm)
    perm_method = NUM.linear_perm;
end
DIAG.perm_method = perm_method;

use_ilu = contains(solver,'_ilu');

% Use cached permutation and ILU factors if valid.
cache_ok  = 0;
cache_idx = 0;
pattern_sig = [];

if use_ilu && ilu_reuse == 1 && ilu_cache_check_pattern == 1
    pattern_sig = PatternSignature_Local(A);
end

if use_ilu && ilu_reuse == 1
    [cache_ok,cache_idx] = CacheOK_Local(ILU_CACHE,N,linear_cache_id,perm_method, ...
        ilu_type,ilu_droptol,ilu_thresh,ilu_milu,ilu_reuse_steps, ...
        pattern_sig,ilu_cache_check_pattern);
end

if cache_ok

    p    = ILU_CACHE(cache_idx).p;
    Lilu = ILU_CACHE(cache_idx).L;
    Uilu = ILU_CACHE(cache_idx).U;

    DIAG.ilu_cache_hit = 1;
    DIAG.ilu_cache_age = ILU_CACHE(cache_idx).age;
    DIAG.ilu_nnz_L     = nnz(Lilu);
    DIAG.ilu_nnz_U     = nnz(Uilu);

else

    p = Permutation_Local(A,perm_method);
    Lilu = [];
    Uilu = [];

end

Ap = A(p,p);
Rp = R(p);

if isempty(x0)
    x0p = [];
else
    x0 = x0(:);
    if numel(x0) ~= N
        error('PF_LinearSolve: x0 must have same length as R.')
    end
    x0p = x0(p);
end

try

    if use_ilu && ~cache_ok

        setup = ILUSetup_Local(ilu_type,ilu_droptol,ilu_thresh,ilu_milu);

        t_prec = tic;
        [Lilu,Uilu] = ilu(Ap,setup);
        DIAG.prec_time = toc(t_prec);

        DIAG.ilu_rebuilt = 1;
        DIAG.ilu_nnz_L   = nnz(Lilu);
        DIAG.ilu_nnz_U   = nnz(Uilu);

        if ilu_reuse == 1
            ILU_CACHE = StoreCache_Local(ILU_CACHE,N,linear_cache_id,p,Lilu,Uilu, ...
                perm_method,ilu_type,ilu_droptol,ilu_thresh,ilu_milu,pattern_sig);
        end

    end

    t_solve = tic;
    [yp,flag,relres,iter] = KrylovSolve_Local( ...
        solver,Ap,Rp,linear_tol,linear_maxit,gmres_restart,Lilu,Uilu,x0p,use_ilu);
    DIAG.solve_time = toc(t_solve);

catch ME

    flag   = 99;
    relres = inf;
    iter   = [0 0];
    yp     = zeros(size(Rp));
    DIAG.error_message = ME.message;

end

% If cached ILU failed, rebuild once with the current matrix and retry.
if use_ilu && ilu_reuse == 1 && DIAG.ilu_cache_hit == 1 && flag ~= 0 && ilu_rebuild == 1

    try

        setup = ILUSetup_Local(ilu_type,ilu_droptol,ilu_thresh,ilu_milu);

        t_prec = tic;
        [Lilu,Uilu] = ilu(Ap,setup);
        DIAG.prec_time = DIAG.prec_time + toc(t_prec);

        DIAG.ilu_rebuild         = 1;
        DIAG.ilu_rebuilt         = 1;
        DIAG.ilu_nnz_L           = nnz(Lilu);
        DIAG.ilu_nnz_U           = nnz(Uilu);

        ILU_CACHE = StoreCache_Local(ILU_CACHE,N,linear_cache_id,p,Lilu,Uilu, ...
            perm_method,ilu_type,ilu_droptol,ilu_thresh,ilu_milu,pattern_sig);

        t_solve = tic;
        [yp,flag,relres,iter] = KrylovSolve_Local( ...
            solver,Ap,Rp,linear_tol,linear_maxit,gmres_restart,Lilu,Uilu,x0p,use_ilu);
        DIAG.solve_time = DIAG.solve_time + toc(t_solve);

    catch ME

        flag = 98;
        relres = inf;
        iter = [0 0];
        yp = zeros(size(Rp));
        DIAG.error_message = ME.message;

    end

end

% Restore original ordering
sol = zeros(N,1);
sol(p) = yp;

% Fallback to direct solve if iterative method failed
if flag ~= 0 && direct_fallback == 1

    DIAG.fallback_used     = 1;
    DIAG.iterative_flag    = flag;
    DIAG.iterative_relres  = relres;
    DIAG.iterative_iter    = IterVec_Local(iter);

    t_solve = tic;
    sol = DirectSolve_Local(A,R,'colamd');
    DIAG.solve_time = DIAG.solve_time + toc(t_solve);

    DIAG.flag   = -abs(flag);
    DIAG.relres = norm(A*sol - R)/max(norm(R),eps);
    DIAG.iter   = [0 0];

else

    DIAG.flag   = flag;
    DIAG.relres = norm(A*sol - R)/max(norm(R),eps);
    DIAG.iter   = IterVec_Local(iter);

end

% Increase the age of only this named cache after one use.
if use_ilu && ilu_reuse == 1
    ILU_CACHE = IncreaseCacheAge_Local(ILU_CACHE,N,linear_cache_id,perm_method,ilu_type, ...
        ilu_droptol,ilu_thresh,ilu_milu,pattern_sig,ilu_cache_check_pattern);
end

DIAG.total_time = toc(t_all);

end


function sol = DirectSolve_Local(A,R,perm_method)
%DIRECTSOLVE_LOCAL Sparse direct solve with selected permutation.

N = size(A,1);

if strcmpi(perm_method,'none')

    sol = A\R;

elseif strcmpi(perm_method,'symamd')

    S = spones(A) + spones(A.');
    p = symamd(S);

    yp = A(p,p)\R(p);
    sol = zeros(N,1);
    sol(p) = yp;

elseif strcmpi(perm_method,'dissect')

    S = spones(A) + spones(A.');

    try
        p = dissect(S);
    catch
        p = symamd(S);
    end

    yp = A(p,p)\R(p);
    sol = zeros(N,1);
    sol(p) = yp;

else

    q = colamd(A);
    y = A(:,q)\R;

    sol = zeros(N,1);
    sol(q) = y;

end

end


function p = Permutation_Local(A,perm_method)
%PERMUTATION_LOCAL Permutation for iterative ILU/Krylov solve.

N = size(A,1);

if strcmpi(perm_method,'none')

    p = (1:N).';

elseif strcmpi(perm_method,'colamd')

    p = colamd(A).';

elseif strcmpi(perm_method,'dissect')

    S = spones(A) + spones(A.');

    try
        p = dissect(S).';
    catch
        p = symamd(S).';
    end

else

    S = spones(A) + spones(A.');
    p = symamd(S).';

end

end


function setup = ILUSetup_Local(ilu_type,ilu_droptol,ilu_thresh,ilu_milu)
%ILUSETUP_LOCAL Build MATLAB ilu setup struct.

setup = struct();
setup.type = ilu_type;

if strcmpi(ilu_type,'nofill')

    % No extra options.

elseif strcmpi(ilu_type,'crout')

    setup.droptol = ilu_droptol;
    setup.milu    = ilu_milu;

elseif strcmpi(ilu_type,'ilutp')

    setup.droptol = ilu_droptol;
    setup.thresh  = ilu_thresh;

else

    error('PF_LinearSolve: unknown NUM.ilu_type = %s',ilu_type)

end

end


function [x,flag,relres,iter] = KrylovSolve_Local( ...
    solver,A,b,tol,maxit,restart,L,U,x0,use_ilu)
%KRYLOVSOLVE_LOCAL Run selected Krylov method.

if strcmpi(solver,'gmres_ilu') || strcmpi(solver,'gmres')

    if use_ilu
        if isempty(x0)
            [x,flag,relres,iter] = gmres(A,b,restart,tol,maxit,L,U);
        else
            [x,flag,relres,iter] = gmres(A,b,restart,tol,maxit,L,U,x0);
        end
    else
        if isempty(x0)
            [x,flag,relres,iter] = gmres(A,b,restart,tol,maxit);
        else
            [x,flag,relres,iter] = gmres(A,b,restart,tol,maxit,[],[],x0);
        end
    end

elseif strcmpi(solver,'bicgstab_ilu') || strcmpi(solver,'bicgstab')

    if use_ilu
        if isempty(x0)
            [x,flag,relres,iter] = bicgstab(A,b,tol,maxit,L,U);
        else
            [x,flag,relres,iter] = bicgstab(A,b,tol,maxit,L,U,x0);
        end
    else
        if isempty(x0)
            [x,flag,relres,iter] = bicgstab(A,b,tol,maxit);
        else
            [x,flag,relres,iter] = bicgstab(A,b,tol,maxit,[],[],x0);
        end
    end

elseif strcmpi(solver,'tfqmr_ilu') || strcmpi(solver,'tfqmr')

    if use_ilu
        if isempty(x0)
            [x,flag,relres,iter] = tfqmr(A,b,tol,maxit,L,U);
        else
            [x,flag,relres,iter] = tfqmr(A,b,tol,maxit,L,U,x0);
        end
    else
        if isempty(x0)
            [x,flag,relres,iter] = tfqmr(A,b,tol,maxit);
        else
            [x,flag,relres,iter] = tfqmr(A,b,tol,maxit,[],[],x0);
        end
    end

else

    error('PF_LinearSolve: unknown NUM.linear_solver = %s',solver)

end

end



function [ok,idx] = CacheOK_Local(CACHE,N,cache_id,perm_method,ilu_type,ilu_droptol,ilu_thresh,ilu_milu,reuse_steps,pattern_sig,check_pattern)
%CACHEOK_LOCAL Check if cached ILU can be reused.

ok  = 0;
idx = 0;

if isempty(CACHE) || ~isstruct(CACHE)
    return
end

for ic = 1:numel(CACHE)

    C = CACHE(ic);

    if ~isfield(C,'N') || C.N ~= N
        continue
    end
    if ~isfield(C,'cache_id') || ~strcmp(C.cache_id,cache_id)
        continue
    end
    if ~isfield(C,'age') || C.age >= reuse_steps
        continue
    end
    if ~isfield(C,'perm_method') || ~strcmpi(C.perm_method,perm_method)
        continue
    end
    if ~isfield(C,'ilu_type') || ~strcmpi(C.ilu_type,ilu_type)
        continue
    end
    if ~isfield(C,'ilu_droptol') || C.ilu_droptol ~= ilu_droptol
        continue
    end
    if ~isfield(C,'ilu_thresh') || C.ilu_thresh ~= ilu_thresh
        continue
    end
    if ~isfield(C,'ilu_milu') || ~strcmpi(C.ilu_milu,ilu_milu)
        continue
    end
    if ~isfield(C,'L') || ~isfield(C,'U') || ~isfield(C,'p')
        continue
    end
    if check_pattern == 1
        if ~isfield(C,'pattern_sig') || isempty(C.pattern_sig) || isempty(pattern_sig)
            continue
        end
        if numel(C.pattern_sig) ~= numel(pattern_sig) || any(C.pattern_sig ~= pattern_sig)
            continue
        end
    end

    ok  = 1;
    idx = ic;
    return

end

end


function CACHE = StoreCache_Local(CACHE,N,cache_id,p,L,U,perm_method,ilu_type,ilu_droptol,ilu_thresh,ilu_milu,pattern_sig)
%STORECACHE_LOCAL Store or replace one named ILU cache.

C               = struct();
C.N             = N;
C.cache_id      = cache_id;
C.p             = p;
C.L             = L;
C.U             = U;
C.perm_method   = perm_method;
C.ilu_type      = ilu_type;
C.ilu_droptol   = ilu_droptol;
C.ilu_thresh    = ilu_thresh;
C.ilu_milu      = ilu_milu;
C.pattern_sig   = pattern_sig;
C.age           = 0;

if isempty(CACHE) || ~isstruct(CACHE)
    CACHE = C;
    return
end

for ic = 1:numel(CACHE)
    if isfield(CACHE(ic),'cache_id') && strcmp(CACHE(ic).cache_id,cache_id)
        CACHE(ic) = C;
        return
    end
end

CACHE(end+1) = C;

end


function CACHE = IncreaseCacheAge_Local(CACHE,N,cache_id,perm_method,ilu_type,ilu_droptol,ilu_thresh,ilu_milu,pattern_sig,check_pattern)
%INCREASECACHEAGE_LOCAL Increase age for the cache used by this call.

[ok,idx] = CacheOK_Local(CACHE,N,cache_id,perm_method,ilu_type,ilu_droptol,ilu_thresh,ilu_milu,inf,pattern_sig,check_pattern);

if ok
    CACHE(idx).age = CACHE(idx).age + 1;
end

end


function CACHE = ClearCacheID_Local(CACHE,cache_id)
%CLEARCACHEID_LOCAL Clear only one named ILU cache.

if isempty(CACHE) || ~isstruct(CACHE)
    return
end

keep = true(1,numel(CACHE));

for ic = 1:numel(CACHE)
    if isfield(CACHE(ic),'cache_id') && strcmp(CACHE(ic).cache_id,cache_id)
        keep(ic) = false;
    end
end

CACHE = CACHE(keep);

end


function sig = PatternSignature_Local(A)
%PATTERNSIGNATURE_LOCAL Cheap sparsity-pattern signature.
%
% This avoids reusing ILU factors between different systems that happen to
% have the same size.  It is much cheaper than an ILU factorization, but can
% be disabled with NUM.ilu_cache_check_pattern = 0.

[i,j] = find(spones(A));

i = double(i(:));
j = double(j(:));

sig = [size(A,1), nnz(A), sum(i), sum(j), sum(i.*j)];

end


function it = IterVec_Local(iter)
%ITERVEC_LOCAL Return iteration as 1x2 diagnostic vector.

if isempty(iter)
    it = [0 0];
elseif isscalar(iter)
    it = [iter 0];
else
    it = iter;
end

end
