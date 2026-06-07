clear; clf

% ------------------------------------------------------------
% 1D Cahn-Hilliard spinodal test
% ------------------------------------------------------------
N      = 256;
Lbox   = 2*pi;
dx     = Lbox/N;
x      = (0:N-1)'*dx;

M      = 1.0;
kappa  = 1e-2;

dt     = 1e-4;
nstep  = 0.4/dt;

rng(1)
c0     = 0.01*randn(N,1);

% ------------------------------------------------------------
% Periodic Laplacian
% ------------------------------------------------------------
e      = ones(N,1);
L      = spdiags([e -2*e e],[-1 0 1],N,N)/dx^2;
L(1,N) = 1/dx^2;
L(N,1) = 1/dx^2;

I      = speye(N);

% ------------------------------------------------------------
% Initial states
% ------------------------------------------------------------
c_raw  = c0;
c_cs   = c0;

E_raw  = zeros(nstep,1);
E_cs   = zeros(nstep,1);
amp_raw= zeros(nstep,1);
amp_cs = zeros(nstep,1);

for it = 1:nstep

    % ========================================================
    % A. Raw semi-implicit tangent
    % ========================================================
    mu_old = c_raw.^3 - c_raw - kappa*(L*c_raw);
    Hraw   = 3*c_raw.^2 - 1;

    % Linearized:
    %   mu_new = mu_old + Hraw*dc - kappa*L*dc
    %
    % CH:
    %   dc/dt = M*L*mu_new
    %
    % Therefore:
    %   [I - dt*M*L*diag(Hraw) + dt*M*kappa*L*L] dc
    %       = dt*M*L*mu_old
    Araw   = I - dt*M*L*spdiags(Hraw,0,N,N) + dt*M*kappa*(L*L);
    rhs    = dt*M*L*mu_old;
    dc_raw = Araw\rhs;
    c_raw  = c_raw + dc_raw;

    % ========================================================
    % B. Convex-split / stabilized semi-implicit tangent
    % ========================================================
    mu_old = c_cs.^3 - c_cs - kappa*(L*c_cs);
    Hraw   = 3*c_cs.^2 - 1;

    % Add a positive shift so Hpos is positive everywhere.
    Hfloor = 1e-2;
    S      = max(0,-min(Hraw)) + Hfloor;
    Hpos   = Hraw + S;

    % Same equation, but use Hpos only in the implicit tangent.
    % Raw mu_old is unchanged and still contains the spinodal force.
    Acs = I - dt*M*L*spdiags(Hpos,0,N,N) + dt*M*kappa*(L*L);
    rhs = dt*M*L*mu_old;

    dc_cs = Acs\rhs;
    c_cs  = c_cs + dc_cs;

    % --------------------------------------------------------
    % Diagnostics
    % --------------------------------------------------------
    E_raw(it) = FreeEnergy(c_raw,L,kappa,dx);
    E_cs(it)  = FreeEnergy(c_cs,L,kappa,dx);

    amp_raw(it) = max(abs(c_raw));
    amp_cs(it)  = max(abs(c_cs));

    if any(~isfinite(c_raw)) || max(abs(c_raw)) > 10
        fprintf('Raw method blew up at step %d\n',it)
        E_raw(it:end)   = NaN;
        amp_raw(it:end) = NaN;
        break
    end
end

% ------------------------------------------------------------
% Plot final profiles
% ------------------------------------------------------------
subplot(3,1,1)
plot(x,c0,'k--'); hold on
plot(x,c_raw,'r-')
plot(x,c_cs,'b-')
legend('initial','raw tangent','convex split')
xlabel('x')
ylabel('c')
title('Final composition profile')

% ------------------------------------------------------------
% Plot free energy
% ------------------------------------------------------------
subplot(3,1,2)
plot(E_raw,'r-'); hold on
plot(E_cs,'b-')
legend('raw tangent','convex split')
xlabel('step')
ylabel('free energy')
title('Free energy evolution')

% ------------------------------------------------------------
% Plot max amplitude
% ------------------------------------------------------------
subplot(3,1,3)
plot(amp_raw,'r-'); hold on
plot(amp_cs,'b-')
legend('raw tangent','convex split')
xlabel('step')
ylabel('max |c|')
title('Amplitude growth')

% ------------------------------------------------------------
% Helper
% ------------------------------------------------------------
function F = FreeEnergy(c,L,kappa,dx)

fbulk = 0.25*(c.^2 - 1).^2;

cx2 = -c.*(L*c);   % approximate |grad c|^2 after integration by parts

F = sum(fbulk + 0.5*kappa*cx2)*dx;

end