% Faber approximation for ellipse-shaped domains
clear;
close all;

% constants
e0 = 1.60217646e-19;
hbar = 1.05457168e-34;

% level count
N = 6;

% static electric field (V/m)
E = 9e9;

% Hamiltonian (diagonal elements)
Hamiltonian = zeros(N, N);
for n = 1:(N - 1)
    Hamiltonian(n + 1, n + 1) = Hamiltonian(n, n) + (1 - 0.1 * (n - 3)) * 2 * pi * 1e13 * hbar;
end

% dipole moment
d = 1e-29;

% Hamiltonian (off-diagonal elements)
for n = 1:(N-1)
    Hamiltonian(n, n + 1) = d * E;
    Hamiltonian(n + 1, n) = d * E;
end

% initial rho
rho_init = zeros(N,N);
rho_init(1, 1) = 1;

% enable higher precisionf
digits(32);
Hamiltonian = vpa(Hamiltonian);

% right-hand side Liouville-von Neumann equation-->Liouvillian
Liouvillian = @(rho) ( -1i/hbar * (Hamiltonian * rho - rho * Hamiltonian));

% time span
% te = 0.1e-12;
te = 0.1e-12;
dt = 1e-16;
t = 0:dt:te;

% state variables
rho_me = rho_init;
setglobal(rho_init);

% result: trace error
trace_me = zeros(size(t));
trace_fab = zeros(size(t));

% result: populations
pop_me = zeros(N, length(t));
pop_fab = zeros(N, length(t));

% prepare matrix exponential
U = expm(-1i * dt/hbar * Hamiltonian);

tic;
% length of time vector is equivalent to number of timesteps
for n = 1:100
    rho_me = U * rho_me * U';
    rho_fab = getglobal
    % Liouvillian changing with every iteration
    L = Liouvillian(rho_fab);
    
    % result of Faber approximation
    rho_fab= faber(L, dt, N);
    
    % evaluate the trace at every iteration
    trace_me(n) = trace(rho_me) - 1;
    pop_me(:, n) = diag(rho_me);
    
    % evaluate the trace at every iteration
    trace_fab(n) = trace(rho_fab) - 1;
    pop_fab(:, n) = diag(rho_fab);
end
toc;

% display minimum population values
disp(['Minimum population values ME: ' num2str(min(min(real(pop_me))))]);
disp(['Minimum population values FABER: ' num2str(min(min(real(pop_fab))))]);

% plot population rho_33
papersize = [ 15 12 ];
fig = figure('units', 'centimeters');
pos = get(gcf, 'pos');
set(gcf, 'pos', [pos(1) pos(2) papersize]);
% plot third line of populations
plot(t/1e-12, real(pop_me(3, :)), '-.', 'Color', [0, 101, 189]/255, 'DisplayName', 'ME');
grid on;
hold on;
% plot third line of populations
plot(t/1e-12, real(pop_fab(3, :)), '-', 'Color', [227, 114, 34]/255, 'DisplayName', 'FABER');
ax = gca;
set(gca, 'FontName', 'Helvetica', 'FontSize', 12);
xlabel('Time/ps');
legend('show', 'Location', 'northeast');
ylabel('Population \rho_{33}/1');
ylim([-0.1 1]);

% create inset
axes(fig, 'Position', [0.25, 0.72, 0.4, 0.2]);
box on;
plot(t/1e-12, real(pop_me(3, :)), '-.', 'LineWidth', 1, 'Color', [0, 101, 189]/255, 'DisplayName', 'ME');
grid on;
hold on;
plot(t/1e-12, real(pop_fab(3, :)), '-', 'LineWidth', 1, 'Color', [227, 114, 34]/255, 'DisplayName', 'FABER');
xlim([0.07 0.1]);
ylim([-0.1 0.6]);

set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', 'PaperSize', papersize);
print(fig, 'fab.pdf', '-dpdf', '-fillpage');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Faber approximation of the matrix exponential
function rho_fab = faber(L, dt, N)
% Compute real part of the field of values
Re_L = (L+L')/2;

% Compute imaginary part of the field of values 
Im_L = (L-L')/(2i);

% Compute factors for rectangle from field of values enclosing spectrum
% and factors of ellipse enclosing the rectangle optimally
% [Xi_1, Xi_2, -l, l] containing the spectrum of L
Xi_1 = min(eig(Re_L));
Xi_2 = max(eig(Re_L));
l = max(eig(Im_L));
c = abs((Xi_2-Xi_1)/2);

% Compute optimal scaling factor to preserve stability
sf = ((l^(2/3)+c^(2/3))^(3/2))/2;

% scaled Liouvillian
L_sc = L/sf;

% Scaled factors of rectangle and optimal ellipse
l_sc = l/sf;
Xi_1_sc = Xi_1/sf;
Xi_2_sc = Xi_2/sf;
c_sc = abs((Xi_2_sc-Xi_1_sc)/2);

% Compute the center of the optimal ellipse, first laurent expansion
% coefficient of PSI
b_0 =(Xi_1_sc+Xi_2_sc)/2;

% Compute the logarithmic capacity of the optimal ellipse
r = (sqrt(c_sc^2+(l_sc*c_sc^2)^(2/3))+sqrt(l_sc^2+(c_sc*l_sc^2)^(2/3)))/2;

% Compute second factor in finite conformal mapping for elliptic domain
b_1 = ((c_sc^(2/3)+l_sc^(2/3))*(c_sc^(4/3)-l_sc^(4/3)))/(4*r);

% Scaled timestep
dt_tilde = sf*dt;

% Faber coefficients for ellipic domain
    c_m = @(dt_tilde, m) (((-1i/sqrt(b_1))^m)*exp(dt_tilde*b_0)*...
    besselj(m, 2*dt_tilde*sqrt(-b_1)));

% establish polynomial truncation order so that M >e*sf*dt
M = 1;
while true
    CM(M) = double(c_m(dt_tilde, M));
    if abs(CM(M))<10e-15
        break
    else
    M = M+1;
    end
end

% Compute time dependent Faber coefficients of order M for elliptic domain
% initialize CM
CM = zeros(M+1,1);
for m = 0:M
    CM(m+1) = c_m(dt_tilde, m);
end

% Compute matrix valued polynomials with initial value density matrix from
% previous iteration
I = ones(size(L_sc));
P = zeros((M+1)*N, N);

% Compute matrix valued Faber polynomial recurrence relation
    % Compute initial value polynomials P_0, P_1, P_2:
    % P_0 = F_0*rho_fab (from previous timestep n-1)
    P(1:N, 1:N) = I;
    
    % P_1 = F1*rho_fab (from previous timestep n-1)
    P((N+1):(2*N), 1:N) = (L_sc-b_0*I);
    
    % P_2 = F_2*rho_fab (from previous timestep n-1)
    P((2*N+1):(3*N), 1:N) = (L_sc-b_0*I)-2*b_1*I;
    rho_fab = getglobal
    % c0*p_0+c1*p_1+c_2*p_2
    setglobal(CM(1)*P(1:N, 1:N)+CM(2)*P(N+1:2*N, 1:N)+CM(3)*P(2*N+1:3*N, 1:N)); 
    
    % P_3 ... P_M = F_2*rho_fab...P_M*rho_fab(from previous timestep n-1)
    for i = 3:M
        % p_{m+1} = L_sc(rho_fab)*p_m - b_0*p_m - b_1*p_{m-1}
        P((i*N)+1:(i+1)*N, 1:N) = L_sc*P((i-1)*N+1:(i)*N, 1:N)...
            -b_0*P((i-1)*N+1:(i)*N, 1:N)...
            -b_1*P((i-1-1)*N+1:(i-1)*N, 1:N);
        
        rho_fab = rho_fab+CM(i+1)*P((i*N)+1:(i+1)*N, 1:N);
    end
end

function setglobal(rho_fab)
    global rho_fab_persistent
    rho_fab_persistent = rho_fab
end

function result = getglobal
    global rho_fab_persistent
    result = rho_fab_persistent
end
