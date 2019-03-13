% Faber approximation for ellipse-shaped domains(recursive polynomial
% computation)
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
H = zeros(N, N);
for n = 1:(N - 1)
    H(n + 1, n + 1) = H(n, n) + (1 - 0.1 * (n - 3)) * 2 * pi * 1e13 * hbar;
end

% dipole moment
d = 1e-29;

% Hamiltonian (off-diagonal elements)
for n = 1:(N-1)
    H(n, n + 1) = d * E;
    H(n + 1, n) = d * E;
end

% initial rho
rho_init = zeros(N, N);
rho_init(1, 1) = 1;

% enable higher precisionf
digits(32);
H = vpa(H);

% right-hand side Liouville-von Neumann equation-->Liouvillian
rhs = @(rho) ( -1i/hbar * (H * rho - rho * H));

% time span
% te = 0.1e-12;
te = 0.1e-12;
dt = 1e-16;
t = 0:dt:te;

% state variables
rho_me = rho_init;
rho_fab = rho_init;

% result: trace error
trace_me = zeros(size(t));
trace_fab = zeros(size(t));

% result: populations
pop_me = zeros(N, length(t));
pop_fab = zeros(N, length(t));



% Compute real part of the field of values
Re_rhs = (rhs(rho_init)+rhs(rho_init)')/2;

% Compute imaginary part of the field of values 
Im_rhs = (rhs(rho_init)-rhs(rho_init)')/(2i);

% Compute factors for rectangle from field of values enclosing spectrum
% and factors of optimal ellipse enclosing the rectangle 
% [Xi_1, Xi_2, -l, l] containing the spectrum of L
Xi_1 = min(eig(Re_rhs));
Xi_2 = max(eig(Re_rhs));
l = max(imag(eig(Im_rhs)));
c = abs((Xi_2-Xi_1)/2);

% Compute optimal scaling factor to preserve stability
sf = ((l^(2/3)+c^(2/3))^(3/2))/2;

% scaled right-hand side Liouville-von Neumann equation-->scaled Liouvillian
rhs_sc = @(rho) ( -1i/hbar * (H * rho - rho * H)/sf);

% Scaled factors of rectangle and optimal ellipse
l_sc = l/sf;
Xi_1_sc = Xi_1/sf;
Xi_2_sc = Xi_2/sf;
c_sc = abs((Xi_2_sc-Xi_1_sc)/2);

% Compute the center of the optimal ellipse, first laurent expansion
% coefficient
b_0 =(Xi_1_sc+Xi_2_sc)/2;

% Compute the logarithmic capacity of the optimal ellipse
r = (sqrt(c_sc^2+(l_sc*c_sc^2)^(2/3))+sqrt(l_sc^2+(c_sc*l_sc^2)^(2/3)))/2;

% Compute second factor in finite conformal mapping for ellipsoid domain
b_1 = (c_sc^(2/3)+l_sc^(2/3)*(c_sc^(4/3)-l_sc^(4/3)))/(4*r);

% Scaled timestep
dt_tilde = sf*dt;

% scaled time vector
t_tilde = 0:dt_tilde:te;

% prepare matrix exponential
U = expm(-1i * dt_tilde/hbar * H);

% establish polynomial truncation order
for n = 1:1000
    c_m = @(dt_tilde, m) ((-1i/sqrt(b_1))^m*exp(dt_tilde*b_0)*...
    besselj(m, 2*dt_tilde*sqrt(-b_1)));
    CM(n) = c_m(dt_tilde, n);
    if abs(CM(n))<10e-15
        break
    end
end

% polynomial truncation order
M = length(CM);

% non recursive solution
I = ones(size(rhs_sc(rho_fab)));

P = zeros((M+1)*N,N);

tic;
% length of time vector is equivalent to number of timesteps
for n = 1:10
    rho_me = U * rho_me * U';
    
    % Faber coefficients dependent on discrete time step
    CM = fab(M, n*dt_tilde, b_0, b_1);
    
    % First compute the polynomial projection P
    % P_0 = F_0*rho_fab (from previous timestep n-1)
    P(1:N, 1:N) = rho_fab;
    
    % P_1 = F1*rho_fab (from previous timestep n-1)
    P(N+1:2*N, 1:N) = (rhs_sc(rho_fab)-b_0*I)*rho_fab;
    
    % P_2 ... P_M = F_2*rho_fab...P_M*rho_fab(from previous timestep n-1)
    for i =2:(M+1)
        P((i*N)+1:(i+1)*N, 1:N) = rhs_sc(rho_fab)*P((i-1)*N+1:i*N, 1:N)-...
            b_0*P((i-1)*N+1:i*N, 1:N)-b_1*P((i-1-1)*N+1:(i-1)*N, 1:N);
    end
    
    % z_0
    rho_fab(1:N, 1:N) = CM(1)*P(1:N, 1:N);
    % z_m+1 = z_m 
    for i = 2:(M+1)
        rho_fab = rho_fab + CM(i)*P((i-1)*N+1:i*N, 1:N);
    end
    
    trace_me(n) = trace(rho_me) - 1;
    pop_me(:, n) = diag(rho_me);

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
plot(t/1e-12, real(pop_me(3, :)), '-.', 'Color', [0, 101, 189]/255, 'DisplayName', 'ME');
grid on;
hold on;
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

% tau = discrete time step
function CM = fab(M, tau, b_0, b_1)
% initialize CM
% Faber coefficients for an elliptic domain
c_m = @(dt_tilde, m) ((-1i/sqrt(b_1))^m*exp(dt_tilde*b_0)*...
     besselj(m, 2*dt_tilde*sqrt(-b_1)));

CM = zeros(M+1,1);

for m = 1:(M+1)
    CM(m) = c_m(tau, m);
end
end
