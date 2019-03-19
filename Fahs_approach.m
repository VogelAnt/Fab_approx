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
rho_fab = rho_init;

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
for n = 1:length(t)
    
    % Matrix exponential via Diagonalization
    rho_me = U * rho_me * U';
    % evaluate the trace at every iteration
    trace_me(n) = trace(rho_me) - 1;
    pop_me(:, n) = diag(rho_me);
end
toc;

tic;
H = -1i*Hamiltonian/hbar;

% Compute real part of the field of values
Re_H = (H+H')/2;

% Compute imaginary part of the field of values 
Im_H = (H-H')/(2i);

% Compute factors for rectangle from field of values enclosing spectrum
% and factors of ellipse enclosing the rectangle optimally
% [Xi_1, Xi_2, -l, l] containing the spectrum of L
Xi_1 = min(eig(Re_H));
Xi_2 = max(eig(Re_H));
l = max(eig(Im_H));
c = abs((Xi_2-Xi_1)/2);

% Compute optimal scaling factor to preserve stability
sf = ((l^(2/3)+c^(2/3))^(3/2))/2;

% scaled Hamiltonian
H_sc = H/sf;

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

U_fab = faber(b_0, b_1, H_sc, dt_tilde, N);

% length of time vector is equivalent to number of timesteps
for n = 1:length(t)
  
    % result of Faber approximation
    rho_fab= U_fab*rho_fab*U_fab';
    
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
% plot third line of populations -.
plot(t/1e-12, real(pop_me(3, :)), 'o', 'Color', [0, 101, 189]/255, 'DisplayName', 'ME');
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
function z = faber(b0, b1, H_sc, dt_tilde, N)
    % Assign time step and laurent expansion coefficient values
    tau = dt_tilde;
    b_0 = b0;
    b_1 = b1;
    
    % Faber coefficients for ellipic domain
    c_m = @(tau, m) (((-1i/sqrt(b_1))^m)*exp(tau*b_0)*...
    besselj(m, 2*tau*sqrt(-b_1)));

    % establish polynomial truncation order so that M >e*sf*dt
    M = 1;
    while true
        orderiterator(M) = double(c_m(tau, M));  
%        if M>(exp(1)*tau)|| 
        if abs(orderiterator(M))<10e-15
            break
        else
        M = M+1;
        end
    end

    % Compute time dependent Faber coefficients of order M for elliptic domain
    % initialize CM
    CM = zeros(M+1,1);
    for m = 0:M
        CM(m+1) = c_m(tau, m);
    end

    % Compute matrix valued polynomials with initial value density matrix from
    % previous iteration
    I = eye(N,N);
    P = zeros((M+1)*N, N);

    % Compute matrix valued Faber polynomial recurrence relation
    % Compute initial value polynomials P_0, P_1, P_2:
    % F_0
    P(1:N, 1:N) = I;
    
    % F_1*rho_n
    P((N+1):(2*N), 1:N) = H_sc-b_0*I;
    
    % F_2*rho_n = (L_sc-b0)*F_1-2*b1*F_0
    P((2*N+1):(3*N), 1:N) = (H_sc-b_0*I)*P((N+1):(2*N), 1:N)...
        -2*b_1*P(1:N, 1:N);
    
    if M < 2
    temp = CM(1)*P(1:N, 1:N)+CM(2)*P(N+1:2*N, 1:N);
    z = temp;
    else
        
    % z = c0*F_0+c1*F_1+c_2*F_2
    temp = CM(1)*P(1:N, 1:N)+CM(2)*P(N+1:2*N, 1:N)+CM(3)*P(2*N+1:3*N, 1:N); 
    
    % F_3*rho_n ... F_M*rho_n
    for i = 3:M
        % F_{m+1} 
        P((i*N)+1:(i+1)*N, 1:N) = H_sc*P((i-1)*N+1:(i)*N, 1:N)...
            -b_0*P((i-1)*N+1:(i)*N, 1:N)...
            -b_1*P((i-1-1)*N+1:(i-1)*N, 1:N);
        % c3*F_3...cm*F_m
        temp = temp + CM(i+1)*P((i*N)+1:(i+1)*N, 1:N);
    end
    z = temp;
    end
end
