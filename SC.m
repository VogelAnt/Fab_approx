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

% Scaled timestep
dt_tilde = sf*dt;
U_fab = faber_poly(H_sc, dt_tilde, N, Xi_1_sc, Xi_2_sc, l_sc);

% length of time vector is equivalent to number of timesteps
for n = 1:length(t)
    
    % Matrix exponential via Diagonalization
    rho_me = U * rho_me * U';
    
    % result of Faber approximation
    rho_fab= U_fab*rho_fab*U_fab';
    
    % evaluate the trace at every iteration
    trace_fab(n) = trace(rho_fab) - 1;
    pop_fab(:, n) = diag(rho_fab);    
    
    % evaluate the trace at every iteration
    trace_me(n) = trace(rho_me) - 1;
    pop_me(:, n) = diag(rho_me);
end

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
plot(t/1e-12, real(pop_me(3, :)), 'o', 'LineWidth', 1, 'Color', [0, 101, 189]/255, 'DisplayName', 'ME');
grid on;
hold on;
plot(t/1e-12, real(pop_fab(3, :)), '-', 'LineWidth', 1, 'Color', [227, 114, 34]/255, 'DisplayName', 'FABER');
xlim([0.07 0.1]);
ylim([-0.1 0.6]);

set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', 'PaperSize', papersize);
print(fig, 'fab.pdf', '-dpdf', '-fillpage');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Faber approximation of the matrix exponential with modification of
% SC-toolbox files (Copyright 1998 by Toby Driscoll)
function z = faber_poly(H_sc, dt_tilde, N, Xi_1_sc, Xi_2_sc, l_sc)
    %   Credit where credit is due:
    %   Copyright 1998 by Toby Driscoll.
    %   $Id: faber.m 298 2009-09-15 14:36:37Z driscoll $
    %   This function follows somewhat closely the procedure outlined in
    %   section 4 of Starke and Varga (Num. Math., 1993), except that no
    %   symmetry of the polygon is assumed.
    
    % create rectangular polygon from field of values with vertices in counter-clockwise order
    % compute map from exterior of unit disc to exterior of polygon
    PSI = extermap(polygon([Xi_1_sc+1i*l_sc; Xi_1_sc-1i*l_sc; Xi_2_sc-1i*l_sc; Xi_2_sc+1i*l_sc]));

    % extract polygon structure from exterior map PSI
    p = polygon(PSI);

    % vector containing interior angles of vertics
    alpha = angle(p);

    % vector containing exterior angles of vertices
    beta = 1-angle(p);

    % Returns structure of the Schwarz-Christoffel map parameters
    r = parameters(PSI);

    % prevertices
    z = r.prevertex;

    % logarithmic capacity of domain
    C = r.constant;

    % length of polygon (4 for our rectangular bounding box)
    n = length(p);

    % Computation of Faber coefficients
    syms u;
    CM = zeros(M+1, 1);
    I_handle = @(omega) (omega^(-2)*(((1-omega/z(1))^(1-alpha(1)))...
        *((1-omega/z(2))^(1-alpha(2)))*((1-omega/z(3))^(1-alpha(3)))...
        *((1-omega/z(4))^(1-alpha(4)))));

    % symbolically integrate inner function
    I_inner = int(I_handle, 0, u);
    f_0 = matlabFunction(I_inner);
    f_1 = exp(dt_tilde*f_0);
    f_2 = matlabFunction(f1);

    % convert result above to anonymous function
    for i=0:M
        temp = f_2/(u^(i+1));
        % Compute the integral across the polygon with Gauss-Kronrod Quadrature
        CM = 1/(2*pi)*(quadgk(temp, Xi_1_sc+1i*l_sc, Xi_2_sc-1i*l_sc,'waypoints', [Xi_1_sc-1i*l_sc])...
        +quadgk(temp, Xi_2_sc-1i*l_sc, Xi_2_sc+1i*l_sc));
    end

    % coefficients gamma_k of the binomial expansion
    gamma = ones(n,M+1);
    for k = 1:M
      gamma(:,k+1) = -gamma(:,k).*(beta-k+1)./(k*z);
    end

    % Compute the coefficients b_m of the Laurent expansion of PSI
    % initialize output vector
    e1 = zeros(M+1,1);
    e1(1) = 1;
    b_m = -c*e1;
    % Coefficients of Laurent expansion polynomials as a product of 
    % lower triangular Toeplitz matrices
    for j = 1:n
      b_m = toeplitz(gamma(j,:).', e1')*b_m;
    end
    b_m = b_m(3:M+1)./(-(1:M-1)');
    x0 = 10^(-10/M);

    % estimation of b_0 by evaluating f near the origin and substracting off
    % the known part of the series
    b_0 = eval(PSI,x0) + C/x0 - x0.^(1:M-1)*b_m;

    % laurent expansion coefficients of PSI
    b_m = [b_0;b_m];

    % Polynomial recurrence relation
    P = zeros((M+1)*N, N);
    I = zeros(size(H_sc));

    % F_0 = I
    P(1:N, 1:N) = rho_n;

    % F_1*rho_n = (L_sc -b0)*rho_n
    P(N+1:2*N, 1:N) = H_sc-b_m(1)*rho_n;

    % F_2*rho_n = L_sc*F_1*rho_n-2*b1*rho_n-c0*F_1*rho_n
    P(2*N+1:3*N, 1:N) = H_sc*P(N+1:2*N, 1:N)-2*b_m(2)*rho_n-b(1)*P(N+1:2*N, 1:N);

    if M < 2
    temp = CM(1)*P(1:N, 1:N)+CM(2)*P(N+1:2*N, 1:N);
    z = temp;
    else

    %
    temp = CM(1)*P(1:N, 1:N)+CM(2)*P(N+1:2*N, 1:N)+CM(3)*P(2*N+1:3*N, 1:N);

    % initialize sum variable
    sigma = zeros(N,N);

    % matrix valued Faber polynomial recurrence relation
    % F_3...F_M
    for i=3:M

        % compute sum with shifted indeces
        for j=0:i
            if (i-j)== 0
            sigma = sigma + b_m(j+1)*P(1:N, 1:N);
            else
            sigma = sigma + b_m(j+1)*P((i-j)*N+1:(i-j+1)*N, 1:N);
            end
        end
        P(((i*N)+1):(i+1)*N, 1:N) = L_sc*P(((i-1)*N+1):i*N, 1:N)...
            -(i-1)*b_m(i)...
            -sigma;

        % c3*F_3...cm*F_m
        temp = temp + CM(i+1)*P((i*N)+1:(i+1)*N, 1:N);

    end
    z = temp;
    end
end