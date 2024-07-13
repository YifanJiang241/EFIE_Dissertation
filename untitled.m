% Define constants
Eamp = 1.0;
Epsilon_0 = 8.854e-12;
Epsilon_d = 2.47912e-11;
Mu_0 = 12.56637061e-7;
c = 1.0 / sqrt(Mu_0 * Epsilon_0);
GrossStep = 10.0;
f = 970e6;
Lambda = c / f;
DeltaX = Lambda / 4.0;
Omega = 2.0 * pi * f;
j = 1i;
GrossNoSteps = 70;
Beta_0 = Omega * sqrt(Mu_0 * Epsilon_0);
Eta_0 = sqrt(Mu_0 / Epsilon_0);
TOL = 10e-15;
NoLinesubs = floor((GrossStep * GrossNoSteps) / DeltaX);
Xsource = 0.0;
Ysource = 442.0;
I = 1.0;
PI = 3.14159265358979323846;
EXP = 2.718281828;
Const = 5;  % New constant

% Load terrain profile data
data = load('X.04');
X = data(:, 1);
Y = data(:, 2);

% Initialize arrays
ModJ = zeros(1, NoLinesubs);
ModEt = zeros(1, NoLinesubs);
J = zeros(1, NoLinesubs);
Et = zeros(1, NoLinesubs);
Sigma = zeros(1, NoLinesubs);  % New Sigma array

% Calculate solution for surface current
% (Note: The following code assumes the functions R_source_p, R_source_obs,
% R_p_q, R_surf_obs, abs, EiRad, H02, Z, Zself, Exp, cplx_Exp have been defined.)

% Forward scattering
J(1) = EiRad(R_source_p(1, DeltaX, GrossStep, Ysource, Xsource, Y), 1, Beta_0, Omega, Epsilon_0) / Zself(1, DeltaX, GrossStep, Y, Beta_0, Epsilon_0, Omega, PI, EXP);

for p = 1:NoLinesubs %
    SUM = 0;
    for q = 1:p-1
        SUM = SUM + R_p_q(q, q+1, DeltaX, GrossStep, Y) * Z(p, q, Beta_0, Epsilon_0, DeltaX, GrossStep, Y, Omega) * J(q);
    end
    J(p) = (EiRad(R_source_p(p, DeltaX, GrossStep, Ysource, Xsource, Y), p, Beta_0, Omega, Epsilon_0) - SUM) / Zself(p, DeltaX, GrossStep, Y, Beta_0, Epsilon_0, Omega, PI, EXP);
    Sigma(p) = Const * abs(J(p)) * abs(H02(Beta_0 * R_source_p(p, DeltaX, GrossStep, Ysource, Xsource, Y)));  % Calculate Sigma
end

% Backscattering
for p = NoLinesubs:-1:2
    SUM = 0;
    for q = NoLinesubs:-1:p+1
        SUM = SUM + R_p_q(q, q+1, DeltaX, GrossStep, Y) * Z(p, q, Beta_0, Epsilon_0, DeltaX, GrossStep, Y, Omega) * J(q);
    end
    J(p) = J(p) + (-1.0 * SUM) / Zself(p, DeltaX, GrossStep, Y, Beta_0, Epsilon_0, Omega, PI, EXP);
end

% Write current values to file
writematrix([DeltaX*(0:NoLinesubs-1); abs(J)], 'J.dat')

figure;
plot(0:DeltaX:DeltaX*(NoLinesubs-1), abs(J));
xlabel('Index');
ylabel('Magnitude of J');
title('Surface Current Density Magnitude');
grid on;

figure;
plot(0:DeltaX:DeltaX*(NoLinesubs-1), Sigma);
xlabel('Index');
ylabel('Sigma');
title('Sigma Calculation');
grid on;

% Calculate total electric field above the surface
coutput = fopen('E.dat', 'w');
for index = 1:NoLinesubs
    Et(index) = 0;
    for n = 1:index
        Et(index) = Et(index) + (J(n) * R_p_q(n, n+1, DeltaX, GrossStep, Y) * Z_1(R_surf_obs(n, index, DeltaX, GrossStep, Y), Beta_0, Omega, Epsilon_0));
    end
    Et(index) = EiRad(R_source_obs(index, DeltaX, GrossStep, Y, Xsource, Ysource), index, Beta_0, Omega, Epsilon_0) - Et(index);
    fprintf(coutput, '%f  %f\n', DeltaX * (index - 1), 20.0 * log10(abs(Et(index)) / sqrt(R_source_obs(index, DeltaX, GrossStep, Y, Xsource, Ysource))));
end
fclose(coutput);

Eplot = load('E.dat');
figure;
Ex = Eplot(:, 1);
Ey = Eplot(:, 2);
plot(Ex, Ey);
xlabel('Distance (meters)');
ylabel('Electric Field Magnitude (dB)');
title('Electric Field Above the Surface');
grid on;

% Function definitions start
function result = x(a, DeltaX)    
    result = double(a) * DeltaX;
end

function result = R_source_p(p, DeltaX, GrossStep, Ysource, Xsource, Y)
    result = sqrt(((Xsource - x(p, DeltaX))^2) + ((Ysource - y(p, DeltaX, GrossStep, Y))^2));
end

function result = y(a, DeltaX, GrossStep, Y)
    Temp = (a * DeltaX) / GrossStep;
    Index = floor(Temp);
    Prop = Temp - Index;
    s = Y(Index+1) + (Prop * (Y(Index+2) - Y(Index+1)));
    result = s;
end

function result = R_source_obs(p, DeltaX, GrossStep, Y, Xsource, Ysource)
    result = sqrt(((Xsource - x(p, DeltaX))^2) + ((Ysource - y(p, DeltaX, GrossStep, Y) - 2.4)^2));
end

function result = R_p_q(p, q, DeltaX, GrossStep, Y)
    result = sqrt(((x(q, DeltaX) - x(p, DeltaX))^2) + ((y(q, DeltaX, GrossStep, Y) - y(p, DeltaX, GrossStep, Y))^2));
end

function result = R_surf_obs(p, q, DeltaX, GrossStep, Y)
    result = sqrt(((x(q, DeltaX) - x(p, DeltaX))^2) + ((y(q, DeltaX, GrossStep, Y) + 2.4 - y(p, DeltaX, GrossStep, Y))^2));
end

function result = EiRad(dist, ~, Beta_0, Omega, Epsilon_0)
    E = -((Beta_0^2) / (4.0 * Omega * Epsilon_0)) * H02(Beta_0 * dist);
    result = E;
end

function result = H02(Arg)
    H = complex(besselj(0, Arg), -bessely(0, Arg));
    result = H;
end

function result = Z(p, q, Beta_0, Epsilon_0, DeltaX, GrossStep, Y, Omega)
    H = ((Beta_0^2) / (4.0 * Omega * Epsilon_0)) * H02(Beta_0 * R_p_q(p, q, DeltaX, GrossStep, Y));
    result = H;
end

function result = Z_1(R, Beta_0, Omega, Epsilon_0)
    H = ((Beta_0^2) / (4.0 * Omega * Epsilon_0)) * H02(Beta_0 * R);
    result = H;
end

function result = Zself(i, DeltaX, GrossStep, Y, Beta_0, Epsilon_0, Omega, PI, EXP)
    Linesubln = R_p_q(i, i+1, DeltaX, GrossStep, Y);
    H = complex(((Beta_0^2) / (4.0 * Omega * Epsilon_0)) * Linesubln, -((Beta_0^2) / (4.0 * Omega * Epsilon_0)) * ((2.0 * Linesubln) / PI) * log((1.781 * Beta_0 * Linesubln) / (4.0 * EXP)));
    result = H;
end

function result = Exp(d)
    result = exp(d);
end

function result = cplx_Exp(d)
    result = complex(cos(d),-sin(d));
end