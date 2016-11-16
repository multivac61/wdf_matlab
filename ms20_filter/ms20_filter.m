% WDF model of a MS20 principal filter, op amp modelled as nullator
% clear all; close all;
close all;

Fs = 96e3; % sample rate (Hz)
N = Fs; % number of samples to simulate
input = [1, zeros(1, N-1)];

% w = 0.5e-3*blackmanharris(64);
% input = [w' input];
%[input, Fs] = audioread('funk.wav');
%input = (input(:,1)+input(:,2))/2;
%output = zeros(1, N);

% Component values
R1 = 2.2e3;
R2 = 47;
R3 = 800e3;
R4 = 800e3;
R5 = 47;
R6 = 5.038e3;
% R6 = 7.5e3;
R7 = 10e3;
R8p = 0.1e3;
R8m = 9.9e3;
RL = 10e3;
C1 = 15e-9;
C2 = 15e-9;

%% Port resistances
% R-type adaptor port resistances
% Ra = R3;
Ra = R3;
Rb = R4;
Rc = 1 / (2 * Fs * C2);
Rd = R6;
Re = R7;
Rf = R8p;
Rg = R8m;
% Rh = 1 / (2 * Fs * C1);
syms Rh;

%% Calculate Scattering matrix
N_dep_sources = 1;
N_nodes = 13;
N_ports = 8;

Ga = 1/Ra;
Gb = 1/Rb;
Gc = 1/Rc;
Gd = 1/Rd;
Ge = 1/Re;
Gf = 1/Rf;
Gg = 1/Rg;
Gh = 1/Rh;

Y = sym(zeros(N_nodes,N_nodes));
Y(1,1) = Y(1,1) + Ga;
Y(1,2) = Y(1,2) - Ga;
Y(2,1) = Y(2,1) - Ga;
Y(2,2) = Y(2,2) + Ga;

Y(3,3) = Y(3,3) + Gb;
Y(3,4) = Y(3,4) - Gb;
Y(4,3) = Y(4,3) - Gb;
Y(4,4) = Y(4,4) + Gb;

Y(5,5) = Y(5,5) + Gc;

Y(7,7) = Y(7,7) + Gd;

Y(9,9) = Y(9,9) + Ge;
Y(9,8) = Y(9,8) - Ge;
Y(8,9) = Y(8,9) - Ge;
Y(8,8) = Y(8,8) + Ge;

Y(11,11) = Y(11,11) + Gf;
Y(11,10) = Y(11,10) - Gf;
Y(10,11) = Y(10,11) - Gf;
Y(10,10) = Y(10,10) + Gf;

Y(11,11) = Y(11,11) + Gg;
Y(11,12) = Y(11,12) - Gg;
Y(12,11) = Y(12,11) - Gg;
Y(12,12) = Y(12,12) + Gg;

Y(13,13) = Y(13,13) + Gh;
Y(13,2) = Y(13,2) - Gh;
Y(2,13) = Y(2,13) - Gh;
Y(2,2) = Y(2,2) + Gh;

B = sym(zeros(N_ports+N_dep_sources, N_nodes));
B(1,1) = 1;

B(2,3) = 1;
B(2,2) = -1;

B(3,5) = 1;
B(3,4) = -1;

B(4,7) = 1;
B(4,6) = -1;

B(5,8) = 1;
B(5,6) = -1;

B(6,10) = 1;
B(6,9) = -1;

B(7,12) = 1;

B(8,13) = 1;
B(8,11) = -1;

A = transpose(B);

% Add nullator models
B(9,6) = -1;
B(9,4) = 1;

A(9,9) = 1;

D = sym(zeros(N_ports+N_dep_sources,N_ports+N_dep_sources));

X = [Y, A; B, D];

%% Port resitances
Rp = diag([Ra, Rb, Rc, Rd, Re, Rf, Rg, Rh]);

%% Scatterig matrix
S = eye(N_ports) + 2*[zeros(N_ports,N_nodes) Rp zeros(N_ports,N_dep_sources)] * inv(X) * [zeros(N_ports,N_nodes) eye(N_ports) zeros(N_ports,N_dep_sources)]';

eqn = S(8,8) == 0;
Rh_sol = solve(eqn, Rh);
R_adp = double(vpa(Rh_sol));
S = subs(S, Rh, R_adp);

S = double(S);

%% initialize wave variables
Rp13 = R_adp;
Rp12 = 1 / (2 * Fs * C1);
Rp11 = 1 / (1/Rp13 + 1/Rp12);

% Parallel adaptor 1
ap11 = 0; bp11 = 0;
ap12 = 0; bp12 = 0;
ap13 = 0; bp13 = 0;

% R-type adaptor
br1 = 0; br2 = 0; br3 = 0; br4 = 0; br5 = 0; br6 = 0; br7 = 0; br8 = 0;
ar1 = 0; ar2 = 0; ar3 = 0; ar4 = 0; ar5 = 0; ar6 = 0; ar7 = 0; ar8 = 0;

% Delay lines
delay_C1 = 0; delay_C2 = 0;

% tmp
output_vin = zeros(1,length(input));

% Diode-specific values
Is = 2.52e-9;     % Saturation current.
Vt = 0.02585;     % Thermal voltage.
n  = 1.792;       % Quality factor
R0 = Rp11;

% Lookup table for nonlinearity

if exist('lookupTable') == 0
    load('lookupTable.mat')
end

% Output variables
N = length(input);
output = zeros(1, N);

%% The simulation loop:
for i = 1:N % run each time sample until N
    if mod(i, N / 10) == 0
        fprintf('i = %f\n', i); 
    end
    ar1 = input(i); ar2 = 0; ar3 = delay_C2; ar4 = 0;
    ar5 = 0; ar6 = 0; ar7 = 0; ar8 = 0;
    
    tmp = S * [ar1; ar2; ar3; ar4; ar5; ar6; ar7; ar8];
    br1 = tmp(1);
    br2 = tmp(2);
    br3 = tmp(3);
    br4 = tmp(4);
    br5 = tmp(5);
    br6 = tmp(6);
    br7 = tmp(7);
    br8 = tmp(8);
    
    ap13 = br8; ap12 = delay_C1;
    [bp11, bp12, bp13] = parallel3(ap11, Rp11, ap12, Rp12, ap13, Rp13);
    
    % Root - AP Diodes
    a = bp11;
%     a_str = num2str(a);
%     ap11 = sign(a) * (abs(a) + 2*R0*Is - 2*Vt*n*lambertw(R0*Is / (Vt*n) * exp((R0*Is + abs(a)) / (Vt*n))));
    ap11 = bp11; % open circuit
%     if lookupTable.isKey(a_str)
%         ap11 = lookupTable(a_str);
%     else
%         ap11 = sign(a) * (abs(a) + 2*R0*Is - 2*Vt*n*lambertw(R0*Is / (Vt*n) * exp((R0*Is + abs(a)) / (Vt*n))));
%         lookupTable(a_str) = ap11;
%     end
    [bp11, bp12, bp13] = parallel3(ap11, Rp11, ap12, Rp12, ap13, Rp13);
    
    ar8 = bp13;
    tmp = S * [ar1; ar2; ar3; ar4; ar5; ar6; ar7; ar8];
    br1 = tmp(1);
    br2 = tmp(2);
    br3 = tmp(3);
    br4 = tmp(4);
    br5 = tmp(5);
    br6 = tmp(6);
    br7 = tmp(7);
    br8 = tmp(8);
    
    output(i) = 0.5 * (ar7 + br7) - 0.5 * (ar6 + br6);
    delay_C1 = bp12; delay_C2 = br3;
end; 

%% Plots
% Plot impulse- and magnitude response.
f = figure();
subplot(2,1,1);
plot(output);
title('Impulse response');

subplot(2,1,2);
[H, w] = freqz(output, 1, N);
semilogx(w / (2*pi) * Fs, 20 * log10(abs(H)));
title('Magnitude response');
axis tight;

% Plot output, before and after write to file.
norm_output = output / (max(output) - min(output));
filename = 'output/wdf_ms20_filter_64bitprecision_96k_32bit.wav'
audiowrite(filename, output, Fs,'BitsPerSample',32)

figure();
subplot(3,1,1);
spectrogram(output,blackmanharris(1024),10,256,Fs,'yaxis');
title('Output');

subplot(3,1,2);
spectrogram(norm_output,blackmanharris(1024),10,256,Fs,'yaxis');
title('Normalized output');

subplot(3,1,3);
[my_wav, my_fs] = audioread(filename);
spectrogram(my_wav,blackmanharris(1024),10,256,my_fs,'yaxis');
title('Re-read .wav output');