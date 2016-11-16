% WDF model of a MS20 principal filter, op amp modelled as nullator
% Backwards-Euler method used for digitization
clear all; close all;

% Fs = 22.05e3; % sample rate (Hz)
% Fs = 44.1e3; % sample rate (Hz)
Fs = 96e3; % sample rate (Hz)
N = 10*Fs; % number of samples to simulate

input = [1, zeros(1, N-1)];

% Component values
R1 = 500;
R2 = 10e6;
RL = 10e3;
Rv = 100;
C1 = 1e-9;
C2 = 1e-9;

%% Port resistances
% R-type adaptor port resistances
% syms Ra;
Ra = Rv;
Rb = 1 / (Fs * C1);
Rc = R1;
Rd = R2;
Re = 1 / (Fs * C2);
Rf = RL;

%% Calculate Scattering matrix
N_dep_sources = 1;
N_nodes = 10;
N_ports = 6;

Ga = 1/Ra;
Gb = 1/Rb;
Gc = 1/Rc;
Gd = 1/Rd;
Ge = 1/Re;
Gf = 1/Rf;

Y = sym(zeros(N_nodes,N_nodes));
Y(5,5) = Y(5,5) + Ga;
Y(1,5) = Y(1,5) - Ga;
Y(5,1) = Y(5,1) - Ga;
Y(1,1) = Y(1,1) + Ga;

Y(6,6) = Y(6,6) + Gb;
Y(3,6) = Y(3,6) - Gb;
Y(6,3) = Y(6,3) - Gb;
Y(3,3) = Y(3,3) + Gb;

Y(7,7) = Y(7,7) + Gc;
Y(3,7) = Y(3,7) - Gc;
Y(7,3) = Y(7,3) - Gc;
Y(3,3) = Y(3,3) + Gc;

Y(8,8) = Y(8,8) + Gd;
Y(2,8) = Y(2,8) - Gd;
Y(8,2) = Y(8,2) - Gd;
Y(2,2) = Y(2,2) + Gd;

Y(9,9) = Y(9,9) + Ge;
Y(3,9) = Y(3,9) - Ge;
Y(9,3) = Y(9,3) - Ge;
Y(3,3) = Y(3,3) + Ge;

Y(10,10) = Y(10,10) + Gf;
Y(4,10)  = Y(4,10) - Gf;
Y(10,4)  = Y(10,4) - Gf;
Y(4,4)   = Y(4,4) + Gf;

B = sym(zeros(N_ports+N_dep_sources, N_nodes));
B(1,5) = 1;
% B(1,1) = -1;

B(2,6) = 1;
B(2,2) = -1;

B(3,7) = 1;
% B(3,1) = -1;

B(4,8) = 1;
B(4,4) = -1;

B(5,9) = 1;
B(5,4) = -1;

B(6,10) = 1;
% B(6,1) = -1;

A = transpose(B);

% Add nullator models
B(7,1) = 1;
B(7,2) = -1;

A(4,7) = 1;

D = sym(zeros(N_ports+N_dep_sources,N_ports+N_dep_sources));

X = [Y, A; B, D];

%% Port resitances
Rp = diag([Ra, Rb, Rc, Rd, Re, Rf]);

%% Scatterig matrix
S = eye(N_ports) + 2*[zeros(N_ports,N_nodes) Rp zeros(N_ports,N_dep_sources)] * inv(X) * [zeros(N_ports,N_nodes) eye(N_ports) zeros(N_ports,N_dep_sources)]';

% eqn = S(1,1) == 0;
% Ra_sol = solve(eqn, Ra);
% R_adp = double(vpa(Ra_sol))
% S = subs(S, Ra, R_adp);

S = double(S);

%% initialize wave variables
% R-type adaptor
br1 = 0; br2 = 0; br3 = 0; br4 = 0; br5 = 0; br6 = 0;
ar1 = 0; ar2 = 0; ar3 = 0; ar4 = 0; ar5 = 0; ar6 = 0;

% Delay lines
delay_C1a = 0; delay_C2a = 0;
delay_C1b = 0; delay_C2b = 0;

% tmp
N = length(input);
output_vin = double(zeros(1,N));

%% The simulation loop:
for i = 1:N % run each time sample until N
    if mod(i, 10000) == 0
        disp(N/i + '%');
    end
    
    ar1 = input(i); ar2 = 0.5*(delay_C1a + delay_C1b); ar3 = 0; 
    ar4 = 0; ar5 = 0.5*(delay_C2a + delay_C2b); ar6 = 0;
    
    tmp = S * [ar1; ar2; ar3; ar4; ar5; ar6];
    br1 = tmp(1);
    br2 = tmp(2);
    br3 = tmp(3);
    br4 = tmp(4);
    br5 = tmp(5);
    br6 = tmp(6);

    output(i) = 0.5 * (ar6 + br6);
    
    delay_C1a = ar2; delay_C1b = br2; 
    delay_C2a = ar5; delay_C2b = br5;
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

% Normalize output before plot and save
norm_output = output / (max(output) - min(output));
filename = 'output/wdf_bridged_t_resonator_bweuler_64bitprecision_96k_32bit_10s.wav'
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