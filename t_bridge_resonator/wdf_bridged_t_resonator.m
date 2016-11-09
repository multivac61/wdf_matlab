% WDF model of a MS20 principal filter, op amp modelled as nullator
% clear all; close all;
clear all; close all;

Fs = 10*96e3; % sample rate (Hz)
Fs = 655350;
N = 10*Fs; % number of samples to simulate
s = Fs;
input = zeros(1,N-1);
% input(1) = 10e-3;
% input(s/2) = 10e-3;
% input(s) = 10e-3;
% % input() = 10e-3;
% input(2*s) = 10e-3;

% d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.3,0.4,0.5,200);
% Hd = design(d,'equiripple');
% input = filter(Hd,input);
% input = input2;
% w = 10e-3*window_2;
% % w = 10e-3*ones(16, 1);
w = 0.5e-3*blackmanharris(64);
input = [w' input];
d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.01,0.05,0.5,200);
Hd = design(d,'equiripple');
input = filter(Hd,input);
% input = [w', input(1:s/2), w', input(s/2:s), w', input(s:3*s/2), w', input(s:end)];

%[input, Fs] = audioread('funk.wav');
%input = (input(:,1)+input(:,2))/2;
N = length(input);
%output = zeros(1, N);

% Component values
R1 = 500;
R2 = 10e6;
RL = 10e3;
RC = 10;
Rv = 100;
C1 = 30e-9;
C2 = 30e-9;

%% Port resistances
% Series port resistances
Rs12 = 1 / (2 * Fs * C1);
Rs13 = RC;
Rs11 = Rs12 + Rs13;

Rs22 = 1 / (2 * Fs * C2);
Rs23 = RC;
Rs21 = Rs22 + Rs23;

% R-type adaptor port resistances
% syms Ra;
Ra = Rv;
Rb = Rs11;
Rc = R1;
Rd = R2;
Re = Rs21;
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
Y(4,10) = Y(4,10) - Gf;
Y(10,4) = Y(10,4) - Gf;
Y(4,4) = Y(4,4) + Gf;


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

S = double(vpa(S));

%% initialize wave variables
% Series adaptor
as11 = 0; bs11 = 0;
as12 = 0; bs12 = 0;
as13 = 0; bs13 = 0;

% Series adaptor
as21 = 0; bs21 = 0;
as22 = 0; bs22 = 0;
as23 = 0; bs23 = 0;

% R-type adaptor
br1 = 0; br2 = 0; br3 = 0; br4 = 0; br5 = 0; br6 = 0;
ar1 = 0; ar2 = 0; ar3 = 0; ar4 = 0; ar5 = 0; ar6 = 0;

% Delay lines
delay_C1 = 0; delay_C2 = 0;

% tmp
output_vin = zeros(1,length(input));

%% The simulation loop:
for i = 1:N % run each time sample until N
    as11 = 0; as12 = delay_C1; as13 = 0;
    [bs11, bs12, bs13] = series3(as11, Rs11, as12, Rs12, as13, Rs13);

    as21 = 0; as22 = delay_C2; as23 = 0;
    [bs21, bs22, bs23] = series3(as21, Rs21, as22, Rs22, as23, Rs23);
    
    ar1 = input(i); ar2 = bs11; ar3 = 0; 
    ar4 = 0; ar5 = bs21; ar6 = 0;
    
    tmp = S * [ar1; ar2; ar3; ar4; ar5; ar6];
    br1 = tmp(1);
    br2 = tmp(2);
    br3 = tmp(3);
    br4 = tmp(4);
    br5 = tmp(5);
    br6 = tmp(6);

    as11 = br2;
    [bs11, bs12, bs13] = series3(as11, Rs11, as12, Rs12, as13, Rs13);

    as21 = br5;
    [bs21, bs22, bs23] = series3(as21, Rs21, as22, Rs22, as23, Rs23);
    
%     output(i) = 0.5 * (ar6 + br6);
    output(i) = 0.5 * (ar3 + br3) + 0.5 * (ar5 + br5);
    delay_C1 = bs12; delay_C2 = bs22;
end; 

% plot the magnitude response
f = figure();
[H, w] = freqz(output, 1, N);
semilogx(w / (2*pi) * Fs, 20 * log10(abs(H)));
% title('output R3p');
f = w / (2*pi) * Fs;

% output = filter(Hd,output);
% g = figure();
% plot(output);
% 
% h = figure();
% plot(input);
% % f = w / (2*pi) * Fs;
% % save('/home/horigome/McGill/mumt618/wdf/wdfplot/wdf_sallen_key.mat', 'H', 'f');
audiowrite('resonator_output_75_vpa_bl_64_Rc.flac', output, Fs, 'BitsPerSample',24);
