% WDF model of a MS20 principal filter, op amp modelled as nullator
% clear all; close all;
clear all; close all;

% Fs = 10*96e3; % sample rate (Hz)
% Fs = 655350;
Fs = 16e3;
N = 5*Fs; % number of samples to simulate
s = Fs;
input = zeros(1,N-1);
input(1) = 10e-3;
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
% w = 0.5e-3*blackmanharris(64);
% input = [w' input];
% d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.01,0.05,0.5,200);
% Hd = design(d,'equiripple');
% input = filter(Hd,input);
% input = [w', input(1:s/2), w', input(s/2:s), w', input(s:3*s/2), w', input(s:end)];

%[input, Fs] = audioread('funk.wav');
%input = (input(:,1)+input(:,2))/2;
N = length(input);
%output = zeros(1, N);

% Component values
R1 = 500;
R2 = 10e6;
RL = 10e3;
C1 = 1e-9;
C2 = 1e-9;

% Macromodel values - NJM2904D
Ri_cm = 5e6;
Ri_d = 3e6;
Ro = 75;
Ci_cm = 2e-12;
Ci_d = 1.4e-12;
Ao_cm = 5.623;
CMRR = 85;
Ao_d = Ao_cm * 10^(CMRR/20);
Ib1 = 27.5e-9;
Ib2 = 22.5e-9;
Ib1 = 0;
Ib2 = 0;
Voff = 2e-3;
Voff = 0;

Rbw = 100e3;
Cbw = 0.2653e-6;

%% Port resistances
% Parallel port resistances
Rpb23 = 2 * 1 / (2 * Fs * Ci_cm);
Rpb22 = 2*Ri_cm;
Rpb21 = 1/(1/Rpb22 + 1/Rpb23);

Rpc13 = 2 * 1 / (2 * Fs * Ci_cm);
Rpc12 = 2*Ri_cm;
Rpc11 = 1/(1/Rpc13 + 1/Rpc12);

Rpd13 = 1 / (2 * Fs * Ci_d);
Rpd12 = Ri_d;
Rpd11 = 1/(1/Rpd13 + 1/Rpd12);

% R-type adaptor port resistances
% syms Ra;
Ra = R1;
syms Rb;
Rc = Rpc11;
Rd = Rpd11;
Re = Rbw;
Rf = 1 / (2 * Fs * Cbw);
Rg = Ro;
Rh = RL;
Ri = 1 / (2 * Fs * C2);
Rj = R2;
Rk = 1 / (2 * Fs * C1);


%% Calculate Scattering matrix
N_dep_sources = 4;
N_nodes = 21;
N_ports = 11;

Ga = 1/Ra;
Gb = 1/Rb;
Gc = 1/Rc;
Gd = 1/Rd;
Ge = 1/Re;
Gf = 1/Rf;
Gg = 1/Rg;
Gh = 1/Rh;
Gi = 1/Ri;
Gj = 1/Rj;
Gk = 1/Rk;

Y = sym(zeros(N_nodes,N_nodes));
Y(2,2) = Y(2,2) + Ga;
Y(11,2) = Y(11,2) - Ga;
Y(2,11) = Y(2,11) - Ga;
Y(11,11) = Y(11,11) + Ga;

Y(12,12) = Y(12,12) + Gb;
Y(3,12) = Y(3,12) - Gb;
Y(12,3) = Y(12,3) - Gb;
Y(3,3) = Y(3,3) + Gb;

Y(4,4) = Y(4,4) + Gc;
Y(13,4) = Y(13,4) - Gc;
Y(4,13) = Y(4,13) - Gc;
Y(13,13) = Y(13,13) + Gc;

Y(14,14) = Y(14,14) + Gd;
Y(4,14) = Y(4,14) - Gd;
Y(14,4) = Y(14,4) - Gd;
Y(4,4) = Y(4,4) + Gd;

Y(15,15) = Y(15,15) + Ge;
Y(6,15) = Y(6,15) - Ge;
Y(15,6) = Y(15,6) - Ge;
Y(6,6) = Y(6,6) + Ge;

Y(1,1) = Y(1,1) + Gf;
Y(16,1) = Y(16,1) - Gf;
Y(1,16) = Y(1,16) - Gf;
Y(16,16) = Y(16,16) + Gf;

Y(17,17) = Y(17,17) + Gg;
Y(8,17) = Y(8,17) - Gg;
Y(17,8) = Y(17,8) - Gg;
Y(8,8) = Y(8,8) + Gg;

Y(18,18) = Y(18,18) + Gh;

Y(19,19) = Y(19,19) + Gi;
Y(18,19) = Y(18,19) - Gi;
Y(19,18) = Y(19,18) - Gi;
Y(18,18) = Y(18,18) + Gi;

Y(18,18) = Y(18,18) + Gj;
Y(20,18) = Y(20,18) - Gj;
Y(18,20) = Y(18,20) - Gj;
Y(20,20) = Y(20,20) + Gj;

Y(21,21) = Y(21,21) + Gk;
Y(4,21) = Y(4,21) - Gk;
Y(21,4) = Y(21,4) - Gk;
Y(4,4) = Y(4,4) + Gk;

B = sym(zeros(N_ports+N_dep_sources, N_nodes));
B(1,11) = 1;
B(1,1) = -1;

B(2,12) = 1;
B(2,1) = -1;

B(3,13) = 1;
B(3,1) = -1;

B(4,14) = 1;
B(4,3) = -1;

B(5,15) = 1;
B(5,5) = -1;

B(6,16) = 1;
B(6,5) = -1;

B(7,17) = 1;
B(7,7) = -1;

%B(8,16) = 1;
B(8,1) = -1;

B(9,19) = 1;
B(9,2) = -1;

B(10,20) = 1;
B(10,4) = -1;

B(11,21) = 1;
B(11,2) = -1;

% VCVS
B(12, 6) = 1;
B(12, 10) = -1;

B(13, 10) = 1;
B(13, 9) = -1;

B(14, 9) = 1;
B(14, 1) = -1;

B(15, 8) = 1;
B(15, 1) = -1;

A = transpose(B);

% VCVS
B(12,3) = -Ao_cm/2;
B(12,1) = Ao_cm/2;

B(13,3) = -Ao_d;
B(13,4) = Ao_d;

B(14,4) = -Ao_cm/2;
B(14,1) = B(14,1) + Ao_cm/2;

B(15,5) = -1;
B(15,1) = B(15,1) + 1;

D = sym(zeros(N_ports+N_dep_sources,N_ports+N_dep_sources));

X = [Y, A; B, D];

%% Port resitances
Rp = diag([Ra, Rb, Rc, Rd, Re, Rf, Rg, Rh, Ri, Rj, Rk]);

%% Scatterig matrix
S = eye(N_ports) + 2*[zeros(N_ports,N_nodes) Rp zeros(N_ports,N_dep_sources)] * inv(X) * [zeros(N_ports,N_nodes) eye(N_ports) zeros(N_ports,N_dep_sources)]';

eqn = S(2,2) == 0;
Rb_sol = solve(eqn, Rb);
R_adp = double(vpa(Rb_sol))
S = subs(S, Rb, R_adp);

S = double(vpa(S));

% Use adapted port
Rpb11 = R_adp;
Rpb13 = Rpb21;
Rpb12 = 1/(1/Rpb11 + 1/Rpb13);

% initialize wave variables
% Parallel adaptors for port B
apb11 = 0; bpb11 = 0;
apb12 = 0; bpb12 = 0;
apb13 = 0; bpb13 = 0;

apb21 = 0; bpb21 = 0;
apb22 = 0; bpb22 = 0;
apb23 = 0; bpb23 = 0;

% Parallel adaptors for port C
apb11 = 0; bpb11 = 0;
apb12 = 0; bpb12 = 0;
apb13 = 0; bpb13 = 0;

% Parallel adaptors for port D
apc11 = 0; bpc11 = 0;
apc12 = 0; bpc12 = 0;
apc13 = 0; bpc13 = 0;

% R-type adaptor
br1 = 0; br2 = 0; br3 = 0; br4 = 0; br5 = 0; br6 = 0; br7 = 0;  br8 = 0;  br9 = 0;  br10 = 0;  br11 = 0;
ar1 = 0; ar2 = 0; ar3 = 0; ar4 = 0; ar5 = 0; ar6 = 0; ar7 = 0;  ar8 = 0;  ar9 = 0;  ar10 = 0;  ar11 = 0;

% Delay lines
delay_C1 = 0; delay_C2 = 0;

delay_Ci_cm_b = 0; delay_Ci_cm_c = 0; 
delay_Ci_d = 0; delay_Cbw = 0;
output_vin = zeros(1,length(input));

% The simulation loop:
for i = 1:N % run each time sample until N
    % Wave up
    apc11 = 0; apc12 = Ib2*Rpc12; apc13 = delay_Ci_cm_c;
    [bpc11, bpc12, bpc13] = parallel3(apc11, Rpc11, apc12, Rpc12, apc13, Rpc13);

    apd11 = 0; apd12 = 0; apd13 = delay_Ci_d;
    [bpd11, bpd12, bpd13] = parallel3(apd11, Rpd11, apd12, Rpd12, apd13, Rpd13);

    ar1 = 0; ar2 = 0; ar3 = bpc11; ar4 = bpd11; ar5 = 0; ar6 = delay_Cbw;
    ar7 = 0; ar8 = 0; ar9 = delay_C2; ar10 = 0; ar11 = delay_C1; 
    
    tmp = S * [ar1; ar2; ar3; ar4; ar5; ar6; ar7; ar8; ar9; ar10; ar11];
    br1 = tmp(1);
    br2 = tmp(2);
    br3 = tmp(3);
    br4 = tmp(4);
    br5 = tmp(5);
    br6 = tmp(6);
    br7 = tmp(7);
    br8 = tmp(8);
    br9 = tmp(9);
    br10 = tmp(10);
    br11 = tmp(11);

    apb21 = 0; apb22 = Ib1*Rpb22; apb23 = delay_Ci_cm_b;
    [bpb21, bpb22, bpb23] = parallel3(apb21, Rpb21, apb22, Rpb22, apb23, Rpb23);

    apb11 = br2; apb12 = 0; apb13 = bpb21;
    [bpb11, bpb12, bpb13] = parallel3(apb11, Rpb11, apb12, Rpb12, apb13, Rpb13);

    % Root element
    apb12 = 2*(input(i) + Voff) - bpb12;

    % Wave-down
    [bpb11, bpb12, bpb13] = parallel3(apb11, Rpb11, apb12, Rpb12, apb13, Rpb13);

    apb21 = bpb13;
    [bpb21, bpb22, bpb23] = parallel3(apb21, Rpb21, apb22, Rpb22, apb23, Rpb23);

    ar2 = bpb11;
    tmp = S * [ar1; ar2; ar3; ar4; ar5; ar6; ar7; ar8; ar9; ar10; ar11];
    br1 = tmp(1);
    br2 = tmp(2);
    br3 = tmp(3);
    br4 = tmp(4);
    br5 = tmp(5);
    br6 = tmp(6);
    br7 = tmp(7);
    br8 = tmp(8);
    br9 = tmp(9);
    br10 = tmp(10);
    br11 = tmp(11);

    apc11 = br3;
    [bpc11, bpc12, bpc13] = parallel3(apc11, Rpc11, apc12, Rpc12, apc13, Rpc13);

    apd11 = br4;
    [bpd11, bpd12, bpd13] = parallel3(apd11, Rpd11, apd12, Rpd12, apd13, Rpd13);
    
%     output(i) = 0.5 * (apb12 + bpb12);
%     output(i) = 0.5*(ar8 + br8);
    output(i) = 0.5*(ar1 + br1) + 0.5*(ar9 + br9);

    delay_C1 = br11; delay_C2 = br9;
    delay_Ci_cm_b = bpb23; delay_Ci_cm_c = bpc13; delay_Ci_d = bpd13;
    delay_Cbw = br6;
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
% f = w / (2*pi) * Fs;
% save('/home/horigome/McGill/mumt618/wdf/wdfplot/wdf_sallen_key.mat', 'H', 'f');
audiowrite('resonator_output_75_vpa_bl_64_Rc.flac', output, Fs, 'BitsPerSample',24);
