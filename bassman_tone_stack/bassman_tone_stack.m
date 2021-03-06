% WDF version of the Fender Bassman
clear all; close all;

Fs = 96e3; % sample rate (Hz)
N = Fs; % number of samples to simulate

input = [1, zeros(1, N-1)];
output = zeros(1,length(input));

% Component values
R1p = 125e3;
R1m = 125e3;
R2 = 1000e3;
R3p = 25e3 * 0.5;
R3m = 25e3 * 0.5;
R4 = 56e3;
C1 = 0.25e-9;
C2 = 20e-9;
C3 = 20e-9;

%% Port resistances
Rs2_3 = R3p;
Rs2_2 = R2;
Rs2_1 = Rs2_2 + Rs2_3;

Rs32_2 = R1m;
Rs32_3 = R1p;
Rs32_1 = Rs32_2 + Rs32_3;

Rs31_3 = Rs32_1;
Rs31_2 = 1/(2*Fs*C1);
Rs31_1 = Rs31_2 + Rs31_3;

% Q ports
RQ4 = 1/(2*Fs*C2);
RQ5 = R4;
RQ6 = 1/(2*Fs*C3);

%% Calculate Scattering matrix
%syms Ga Gb Gc Gd Ge Gf;
syms Ga
Gb = 1/Rs2_1;
Gc = 1/Rs31_1;
Gd = 1/RQ4;
Ge = 1/RQ5;
Gf = 1/RQ6;

Nnodes = 12;
Y = sym(zeros(9, 9));

Y(1,1) = Ga;
Y(1,2) = -Ga;
Y(2,1) = -Ga;
Y(2,2) = Ga + Gb;
Y(2,4) = -Gb;
Y(3,3) = Gc + Gd;
Y(3,5) = -Gc;
Y(3,7) = -Gd;
Y(4,2) = -Gb;
Y(4,4) = Gb;
Y(5,3) = -Gc;
Y(5,5) = Gc;
Y(6,6) = Gf;
Y(6,9) = -Gf;
Y(7,3) = -Gd;
Y(7,7) = Gd;
Y(8,8) = Ge;
Y(9,6) = -Gf;
Y(9,9) = Gf;

A = sym(zeros(9,6));
A(1,1) = 1;
A(2,6) = -1;
A(3,2) = -1;
A(4,2) = 1;
A(5,3) = 1;
A(6,4) = -1;
A(6,5) = -1;
A(7,4) = 1;
A(8,5) = 1;
A(9,6) = 1;

B = transpose(A);

D = zeros(6,6);

X = [Y, A; B, D];

%% Port resitances
Rp1 = 1/Ga;
Rp2 = Rs2_1;
Rp3 = Rs31_1;
Rp4 = RQ4;
Rp5 = RQ5;
Rp6 = RQ6;
Rp = diag([Rp1, Rp2, Rp3, Rp4, Rp5, Rp6]);

S = eye(6) + 2*[zeros(6,9) Rp] * inv(X) * [zeros(6,9) eye(6,6)]';

eqn = S(1,1) == 0;
Ga_sol = solve(eqn, Ga)
R_adp = double(vpa(1/Ga_sol));

S = subs(S, Ga, Ga_sol);
S = double(vpa(S));

Rs1_3 = R_adp;
Rs1_2 = R3m;
Rs1_1 = Rs1_2 + Rs1_3;

%% initialize wave variables
% Series adaptor 1
as1_1 = 0; bs1_1 = 0;
as1_2 = 0; bs1_2 = 0;
as1_3 = 0; bs1_3 = 0;

% Series adaptor 2
as2_1 = 0; bs2_1 = 0;
as2_2 = 0; bs2_2 = 0;
as2_3 = 0; bs2_3 = 0;

% Series adaptor 31
as31_1 = 0; bs31_1 = 0;
as31_2 = 0; bs31_2 = 0;
as31_3 = 0; bs31_3 = 0;

% Series adaptor 32
as32_1 = 0; bs32_1 = 0;
as32_2 = 0; bs32_2 = 0;
as32_3 = 0; bs32_3 = 0;

% Component 4 (Q4)
aq4 = 0; bq4 = 0;
aq5 = 0; bq5 = 0;
aq6 = 0; bq6 = 0;

% R-type adaptor
br1 = 0; br2 = 0; br3 = 0;
br4 = 0; br5 = 0; br6 = 0;

ar1 = 0; ar2 = 0; ar3 = 0;
ar4 = 0; ar5 = 0; ar6 = 0;

% delay line inside capacitors
delay_c1 = 0;
delay_c2 = 0;
delay_c3 = 0;

% tmp
output_R1p = zeros(1,length(input));
output_R3m = zeros(1,length(input));
output_R3p = zeros(1,length(input));
output_R2 = zeros(1,length(input));

%% Simulation loop
for i = 1:N % run each time sample until N    
    % Wave up
    as2_2 = 0; as2_3 = 0; 
    bs2_1 = -as2_2 -as2_3;
    
    as32_2 = 0; as32_3 = 0; 
    bs31_1 = -as32_2 - as32_3;
    
    as31_2 = delay_c1; as31_3 = bs31_1;
    
    ar2 = bs2_1;
    ar3 = -as31_2 - as31_3;
    ar4 = delay_c2;
    ar5 = 0;
    ar6 = delay_c3;
    ar1 = bs1_3;
    
    tmp = S * [ar1; ar2; ar3; ar4; ar5; ar6];
    br1 = tmp(1);
    br2 = tmp(2);
    br3 = tmp(3);
    br4 = tmp(4);
    br5 = tmp(5);
    br6 = tmp(6);
    
    % Root - This might be the problem...
    as1_2 = 0; as1_3 = br1;
    bs1_1 = -as1_2 - as1_3;
    as1_1 = 2*input(i) - bs1_1;
    
    % Wave down
    [bs1_1, bs1_2, bs1_3] = series3(as1_1, Rs1_1, as1_2, Rs1_2, as1_3, Rs1_3);
    ar1 = bs1_3;
    
    tmp = S * [ar1; ar2; ar3; ar4; ar5; ar6];
    br1 = tmp(1);
    br2 = tmp(2);
    br3 = tmp(3);
    br4 = tmp(4);
    br5 = tmp(5);
    br6 = tmp(6);
    
    as2_1 = br2;
    [bs2_1, bs2_2, bs2_3] = series3(as2_1, Rs2_1, as2_2, Rs2_2, as2_3, Rs2_3);
    
    as31_1 = br3;
    [bs31_1, bs31_2, bs31_3] = series3(as31_1, Rs31_1, as31_2, Rs31_2, as31_3, Rs31_3);

    as32_1 = bs31_3;
    [bs32_1, bs32_2, bs32_3] = series3(as32_1, Rs32_1, as32_2, Rs32_2, as32_3, Rs32_3);
    
    % Store capacitor values
    delay_c1 = bs31_2;
    delay_c2 = br4;
    delay_c3 = br6;
    
    % Store output values
    output_R1p(i) = 0.5*(as32_3 + bs32_3);
    output_R3m(i) = 0.5*(as1_2 + bs1_2);
    output_R3p(i) = 0.5*(as2_3 + bs2_3);
    output_R2(i) = 0.5*(as2_2 + bs2_2);
    output(i) = output_R3m(i) - output_R3p(i) - output_R2(i) + output_R1p(i);
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
filename = 'output/bassman_tone_stack_double_precision_96k_32bit.wav'
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