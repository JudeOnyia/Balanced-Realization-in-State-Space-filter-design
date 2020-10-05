close all
clear
clc

% Filter specifications
N = 28; % High-order of stable FIR filter
r = 12; % Low-order of obtained IIR filter
D = 7; % Group Delay
fc = 0.7; % cut-off frequency for High-pass filter
% upper and lower cut-off frequencies for band-pass and band-stop
fc1 = 0.2;
fc2 = 0.6;

% Display filter coefficients of High-pass IIR filter
[a,b] = High_Pass_bt_iir(N,r,D,fc);
format long
fprintf('Denominator of High Pass Filter (a):\n');
display(a);
fprintf('Numerator of High Pass Filter (b):\n');
display(b);

% Display filter coefficients of Band-pass IIR filter
[a,b] = Band_Pass_bt_iir(N,r,D,fc1,fc2);
format long
fprintf('Denominator of Band Pass Filter (a):\n');
display(a);
fprintf('Numerator of Band Pass Filter (b):\n');
display(b);

% Display filter coefficients of Band-stop IIR filter
[a,b] = Band_Stop_bt_iir(N,r,D,fc1,fc2);
format long
fprintf('Denominator of Band Stop Filter (a):\n');
display(a);
fprintf('Numerator of Band Stop Filter (b):\n');
display(b);