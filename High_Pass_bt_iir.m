% Design highpass IIR filter of order r with constant
% passband group delay using the balanced truncation technique.
% Inputs:
% N: order of the high-order FIR filter
% r: order of the IIR filter.
% D: desired group delay.
% fc: normalized cutoff frequency between 0 and 1.
% Outputs:
% a,b: denominator and numerator of the IIR filter.
% Written by W.-S. Lu, University of Victoria.
% Example: [a,b] = bt_iir(28,12,7,0.7);
function [a,b] = High_Pass_bt_iir(N,r,D,fc)
ep = 0.09;
h = bandpass_fir(N,D,1,fc-ep,fc+ep,0.9999,0.9999);

% plot FIR frequeny response
figure(1)
fp1 = fc+ep;
fp2 = 0.9999;
j = sqrt(-1);
f = 0:1/1023:1;
ff = f*pi;
ff = ff(:);
F = freqz(h,1,ff);
subplot(211)
amp = abs(F);
plot(f,20*log10(amp))
axis([0 1 -80 10])
title('Amplitude response of the Highpass fir filter')
xlabel('Normalized frequency')
ylabel('Amplitude response in dB')
grid
subplot(212)
i1 = ceil(1024*fp1) + 1;
i2 = floor(1024*fp2)+ 1;
Hdp = exp(-j*D*ff(i1:i2));
er_p = abs(F(i1:i2)-Hdp);
plot(f(i1:i2),er_p)
axis([fp1 fp2 0 2*max(er_p)])
grid
title('Passband Ripple of the Highpass fir filter')
xlabel('Normalized frequency')
ylabel('Amplitude response in passband') 


% State-space realization {A, b, c, d}
d = h(1);
b = h(2:end);
I = eye(N-1);
z1 = zeros(N-1,1);
A = [z1 I; 0 z1'];
c = [1 z1'];

% Finding Balanced realization {Ab, bb, cb, d}
Q = b*b';
K = dlyap(A,Q); % The controllability Gramian 
[U,S,V] = svd(K); % Orthogonal matrix U
s = diag(S);
s = s.^0.25;
V = diag(s); % Diagonal matrix V
T = U*V;
Ti = inv(T);
Ab = Ti*A*T;
bb = Ti*b;
cb = c*T;

% Finding reduced order of filter properties {Ar, br, cr, d} using
% the balanced truncation technique
Ar = Ab(1:r,1:r);
br = bb(1:r);
cr = cb(1:r);
[b,a] = ss2tf(Ar,br,cr,d,1);
a = a(:);
b = b(:);
[H,w] = freqz(b,a,1024);

% plot IIR frequeny response
figure(2)
plot(w/pi,20*log10(abs(H)))
grid
axis([0 1 -60 10])
xlabel('Normalized frequency')
ylabel('dB')
title("Amplitude response of the Highpass IIR filter")
figure(3)
[gd,w] = grpdelay(b,a,1024);
np = floor(1024*fc);
plot(w(np:1024)/pi,gd(np:1024));
axis([fc 1 2 20])
grid
xlabel('Normalized frequency')
ylabel('samples')
title("Passband group delay of the Highpass IIR filter")