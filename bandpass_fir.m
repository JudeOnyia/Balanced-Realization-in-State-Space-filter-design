% This function implements the least-squares design of nonlinear-phase
% lowpass FIR filters.
% Inputs:
% N: Order of the FIR filter, N must be an even integer.
% d: group-dealy in passband
% w: weight for the stopbands
% (assuming the weight in passband is 1)
% fa1: normalized lower stopband edge between 0 and pi
% fp1: normalized lower passband edge between 0 and pi
% with omi_p1 > omi_a1
% fp2: normalized higher passband edge between 0 and pi
% with omi_p2 > omi_p1
% fa2: normalized higher cutoff frequency between 0 and pi
% with omi_c2 > omi_p2.
% Output:
% h: impulse response of the bandpass FIR filter.
% Written by W.-S. Lu, University of Victoria
function h = bandpass_fir(N,d,w,fa1,fp1,fp2,fa2)
a1 = pi*fa1;
p1 = pi*fp1;
p2 = pi*fp2;
a2 = pi*fa2;
q = zeros(N+1,N+1);
b = zeros(N+1,1);
c = zeros(N+1,1);
c(1) = w*(a1-a2+pi)+p2-p1;
b(1) = (sin(d*p2)-sin(d*p1))/d;
for i = 1:N,
 z1 = sin(i*a1) - sin(i*a2);
 z2 = sin(i*p2 )- sin(i*p1);
 c(i+1) = (w*z1 + z2)/i;
 if i == d,
 b(i+1) = p2 - p1;
 else
 b(i+1) = (sin((i-d)*p2) -sin((i-d)*p1))/(i-d);
 end
end
Q = toeplitz(c);
h = inv(Q)*b;