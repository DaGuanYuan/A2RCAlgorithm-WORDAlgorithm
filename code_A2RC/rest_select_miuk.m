N = 16;     %the number of the elements
c = physconst('LightSpeed');    %the speed of light
Lambda = 10;    %wavelength
freq = c / Lambda;  %frequency
Rho_k_dB = -30;     %target in dB
Rho_k = 10^(Rho_k_dB/10);   %target
Theta_0 = 0;
Theta_k = pi * 20 /180;     %in radium
ElementSpace = Lambda / 2;
Phi_k = 2 * pi * ElementSpace * sin(Theta_k) / Lambda;
Phi_0 = 2 * pi * ElementSpace * sin(Theta_0) / Lambda;

a_k = zeros(16 ,1);
a_0 = zeros(16, 1);

for j = 1 : 1 : N
    a_k(j) = exp(Phi_k * 1i * (N - j));     %the standard point refers to the last point
    a_0(j) = exp(Phi_0 * 1i * (N - j));
end
w_0 = a_0;

P_k = (norm(w_0' * a_k, 2))^2;     %P_star_k-1 directed at theta_k
P_0 = (norm(w_0' * a_0, 2))^2;     %P_star_k-1 directed at theta_0
v = a_0' * a_k;     %inner product of a_0 and a_k
d_k = w_0'*a_k*(norm(a_k, 2))^2 - Rho_k*w_0'*a_0*v;     %d_k
Q_k = [P_k - Rho_k * P_0, d_k; conj(d_k), (norm(a_k, 2))^4 - Rho_k*(norm(v, 2))^2];    %Qk
