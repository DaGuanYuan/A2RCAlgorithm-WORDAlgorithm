N = 12;     %number of the elments 
Lambda = 10;    %wavelength
c = physconst('LightSpeed');    %the speed of light
freq = c / Lambda;   %frequency
ElementSpacing = Lambda / 2;    %isotropic elements spaces
Theta_0 = pi * 0 / 180;     %This angle could determine which position the pattern steers at
Theta_i = pi * 45 / 180;
phi_0 = 2 * pi * ElementSpacing * sin(Theta_0) / Lambda;
phi_i = 2 * pi * ElementSpacing * sin(Theta_i) / Lambda;
a_0 = zeros(N, 1);
a_i = zeros(N, 1); 
for j = 1 : 1 : N
   a_0(j) = exp(phi_0 * 1i * (N - j));    %reference point would be set at the last sensor 
   a_i(j) = exp(phi_i * 1i * (N - j));    %reference point would be set at the last sensor 
end
w_0 = a_0;

v = a_0'* a_i;  %inner product of a0 and ai
INR = -0.2 : 0.0001 : 0.2; 
L_0 = zeros(1, length(INR));
L_star = zeros(1, length(INR));     %initialization
for j = -0.2 : 0.0001 : 0.2
%     L_star(round(j/0.01 + 21)) = 
    L_star(round(j/0.001 + 2001)) = (norm(v,2))^2/(((norm(a_i,2)^2)*(norm(a_0,2)^2)-(norm(v,1)^2))*j + (norm(a_0,2)^2))^2  ; %L_star, i = INR
    L_0(round(j/0.001 + 2001)) = (norm(v,1))^2/(norm(a_0,2)^4);    %L0

end



semilogy(INR, L_0, 'b--');
hold on;
semilogy(INR, L_star, 'r');
title('Fig.1. Curves of Lбя and L0 versus INR for a ULA');
xlabel('INR');
ylabel('L');
legend('L_0','L_бя');
hold off;


