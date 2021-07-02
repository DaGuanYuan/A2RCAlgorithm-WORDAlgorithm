N = 12; %number of the elements 
Lambda = 10; %wavelength
c = 3e8; %speed of light 
freq = c / Lambda; %frequency
Theta_0 = 0;
Theta_i = pi * 45/180; 
% angle_0 = [Theta_0; 20]; %theta_0
% angle_i = [Theta_i; 20]; %theta_i
ElementSpacing = Lambda / 2; %isotropic elements spaces
phi_0 = 2 * pi * ElementSpacing * sin(Theta_0) / Lambda;
phi_i = 2 * pi * ElementSpacing * sin(Theta_i) / Lambda;

a_0 = zeros(N, 1);
a_i = zeros(N, 1);

for j = 1 : 1 : N
   a_0(j) = exp(phi_0 * 1i * (N-j));      %the standard point refers to the last point
   a_i(j) = exp(phi_i * 1i * (N-j));
end

P_ai = (a_i * a_i')/((norm(a_i))^2);    %projection ai
P_ai_vtc = eye(N) - P_ai;    %vertical
w0 = a_0;   %initial weight vector
w0_prl = P_ai * w0;     %parallel
w0_vtc = P_ai_vtc * w0;     %vertical
beta = -700 : 1 : 700;
w_star = zeros(N, 1, length(beta));
L_star = zeros(1, length(beta));    %initialization
for beta = -700 : 0.1 : 700
    w_star(:,:,round(beta/0.1+7001)) = w0_vtc + beta * w0_prl;
    phi_i = acos( norm(w_star(:,:,round(beta/0.1+7001))'*a_i,1 )/(norm(  w_star(:,:,round(beta/0.1+7001)))*norm(a_i)));
    phi_0 = acos( norm(w_star(:,:,round(beta/0.1+7001))'*a_0,1 )/(norm( w_star(:,:,round(beta/0.1+7001)))*norm(a_0)));
    L_star(1, round(beta/0.1+7001)) = ((norm(a_i))^2 * (cos(phi_i))^2)/((norm(a_0))^2 * (cos(phi_0))^2);
end
L_star = 10*log10(L_star);  %dB

beta = -700 : 0.1 : 700;
plot(beta, L_star, 'r');
xlabel('汕');
ylabel('L(成_i,成_0)dB');
title('Fig.2. Curve of L_∴(成,成_0) versus 汕 for a ULA.');
% set(gca,'YLim',[-150 120]);
% Resp = step(array, freq, angle);


beta_p = -(norm(w0_vtc, 2)^2)/(norm(w0_prl, 2)^2);
L_q = norm((w0' * a_i), 2)^2/norm((w0' * a_0), 2)^2;
L_inf = norm((w0_prl' * a_i), 2)^2/norm((w0_prl' * a_0), 2)^2;
L_q_dB = 10*log10(L_q);
L_inf_dB = 10*log10(L_inf);
a = norm(a_i,2)^2/(L_q*norm(a_0,2)^2);