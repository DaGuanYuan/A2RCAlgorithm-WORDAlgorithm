N = 16; %number of the elements 
Lambda = 10; %wavelength
c = physconst('LightSpeed'); %speed of light 
freq = c / Lambda; %frequency
Rho_k_0 = -30;    %¦Ñ_k = -30dB
Rho_k = 10^(Rho_k_0/10);
Theta_0 = 0;
Theta_k = 20; 
angle_0 = [Theta_0; 0]; %theta_0
angle_k = [Theta_k; 0]; %theta_k
ElementSpacing = Lambda / 2; %isotropic elements spaces
array = phased.ULA(...      
 'NumElements', N, ...
 'ElementSpacing', ElementSpacing);
steervec = phased.SteeringVector('SensorArray', array);
a_0 = step(steervec, freq, angle_0);  %a_0
a_k = step(steervec, freq, angle_k);  %a_k
w_0 = a_0;

P_k = (norm(w_0' * a_k, 2))^2;     %P_star_k-1 directed at theta_k
P_0 = (norm(w_0' * a_0, 2))^2;     %P_star_k-1 directed at theta_0
v = a_0' * a_k;     %inner product of a_0 and a_k
d_k = w_0'*a_k*(norm(a_k, 2))^2 - Rho_k*w_0'*a_0*v;     %d_k
Q_k = [P_k - Rho_k * P_0, d_k; conj(d_k), (norm(a_k, 2))^4 - Rho_k*(norm(v, 2))^2];    %Qk

w_k = zeros(length(w_0), 9001); 
a = zeros(length(w_0), 9001);
L_k = zeros(9001, 9001);
L_0 = zeros(1, 9001);
Theta_Phi = zeros(1, 9001);
J_k = zeros(9001, 9001);
miu_k = zeros(1, 9001);     %initialization
for i = 1 : 1 : 9001
    Phi = 2*pi*(i-1)/9000;
    miu_k(i) = -conj(Q_k(1, 2))/Q_k(2, 2) + exp(Phi*1i)*sqrt(-det(Q_k))/norm(Q_k(2, 2), 2);    %miu_k    
    w_k(:, i) = w_0 + miu_k(i) * a_k;
end

for omiga = -90 : 0.02 : 90
    Theta_Phi(round((omiga+90)/0.02+1)) = omiga;
    a(:, round((omiga+90)/0.02+1)) = step(steervec, freq, [Theta_Phi(round((omiga+90)/0.02+1)); 0]);
%     L_k(round((Phi+90)/0.02+1)) = norm(w_k(:, round((Phi+90)/0.02+1))'* a(:, round((Phi+90)/0.02+1)) ,2)/norm(w_k(:, round((Phi+90)/0.02+1))'*a_0, 2);
    L_0(round((omiga+90)/0.02+1)) = (norm(w_0'* a(:, round((omiga+90)/0.02+1)) ,2))^2/(norm(w_0'*a_0, 2))^2;    %L_¡ï(k-1)(¦Èi£¬¦È0)
%     J(round((Phi+90)/0.02+1)) = sum(abs(L_k(round((Phi+90)/0.02+1)) -  L_0(round((Phi+90)/0.02+1))))/9001;
end

for i = 1 : 1 : 9001
   for  j = 1 : 1 : 9001
      L_k(i, j) = (norm(w_k(:, i)' * a(:, j), 2))^2/(norm(w_k(:, i)'*a_0, 2))^2;    %L_¡ï(k)(¦Èi£¬¦È0)
   end
   J_k(i, :) = abs(L_k(i, :) - L_0);      %every single J regarding the average cost function
end

J = sum(J_k, 2)/9001;
Phi = 2*pi*(0:1:9000)/9000;
figure(1);
semilogy(Phi, J, 'b--');
title('Fig.5. J versus ¦µ£¨L_¡ï(20¡ã, 0¡ã) = -30dB£©.');
xlabel('¦µ(rad)');
ylabel('J');