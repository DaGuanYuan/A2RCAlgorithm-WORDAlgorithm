%% initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 18;     %number of the elments 
Lambda = 10;    %wavelength
c = physconst('LightSpeed');    %the speed of light
freq = c / Lambda;   %frequency
ElementSpacing = Lambda / 2;    %isotropic elements spaces
Theta_0 = pi * 45 / 180;     %This angle could determine which position the pattern steers at
phi_0 = 2 * pi * ElementSpacing * sin(Theta_0) / Lambda;
a_0 = zeros(N, 1);
for j = 1 : 1 : N
   a_0(j) = exp(phi_0 * 1i * (N - j));    %reference point would be set at the last sensor 
end
w_0 = a_0;
Theta = zeros(1, 9001);     %every singel angle of arrival
phi = zeros(1, 9001);
a = zeros(N, 9001);    %steering vectors
L_0 = zeros(1, 9001);   %initial response

for angle_arr = -90 : 0.02 :90
   Theta(round((angle_arr+90)/0.02 + 1)) = pi * angle_arr / 180;  %every single angle of arrival
   phi(round((angle_arr+90)/0.02 + 1)) = 2 * pi * ElementSpacing * sin(Theta(round((angle_arr+90)/0.02 + 1)))/Lambda;
   
   for j = 1 : 1 : N
      a(j, round((angle_arr+90)/0.02 + 1)) = exp(phi(round((angle_arr+90)/0.02 + 1)) * 1i * (N - j));    %reference point would be set at the last sensor
   end
   
   L_0(round((angle_arr+90)/0.02 + 1)) = (norm((w_0' * a(:, round((angle_arr+90)/0.02 + 1))), 2))^2/(norm(w_0' * a_0, 2))^2;
end
% L_0_dB = 10*log10(L_0);   %initial response in dB
% angle_arr = -90 : 0.02 :90;
% figure(1);
% plot(angle_arr, L_0_dB, 'r-');
% set(gca, 'YLim', [-40, 0]);
% hold on;

%% Add First Main-lobe steering at 45¡ã
% % find theta_k
% % omiga_m = -7 : 0.02 : 7;     %omiga_m
% % omiga_s = angle_arr;
% % omiga_s(4151 : 1 : 4850) = [];      %omiga_s
% Rho_k_0 = 0;      %-25dB
% Rho_k = 10^(Rho_k_0/10);
% % srl_nmb = round((omiga_s + 90)/0.02 + 1);
% % [max_L_0, index_0] = max(L_0_dB(srl_nmb));
% % if index_0 >4150
% %    index_0 = index_0 + 700; 
% % end
% Theta_0=45;
% theta_0 = pi * Theta_0 / 180;
% phi_theta0 = 2 * pi *ElementSpacing * sin(theta_0) / Lambda;
% a_theta0 = zeros(N, 1);
% for j = 1 : 1 : N
%     a_theta0(j) = exp(phi_theta0 * 1i * (N-j));
% end
% 
% % find beta_k
% P_prl = a_theta0*a_theta0'/(norm(a_theta0, 2)^2);     %P_star_k-1 directed at theta_k
% P_vtc = eye(N) - P_prl;     %P_star_k-1 directed at theta_0
% w_prl = P_prl * w_0;
% w_vtc = P_vtc * w_0;
% B = [-Rho_k*norm(w_vtc' * a_0)^2, -Rho_k*w_vtc'*(a_0*a_0')*w_prl; -Rho_k*w_prl'*(a_0*a_0')*w_vtc, norm(w_prl'*a_theta0, 2)^2 - Rho_k*norm(w_prl'*a_0, 2)^2];
% d = sqrt((real(B(1,2)))^2 - B(1,1)*B(2,2));
% beta_a = (-real(B(1,2))+d)/B(2,2);
% beta_b = (-real(B(1,2))-d)/B(2,2);
% 
% % w_k and L_k
% w_a = [w_vtc, w_prl]*[1, beta_a].';
% w_b = [w_vtc, w_prl]*[1, beta_b].';
% P_0 = w_0*w_0'/(norm(w_0, 2)^2);
% F_a = norm(P_0 * w_a/norm(w_a, 2), 2)^2;
% F_b = norm(P_0 * w_b/norm(w_b, 2), 2)^2;
% 
% w_k = w_a;
% if (F_a >F_b)
%    w_k = w_b; 
% end
% 
% temp2 = zeros(1, 9001);
% L_theta0 = zeros(1, 9001);
% for angle_arr = -90 : 0.02 :90
%    Theta(round((angle_arr+90)/0.02 + 1)) = pi * angle_arr / 180;  %every single angle of arrival
%    phi(round((angle_arr+90)/0.02 + 1)) = 2 * pi * ElementSpacing * sin(Theta(round((angle_arr+90)/0.02 + 1)))/Lambda;
%    
%    for j = 1 : 1 : N
%       a(j, round((angle_arr+90)/0.02 + 1)) = exp(phi(round((angle_arr+90)/0.02 + 1)) * 1i * (N - j));    %reference point would be set at the last sensor
%    end
%    temp2(round((angle_arr+90)/0.02 + 1)) = (norm((w_k' * a(:, round((angle_arr+90)/0.02 + 1))), 2))^2;
% end
% 
% for angle_arr = -90 : 0.02 : 90
%    L_theta0(round((angle_arr+90)/0.02 + 1)) = (norm((w_k' * a(:, round((angle_arr+90)/0.02 + 1))), 2))^2/max(temp2);    
% end
% 
% [~, index_0] = max(temp2);
% theta_0 = (index_0-1) * 0.02 - 90;
% Theta_0 = pi * theta_0 / 180;
% phi_0 = 2 * pi *ElementSpacing * sin(Theta_0) / Lambda;
% a_0 = zeros(N, 1);
% 
% for j = 1 : 1 : N
%     a_0(j) = exp(phi_0 * 1i * (N-j));
% end
% 
% L_theta0_dB = 10*log10(L_theta0);   %initial response in dB
% angle_arr = -90 : 0.02 :90;
% figure(1);
% plot(angle_arr, L_theta0_dB, 'b');
% set(gca, 'YLim', [-40, 0]);
% title('Fig.3. Single Point Control of Sidelobe using WORD Algorithm');
% xlabel('Angle of Arrival(¡ã)');
% ylabel('L(dB)');
% legend('Initial','The 1st step of WORD Algorithm');
% hold on;

%% Add First Main-lobe steering at -45¡ã
% find theta_k
% omiga_m = -7 : 0.02 : 7;     %omiga_m
% omiga_s = angle_arr;
% omiga_s(4151 : 1 : 4850) = [];      %omiga_s
Rho_k_0 = 0;      %-25dB
Rho_k = 10^(Rho_k_0/10);
% srl_nmb = round((omiga_s + 90)/0.02 + 1);
% [max_L_0, index_0] = max(L_0_dB(srl_nmb));
% if index_0 >4150
%    index_0 = index_0 + 700; 
% end
Theta_0 = -45;
theta_0 = pi * Theta_0 / 180;
phi_theta0 = 2 * pi *ElementSpacing * sin(theta_0) / Lambda;
a_theta0 = zeros(N, 1);
for j = 1 : 1 : N
    a_theta0(j) = exp(phi_theta0 * 1i * (N-j));
end

% find beta_k
P_prl = a_theta0*a_theta0'/(norm(a_theta0, 2)^2);     %P_star_k-1 directed at theta_k
P_vtc = eye(N) - P_prl;     %P_star_k-1 directed at theta_0
w_prl = P_prl * w_0;
w_vtc = P_vtc * w_0;
B = [-Rho_k*norm(w_vtc' * a_0)^2, -Rho_k*w_vtc'*(a_0*a_0')*w_prl; -Rho_k*w_prl'*(a_0*a_0')*w_vtc, norm(w_prl'*a_theta0, 2)^2 - Rho_k*norm(w_prl'*a_0, 2)^2];
d = sqrt((real(B(1,2)))^2 - B(1,1)*B(2,2));
beta_a = (-real(B(1,2))+d)/B(2,2);
beta_b = (-real(B(1,2))-d)/B(2,2);

% w_k and L_k
w_a = [w_vtc, w_prl]*[1, beta_a].';
w_b = [w_vtc, w_prl]*[1, beta_b].';
P_0 = w_0*w_0'/(norm(w_0, 2)^2);
F_a = norm(P_0 * w_a/norm(w_a, 2), 2)^2;
F_b = norm(P_0 * w_b/norm(w_b, 2), 2)^2;

w_k = w_a;
if (F_a >F_b)
   w_k = w_b; 
end

temp2 = zeros(1, 9001);
L_theta0 = zeros(1, 9001);
for angle_arr = -90 : 0.02 :90
   Theta(round((angle_arr+90)/0.02 + 1)) = pi * angle_arr / 180;  %every single angle of arrival
   phi(round((angle_arr+90)/0.02 + 1)) = 2 * pi * ElementSpacing * sin(Theta(round((angle_arr+90)/0.02 + 1)))/Lambda;
   
   for j = 1 : 1 : N
      a(j, round((angle_arr+90)/0.02 + 1)) = exp(phi(round((angle_arr+90)/0.02 + 1)) * 1i * (N - j));    %reference point would be set at the last sensor
   end
   temp2(round((angle_arr+90)/0.02 + 1)) = (norm((w_k' * a(:, round((angle_arr+90)/0.02 + 1))), 2))^2;
end

for angle_arr = -90 : 0.02 : 90
   L_theta0(round((angle_arr+90)/0.02 + 1)) = (norm((w_k' * a(:, round((angle_arr+90)/0.02 + 1))), 2))^2/max(temp2);    
end

L_theta0_dB = 10*log10(L_theta0);   %initial response in dB
angle_arr = -90 : 0.02 :90;
figure(1);
plot(angle_arr, L_theta0_dB, 'r');
% set(gca, 'YLim', [-40, 0]);
% title('Fig.6. Pattern with 2 mainlobes');
% xlabel('Angle of Arrival(¡ã)');
% ylabel('L(dB)');
% legend('Initial','Pattern with 2 mainlobes');
hold on;

%% Synthesize of Mainlobe by iterating
L_0_dB = L_theta0_dB;   %set pattern with two mainlobes as the initial pattern
L_temp = abs(L_0_dB - (-3));    %to change the value
index_M = zeros(1, 4);
for j = 1 : 1 : 4
    for angle_arr = -90 : 0.02 :90
       [temp, index_M(j)] = min(L_temp);
    end
    
    L_temp(index_M(j)) = 100;    %eliminate the first 3 dB value
end
index_M = sort(index_M);


Rho_k_0 = 0;      %0dB
dn = 1;      %-1dB deep null
Rho_k = 10^(Rho_k_0/10);
index_m = [index_M(1) : 1 : index_M(2), index_M(3) : 1 :  index_M(4)];


% srl_nmb = round((omiga_s + 90)/0.02 + 1);   %index of omiga_s

L_dB = L_0_dB;  %for the sake of clarity ¡ª¡ª>iteration
[min_M, index_M2] = min(L_dB(index_m));

while(abs(min_M - 0) >= dn)
%  for i = 1 : 1 :300
    if index_M2 <= (index_M(2) - index_M(1) + 1)
        index_M2 = index_M2 + index_M(1) - 1; 
    else
        index_M2 = index_M2 + index_M(3) - (index_M(2) - index_M(1) + 1) - 1;
    end
    
    theta_k = (index_M2-1) * 0.02 - 90;
    Theta_k = pi * theta_k / 180;
    phi_k = 2 * pi *ElementSpacing * sin(Theta_k) / Lambda;
    a_k = zeros(N, 1);
    
    for j = 1 : 1 : N
        a_k(j) = exp(phi_k * 1i * (N-j));
    end

    %find beta_k
    P_prl = a_k*a_k'/(norm(a_k, 2)^2);     %P_star_k-1 directed at theta_k
    P_vtc = eye(N) - P_prl;     %P_star_k-1 directed at theta_0
    w_prl = P_prl * w_0;
    w_vtc = P_vtc * w_0;
    B = [-Rho_k*norm(w_vtc' * a_0)^2, -Rho_k*w_vtc'*(a_0*a_0')*w_prl; -Rho_k*w_prl'*(a_0*a_0')*w_vtc, norm(w_prl'*a_k, 2)^2 - Rho_k*norm(w_prl'*a_0, 2)^2];
    d = sqrt((real(B(1,2)))^2 - B(1,1)*B(2,2));
    beta_a = (-real(B(1,2))+d)/B(2,2);
    beta_b = (-real(B(1,2))-d)/B(2,2);

    %w_k and L_k
    w_a = [w_vtc, w_prl]*[1, beta_a].';
    w_b = [w_vtc, w_prl]*[1, beta_b].';
    P_0 = w_0*w_0'/(norm(w_0, 2)^2);
    F_a = norm(P_0 * w_a/norm(w_a, 2), 2)^2;
    F_b = norm(P_0 * w_b/norm(w_b, 2), 2)^2;

    w_k = w_a;
    if (F_a >F_b)
       w_k = w_b; 
    end

    temp2 = zeros(1, 9001);
    L_k = zeros(1, 9001);
    for angle_arr = -90 : 0.02 :90
       Theta(round((angle_arr+90)/0.02 + 1)) = pi * angle_arr / 180;  %every single angle of arrival
       phi(round((angle_arr+90)/0.02 + 1)) = 2 * pi * ElementSpacing * sin(Theta(round((angle_arr+90)/0.02 + 1)))/Lambda;

       for j = 1 : 1 : N
          a(j, round((angle_arr+90)/0.02 + 1)) = exp(phi(round((angle_arr+90)/0.02 + 1)) * 1i * (N - j));    %reference point would be set at the last sensor
       end
       temp2(round((angle_arr+90)/0.02 + 1)) = (norm((w_k' * a(:, round((angle_arr+90)/0.02 + 1))), 2))^2;
    end

    for angle_arr = -90 : 0.02 : 90
       L_k(round((angle_arr+90)/0.02 + 1)) = (norm((w_k' * a(:, round((angle_arr+90)/0.02 + 1))), 2))^2/max(temp2);    
    end

    L_k_dB = 10*log10(L_k);   %initial response in dB
    
    %iteration
    L_dB = L_k_dB;
    w_0 = w_k;
    [min_M, index_M2] = min(L_dB(index_m));

end

angle_arr = -90 : 0.02 :90;
figure(1);
plot(angle_arr, L_dB, 'b','LineWidth', 1);
set(gca, 'YLim', [-40, 0]);
title('Fig.7.Mainlobe Control');
xlabel('Angle of Arrival(¡ã)');
ylabel('L(dB)');
legend('Initial', 'Result of WORD Algorithm');
hold on;
