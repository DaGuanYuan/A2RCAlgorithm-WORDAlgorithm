%% initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 18;     %number of the elments 
Lambda = 10;    %wavelength
c = physconst('LightSpeed');    %the speed of light
freq = c / Lambda;   %frequency
ElementSpacing = Lambda / 2;    %isotropic elements spaces
Theta_0 = pi * 20 / 180;     %This angle could determine which position the pattern steers at
phi_0 = 2 * pi * ElementSpacing * sin(Theta_0) / Lambda;
a_0 = zeros(N, 1);
for j = 1 : 1 : N
   a_0(j) = exp(phi_0 * 1i * (N - j));    %reference point would be set at the last sensor 
end

a_0(3) = [];    %delete the third sensor
a_0(7) = [];    %delete the third sensor
a_0(9) = [];    %delete the third sensor
a_0(13) = [];    %delete the third sensor
a_0(14) = [];    %delete the third sensor
% a_0(15) = [];    %delete the third sensor
% a_0(19) = [];    %delete the third sensor

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
   
end

a(3,:) = [];
a(7,:) = [];
a(9,:) = [];
a(13,:) = [];
a(14,:) = [];
% a(15,:) = [];
% a(19,:) = [];

for angle_arr = -90 : 0.02 :90

   L_0(round((angle_arr+90)/0.02 + 1)) = (norm((w_0' * a(:, round((angle_arr+90)/0.02 + 1))), 2))^2/(norm(w_0' * a_0, 2))^2;
end
% L_0_dB = 10*log10(L_0);   %initial response in dB
% angle_arr = -90 : 0.02 :90;
% figure(1);
% plot(angle_arr, L_0_dB, 'r-');
% set(gca, 'YLim', [-40, 0]);
% hold on;

% %% Single Point Control of Sidelobe
% %find theta_k
% omiga_m = -7 : 0.02 : 7;     %omiga_m
% omiga_s = angle_arr;
% omiga_s(4151 : 1 : 4850) = [];      %omiga_s
% Rho_k_0 = -25;      %-25dB
% Rho_k = 10^(Rho_k_0/10);
% srl_nmb = round((omiga_s + 90)/0.02 + 1);
% [max_L_0, index_0] = max(L_0_dB(srl_nmb));
% if index_0 >4150
%    index_0 = index_0 + 700; 
% end
% theta_0 = pi * ((index_0-1) * 0.02 - 90) / 180;
% phi_theta0 = 2 * pi *ElementSpacing * sin(theta_0) / Lambda;
% a_theta0 = zeros(N, 1);
% for j = 1 : 1 : N
%     a_theta0(j) = exp(phi_theta0 * 1i * (N-j));
% end
% 
% %find miu_k
% P_theta0 = (norm(w_0' * a_theta0, 2))^2;     %P_star_k-1 directed at theta_k
% P_0 = (norm(w_0' * a_0, 2))^2;     %P_star_k-1 directed at theta_0
% v_theta0 = a_0' * a_theta0;     %inner product of a_0 and a_k
% d_theta0 = w_0'*a_theta0*(norm(a_theta0, 2))^2 - Rho_k*w_0'*a_0*v_theta0;     %d_k
% Q_theta0 = [P_theta0 - Rho_k * P_0, d_theta0; conj(d_theta0), (norm(a_theta0, 2))^4 - Rho_k*(norm(v_theta0, 2))^2];    %Q_theta_0
% 
% c_miu0 = (1/Q_theta0(2,2))*[-real(Q_theta0(1, 2)), imag(Q_theta0(1, 2))];
% r_miu0 = sqrt(-det(Q_theta0))/norm(Q_theta0(2,2), 2);
% temp = c_miu0 * (norm(c_miu0, 2)-r_miu0)/norm(c_miu0, 2);
% miu0 = temp(1) + temp(2)*1i;
% 
% %w_k and L_k
% w_new = w_0 + miu0 * a_theta0;
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
%    temp2(round((angle_arr+90)/0.02 + 1)) = (norm((w_new' * a(:, round((angle_arr+90)/0.02 + 1))), 2))^2;
% end
% 
% for angle_arr = -90 : 0.02 : 90
%    L_theta0(round((angle_arr+90)/0.02 + 1)) = (norm((w_new' * a(:, round((angle_arr+90)/0.02 + 1))), 2))^2/max(temp2);    
% end
% 
% L_theta0_dB = 10*log10(L_theta0);   %initial response in dB
% angle_arr = -90 : 0.02 :90;
% figure(1);
% plot(angle_arr, L_theta0_dB, 'b');
% set(gca, 'YLim', [-40, 0]);
% title('Fig.1. Single Point Control of Sidelobe using A^2RC Algorithm');
% xlabel('Angle of Arrival(¡ã)');
% ylabel('L(dB)');
% legend('Initial','The 1st step of A^2RC');
% hold on;

%% Synthesize of Sidebole by iterating
omiga_s = angle_arr;
omiga_s(4001 : 1 : 5000) = [];      %omiga_s
Rho_k_0 = -25.1;      %-25dB
Rho_k = 10^(Rho_k_0/10);
srl_nmb = round((omiga_s + 90)/0.02 + 1);   %index of omiga_s

L_dB = L_0_dB;  %for the sake of clarity ¡ª¡ª>iteration
[max_L, index_L] = max(L_dB(srl_nmb));
while(max_L >= -25)
    if index_L >4000
       index_L = index_L + 1000; 
    end
    theta_k = (index_L-1) * 0.02 - 90;
    Theta_k = pi * theta_k / 180;
    phi_k = 2 * pi *ElementSpacing * sin(Theta_k) / Lambda;
    a_k = zeros(N, 1);
    
    for j = 1 : 1 : N
        a_k(j) = exp(phi_k * 1i * (N-j));
    end
    
    %find miu_k
    P_k = (norm(w_0' * a_k, 2))^2;     %P_star_k-1 directed at theta_k
    P_0 = (norm(w_0' * a_0, 2))^2;     %P_star_k-1 directed at theta_0
    v_k = a_0' * a_k;     %inner product of a_0 and a_k
    d_k = w_0'*a_k*(norm(a_k, 2))^2 - Rho_k*w_0'*a_0*v_k;     %d_k
    Q_k = [P_k - Rho_k * P_0, d_k; conj(d_k), (norm(a_k, 2))^4 - Rho_k*(norm(v_k, 2))^2];    %Q_theta_0

    c_k = (1/Q_k(2,2))*[-real(Q_k(1, 2)), imag(Q_k(1, 2))];
    r_k = sqrt(-det(Q_k))/norm(Q_k(2,2), 2);
    temp = c_k * (norm(c_k, 2)-r_k)/norm(c_k, 2);
    miu_k = temp(1) + temp(2)*1i;

    %w_k and L_k
    w_k = w_0 + miu_k * a_k;

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
    [max_L, index_L] = max(L_dB(srl_nmb));
    %update a_0
    [~, index_0] = max(temp2);
    theta_0 = (index_0-1) * 0.02 - 90;
    Theta_0 = pi * theta_0 / 180;
    phi_0 = 2 * pi *ElementSpacing * sin(Theta_0) / Lambda;
    a_0 = zeros(N, 1);
    
    for j = 1 : 1 : N
        a_0(j) = exp(phi_0 * 1i * (N-j));
    end
end

angle_arr = -90 : 0.02 :90;
figure(1);
plot(angle_arr, L_dB, 'b','LineWidth', 1);
set(gca, 'YLim', [-40, 0]);
title('Fig.2. Result of A^2RC Algorithm');
xlabel('Angle of Arrival(¡ã)');
ylabel('L(dB)');
legend('Initial','Result of A^2RC Algorithm');
hold off;