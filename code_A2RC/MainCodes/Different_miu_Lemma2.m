N = 16; %number of the elements 
Lambda = 10; %wavelength
c = 3e8; %speed of light 
freq = c / Lambda; %frequency
angle = -90 : 90;
Theta_0 = 0;
Theta_i = 20; 
angle_0 = [Theta_0; 0]; %theta_0
angle_i = [Theta_i; 0]; %theta_i
ElementSpacing = Lambda / 2; %isotropic elements spaces
array = phased.ULA(...      
 'NumElements', N, ...
 'ElementSpacing', ElementSpacing);
steervec = phased.SteeringVector('SensorArray', array);
a_0 = step(steervec, freq, angle_0);   %a_0
a_i = step(steervec, freq, angle_i);   %a_i
mu1 = 0.0119 + 0.0573i;
mu2 = -0.2974 + 0.0794i;
w_star1 = a_0 + mu1 * a_i;  %for mu1 = 0.0119 + j0.0573
w_star2 = a_0 + mu2 * a_i;  %for mu2 = -0.2974 + j0.0794
AR1 = phased.ArrayResponse( ...
      'SensorArray', array, ...
      'WeightsInputPort', true);     %array response 1
AR2 = phased.ArrayResponse( ...
      'SensorArray', array, ...
      'WeightsInputPort', true);     %array response 2
AR_init = phased.ArrayResponse( ...
          'SensorArray', array, ...
          'WeightsInputPort', true);    %initial array response
      
%%我取了个1-范数=-=
resp1 = abs(step(AR1, freq, angle, w_star1));     %response 1
resp2 = abs(step(AR2, freq, angle, w_star2));     %response 2
resp_init = abs(step(AR_init, freq, angle, a_0));    %initial pattern

resp1 = 10*log10(resp1);    %dB
resp2 = 10*log10(resp2);    %dB
resp_init = 10*log10(resp_init);    %dB

figure(1);  %Fig. 2. Adjust L_★(20°,0°) to -30dB
plot(angle, resp1, 'r');
xlabel('Angle of Arrival (degrees)')
ylabel('Array Gain (dB)');
title('Fig. 2. Adjust L_★(20°,0°) to -30dB');
hold on;
plot(angle, resp_init, 'b--');
legend('modified pattern', 'initial pattern');
hold off;
figure(2);
plot(angle, resp2, 'r');
xlabel('Angle of Arrival (degrees)');
ylabel('Array Gain (dB)');
title('Fig. 3. Adjust L_★(20°,0°) to -10dB');
hold on;
plot(angle, resp_init, 'b--');
legend('modified pattern', 'initial pattern');
hold off;