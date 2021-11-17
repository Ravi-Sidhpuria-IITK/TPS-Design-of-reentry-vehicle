clc;
clear;
close all;


%% Calculation For constant properties

for m =671%:1001; %optimum no. of grid points is obtained here

T_i = 300; % Temperature of TPS at 0 seconds
T_s = 300; % Temperature of face plate at 0 second
T_back = 305; % Max temp. at bac face of TPS
T_inf = 0; % Space temperature
sigma = 5.67 * 10^-8; % Stefan Boltzmann constant
epsl = 0.9; % Emmisivity can be changed here


for x = 0.185  % x is varied here to obtain optimum thickness of TPS

T(1:m) = T_i; % Initial Temperatures
T_old = T; 

% Reinforced Carbon Carbon (RCC) Composite properties

rho = 1580;     % Constant Density
k(1:m) = 7.72;  %constant k & c
c(1:m) = 1485;
alpha = k(1)/(rho*c(1));
T_max = 2030; % Melting/max allowable Temperature of RCC

dx = x/(m-1); % Grid size
for j = 1:m
tvar(j) = ((dx)^2)/((2*(k(j)/(rho*c(j))))*(1 + (sigma*epsl*dx*T(j)^3)/k(j))); % calculation of variable time step due to variable properties
end
count = 1;
dt = min(tvar); % Smallest time step is taken for each iteration
t = dt; % calculation of cumulative time

while t < 600 && max(T) < T_max && min(T) < T_back % requirements of given TPS can be changed here
   if t < 300
       qs = (((8*t)/300) + 2) * 10^5; %variable heat flux
   else
       qs = 10^6;
   end
   k(1:m) = 7.72;  % constant k & c
   c(1:m) = 1485;  
   
   for j = 1:m
       tvar(j) = ((dx)^2)/((2*(k(j)/(rho*c(j))))*(1 + (sigma*epsl*dx*T(j)^3)/k(j)));% time step is calculated at each time step
   end
   dt = min(tvar);
   for i = 2:m-1 % Governing Eqn. for interior points
       T_new(i) = (((k(i) + 0.25*k(i+1) - 0.25*k(i-1))*T_old(i+1)) + ((k(i) - 0.25*k(i+1) + 0.25*k(i-1))*T_old(i-1)))*(dt/(dx*dx*rho*c(i))) + (1 - ((2*k(i)*dt)/(dx*dx*rho*c(i))))*T_old(i);
   end
   T_new(m) = (T_old(m-1) - T_old(m))*(((k(m) + k(m-1))*0.5*dt)/(rho*c(m)*dx*dx)) + T_old(m);% Governing eqn at backend

   T_new(1) = (( epsl*sigma*(T_inf^4 - T_old(1)^4) + ((((k(1) + k(2))/2)/dx)*(T_old(2) - T_old(1))) + qs ) * ((2*dt)/(rho*c(1)*dx))) + T_old(1);   
   %Governing Eqn at exposed end
   
   T_old = T_new; %update temperature at each time step
   T = T_new;
   time(count) = t;
   Texp(count) = T(1);
   Tback(count) = T(end);
   count = count + 1;
   t = t + dt; 

   if int32(t) == 100  % At 100 sec
      fig1 = T;
   end   
   if int32(t) == 200  % At 200 sec
       fig2 = T;
   end
   if int32(t) == 300  % At 300 sec
       fig3 = T;
   end
   if int32(t) == 450  % At 450 sec
       fig4 = T;
   end
   if int32(t) == 600  % At 600 sec
       fig5 = T;
   end

   
end
end
end

% Results at constant properties
figure(1) %plotting temp. variation at different times
subplot(1,2,1)
sgt1 = sgtitle('Results at Fixed Property for 18.5 cm TPS');
sgt1.FontSize = 20;
plot(100*(0:dx:x),fig1)
axis square
grid on
title('Temperature distribution through the TPS material')
ylabel("Temperature('c)")
xlabel("x coordinate at any y location (cm)")
hold on
grid on
plot(100*(0:dx:x),fig2)
plot(100*(0:dx:x),fig3)
plot(100*(0:dx:x),fig4)
plot(100*(0:dx:x),fig5)
legend('$\tau$ = 100 s','$\tau$ = 200 s','$\tau$ = 300 s','$\tau$ = 450 s','$\tau$ = 600 s','Interpreter','latex')

%plotting temperature at exposed and back end during 10 minutes
subplot(1,2,2)
plot(time,Texp)
hold on
plot(time,Tback)
axis square
grid on
title("Temperature variation with time at the exposed and back ends")
ylabel("Temperature(K)")
xlabel("Time(s)")
legend("Exposed End","Back End")

fprintf('\nMaterial Selected : RCC (reinforced Carbon Carbon Composite)\n \n');
fprintf('Thermal conductivity :%6.2f W/mK\n',k(1));
fprintf('Density :%6.2f Kg/m^3\n',rho);
fprintf('Specific heat :%6.2f J/Kg.k\n',c(1));
fprintf('Average Thermal diffusivity :%.9f m^2/s\n',alpha);
fprintf('Maximum allowable temperature :%6.3f K\n',T_max);
fprintf('Thickness :%6.2f cm\n',100*x);
fprintf('Weight :%6.2f kg/m^2\n',rho*x);


%% Calculations for variable properties
clear
% Least square curve fitting using 3rd order polynomial
X = [144.4 255.6 366.7 477.8 533.3 811.1 1088.9 1366.7 1644.4 1811.1 2200.0 ]; %Temperature in kalvin
Y =  [2.30  3.89  5.05  6.06  6.35  7.36 7.65 7.79  7.65  7.58  7.49 ]; %Thermal conductivity in W/mK
Z =  [5.02 * 10^2 7.12 * 10^2 8.79 * 10^2 1.00 *10^3 1.09 * 10^3 1.30 * 10^3 1.42 * 10^3 1.55 * 10^3 1.67 * 10^3 1.72 * 10^3 1.84 * 10^3]; %Specific heat in J/KgK
p1 = polyfit(X,Y,3); %Polyfit is inbuilt matlab command to approximate polynomial of given degree to fit data using least square method
p2 = polyfit(X,Z,3); % cubic polynomial is approximated

for m =671%:1001; %optimum no. of grid points is obtained here

T_i = 300; % Temperature of TPS at 0 seconds
T_s = 300; % Temperature of base material at 0 second
T_back = 305; % Max allowable temp. at back face of TPS
T_inf = 0; % Space temperature
sigma = 5.67 * 10^-8; % Stefan Boltzmann constant
epsl = 0.9; % Emmisivity can be changed here


for x = 0.204  % x is varied here to obtain optimum thickness of TPS

T(1:m) = T_i; % Initial Temperatures
T_old = T; 

% Reinforced Carbon Carbon (RCC) Composite properties

rho = 1580;     % Constant Density

k = p1(1)*T.^3 + p1(2)*T.^2 + p1(3)*T + p1(4);   %variable k & c
c = p2(1)*T.^3 + p2(2)*T.^2 + p2(3)*T + p2(4);

T_max = 2030; % Melting/max allowable Temperature of RCC

dx = x/(m-1); % Grid size
for j = 1:m
tvar(j) = ((dx)^2)/((2*(k(j)/(rho*c(j))))*(1 + (sigma*epsl*dx*T(j)^3)/k(j))); % calculation of variable time step due to variable properties
end
count = 1;
dt = min(tvar); % Smallest time step is taken for each iteration
t = dt; % calculation of cumulative time

while t < 600 && max(T) < T_max && min(T) < T_back % requirements of given TPS can be changed here
   if t < 300
       qs = (((8*t)/300) + 2) * 10^5; %variable heat flux
   else
       qs = 10^6;
   end

   k = p1(1)*T.^3 + p1(2)*T.^2 + p1(3)*T + p1(4);% k & c are calculated at each iteration
   c = p2(1)*T.^3 + p2(2)*T.^2 + p2(3)*T + p2(4);
   for j = 1:m
       tvar(j) = ((dx)^2)/((2*(k(j)/(rho*c(j))))*(1 + (sigma*epsl*dx*T(j)^3)/k(j)));% time step is calculated at each time step
   end
   dt = min(tvar);
   for i = 2:m-1 % Governing Eqn. for interior points
       T_new(i) = (((k(i) + 0.25*k(i+1) - 0.25*k(i-1))*T_old(i+1)) + ((k(i) - 0.25*k(i+1) + 0.25*k(i-1))*T_old(i-1)))*(dt/(dx*dx*rho*c(i))) + (1 - ((2*k(i)*dt)/(dx*dx*rho*c(i))))*T_old(i);
   end
   T_new(m) = (T_old(m-1) - T_old(m))*(((k(m) + k(m-1))*0.5*dt)/(rho*c(m)*dx*dx)) + T_old(m);% Governing eqn at backend

   T_new(1) = (( epsl*sigma*(T_inf^4 - T_old(1)^4) + ((((k(1) + k(2))/2)/dx)*(T_old(2) - T_old(1))) + qs ) * ((2*dt)/(rho*c(1)*dx))) + T_old(1);   
   %Governing Eqn at exposed end
   T_old = T_new; %update temperature at each time step
   T = T_new;
   time(count) = t;
   Texp(count) = T(1);
   Tback(count) = T(end);
   count = count + 1;
   t = t + dt; 

   if int32(t) == 100  % At 100 sec
      fig1 = T;
   end   
   if int32(t) == 200  % At 200 sec
       fig2 = T;
   end
   if int32(t) == 300  % At 300 sec
       fig3 = T;
   end
   if int32(t) == 450  % At 450 sec
       fig4 = T;
   end
   if int32(t) == 600  % At 600 sec
       fig5 = T;
   end

   
end
end
end

% Results at variable properties
figure(2) %plotting temp. variation at different times
subplot(1,2,1)
sgt2 = sgtitle('Results at Variable Property for 20.4 cm TPS');
sgt2.FontSize = 20;
plot(100*(0:dx:x),fig1)
axis square
grid on
title('Temperature distribution through the TPS material')
ylabel("Temperature('c)")
xlabel("x coordinate at any y location (cm)")
hold on
grid on
plot(100*(0:dx:x),fig2)
plot(100*(0:dx:x),fig3)
plot(100*(0:dx:x),fig4)
plot(100*(0:dx:x),fig5)
legend('$\tau$ = 100 s','$\tau$ = 200 s','$\tau$ = 300 s','$\tau$ = 450 s','$\tau$ = 600 s','Interpreter','latex')

%plotting temperature at exposed and back end during 10 minutes
subplot(1,2,2)
plot(time,Texp)
hold on
plot(time,Tback)
axis square
grid on
title("Temperature variation with time at the exposed and back ends")
ylabel("Temperature(K)")
xlabel("Time(s)")
legend("Exposed End","Back End")


