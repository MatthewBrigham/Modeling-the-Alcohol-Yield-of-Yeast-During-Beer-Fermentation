function Ferment1(T, sub_mass, cell_conc, yeast_mass, V)
%Runge Kutta Attempt
%function Fermentation(T, W, cellconc.)
%
%Background:
%This function considers the number of yeast cells of S. cerevisiae
%inoculated into a Beer fermenting system. The function will output the
%evolution of the yeast population, the consumption of the sugar substrate (table sugar),
%and the production of ethanol with respect to time. Many calculations are
%made using the Euler approximation for first order differential equations
%
%Assumptions:
%(1)the substrate has been heated up to form the mash (wort) and has cooled down to the desired temperature of fermentation.
%(2)the yeast cells are either dead or active, none are in the lag phase. 
%(3)the yeast cells used are brand new (they havent been sitting for a long (time)
%
%Inputs:
%T: is the temperature in Celsius
%sub_mass: is the weight of the sugar used for fermentation in lbs. 
%cell_conc: the the approximate number of cells per gram of yeast powder
%yeast_mass: the number of grams of yeast added to the fermentor
%V is volume in Liters of water used in the fermentation. This value should
%be slected in determining how much beer is desired

T = T + 273.15;                     %convert input temperature of Celsius to Kelvin
N = cell_conc*yeast_mass/V;         %Number of inoculated cells-Computes the approximate number of active yeast cells inoculated
S0 = 130; %sub_mass*453.592/V;              %Initial Substrate Conc. Converts lbs to grams/L of initial substrate

mu_x0 = exp(108.31-31934.09/T);     %max cell growth rate
mu_s0 = exp(-41.92+11654.64/T);     %max substrate consumption rate (achieved at substrate saturation
mu_lag = exp(30.72-9501.54/T);      %specific cell activation rate
mu_DT = exp(130.16-38313/T);        %specific cell death rate
mu_e0 = exp(3.27-1267.24/T);        %max EtOH production rate
k_e = exp(-119.63+34203.95/T);      %affinity constant
k_s = k_e;                          %affinity constant
k_x = k_e;

deltat = 0.01;
t = 0:3:200;                    %time interval of fermentation 0 to 4 days in units of hours with data taken at every half hour
X_act = zeros(1, length(t));        %vector for number of active cells at a given time, t
X_act(1) = 0.02*N;                       %initial number of acitve yeast cells 

C_s = zeros(1, length(t));
C_s(1) = 130; %S0/V;                      %initial concentration of sugar substratein g/L
mu_s = zeros(1, length(t));
mu_s(1) = mu_s0*C_s(1)/(k_s + C_s(1));

C_e = zeros(1, length(t));
C_e(1) = 0;
mu_e = zeros(1, length(t));
mu_e(1) = mu_e0*C_s(1)/(k_e+C_s(1));
%keyboard
mu_x = zeros(1, length(t));
mu_x(1) = mu_x0*C_s(1)/(k_x+C_e(1));

f = zeros(1, length(t));
f(1) = 1-C_e(1)/(0.5*S0);

X_lag = zeros(1, length(t));
X_lag(1) = 0.98*N;

%keyboard
for n = 2:length(t)-2                           
    X_act(n) = X_act(n-1) + deltat.*(X_act(n-1).*(mu_x(n-1)-mu_DT) + X_lag(n)*mu_lag);
    X_lag(n) = mu_lag*(N-X_act(n));
    C_s(n) = C_s(n-1) - deltat*(mu_s(n-1)*X_act(n-1));
    C_e(n) = C_e(n-1) + deltat*(f(n-1)*mu_e(n-1)*X_act(n-1));
    
    mu_x(n) = mu_x0*C_s(n)/(k_x+C_s(n));        
    mu_s(n) = mu_s0*C_s(n)/(k_s+C_s(n));
    mu_e(n) = mu_e0*C_s(n)/(k_e+C_s(n));
    f(n) = 1-C_e(n)/(0.5*S0);                      
    
    if C_s(n)<0
       C_s(n) = 0;
    end
end
%keyboard

figure(1);
plot(t(1:n-2), 0.01*X_act(1:n-2));
title('Cell Biomass');
xlabel('t [hrs]');
ylabel('X_act [g/L]');

figure(2);
plot(t(1:n-2), C_s(1:n-2));
title('Sugar Concentration');
xlabel('t [hrs]');
ylabel('C_s [g/L]');

figure(3);
plot(t(1:n-2), C_e(1:n-2));
title('Ethanol Concentration');
xlabel('t [hrs]');
ylabel('C_e [g/L]');

%keyboard
