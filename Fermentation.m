function Fermentation(N, S0, V)

%function Fermentation(N, S0, V)
%
%Background:
%This function considers the number of Brewers yeast cells (mostly S.
%cervisiae inoculated into a Beer fermenting system. The function will 
%output the evolution of the yeast population, the consumption of the 
%sugar substrate, and the production of ethanol with respect 
%to time. Many calculations are made using the Euler approximation for
%interdependent first order differential equations. The optimal time to
%stop the fementation is outputted as an fprintf function.
%
%Assumptions:
%(1)The substrate has been heated up to form the mash (wort) and has cooled 
%   down to the desired temperature of fermentation.
%(2)Temperature is held constant at room temperature. 
%
%Inputs:
%N: number of yeast cells inoculated
%S0: initial sugar substrate concentration (good values are between 80 and 200)
%V: volume in Liters of water used in the fermentation (good values are
%in the range of 20L {about 5 gallons})


%============================== Calculate the starting parameters

T = 24 + 273.15;                         %convert Celsius to Kelvin

mu_x0 = exp(108.31-31934.09/T);         %max cell growth rate
mu_s0 = exp(-41.92+11654.64/T);         %max substrate consumption rate (achieved at substrate saturation
mu_lag = exp(30.72-9501.54/T);          %specific cell activation rate
mu_DT = exp(130.16-38313/T);            %specific cell death rate
mu_e0 = exp(3.27-1267.24/T);            %max EtOH production rate
k_e = exp(-119.63+34203.95/T);          %affinity constant
k_s = k_e;                              %affinity constant
k_x = k_e;

deltat = 0.01;
t = 0:3:200;                            %time interval of fermentation

C_s = zeros(1, length(t));
C_s(1) = S0;                            %initial concentration of sugar substratein g/L
mu_s = zeros(1, length(t));
mu_s(1) = mu_s0*C_s(1)/(k_s + C_s(1));

C_e = zeros(1, length(t));
C_e(1) = 0;
mu_e = zeros(1, length(t));
mu_e(1) = mu_e0*C_s(1)/(k_e+C_s(1));

X_act = zeros(1, length(t));            %vector for number of active cells at a given time, t
X_act(1) = 0.05*N;                      %initial number of acitve yeast cells 

mu_x = zeros(1, length(t));
mu_x(1) = mu_x0*C_s(1)/(k_x+C_e(1)); 

f = zeros(1, length(t));
f(1) = 1-C_e(1)/(0.5*S0);

X_lag = zeros(1, length(t));
X_lag(1) = 0.95*N;

%============================== Calculate the concentrations and rates with the Euler Method

for n = 2:length(t)-2                           
    X_act(n) = X_act(n-1) + deltat.*(X_act(n-1).*(mu_x(n-1)-mu_DT)+ X_lag(n)*mu_lag);
    X_lag(n) = mu_lag*(N-X_act(n));
    C_s(n) = C_s(n-1) - deltat.*(mu_s(n-1).*X_act(n-1));
    C_e(n) = C_e(n-1) + deltat.*(f(n-1).*mu_e(n-1).*X_act(n-1));
    
    mu_x(n) = mu_x0.*C_s(n)./(k_x+C_s(n));       
    mu_s(n) = mu_s0.*C_s(n)./(k_s+C_s(n));
    mu_e(n) = mu_e0.*C_s(n)./(k_e+C_s(n));
    f(n) = 1-C_e(n)./(0.5*S0);                      
    
    if C_s(n)<0
       C_s(n) = 0;
    end
end



X_act(1) = N; 
X_act(2:length(t)-2) = X_act(2:length(t)-2)*0.004;
figure(1);
plot(t(1:n-2), X_act(1:n-2));
title('Cell Biomass for T = %d C');
xlabel('t [hrs]');
ylabel('Active Cells [g/L]');

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


%============================== Non-recursive probability calculations
    
    %*************************************************************************%
    %Calculating W, the initial cell mass distribution of the inoculated cells%
    %*************************************************************************%
    m = 0:0.01E-12:1E-11;
    epsilon = 3*sqrt(2)*10^-13;
    t = 0;
    W = 0.1.*(m./epsilon).^5.*exp(-m./epsilon).*10.^24;

    figure(4);
    plot(m.*1E12, W);
    title(sprintf('Cell Mass Distribution at a t = %d sec', t))
    xlabel('mass [picograms]')
    ylabel('number of yeast cells')

    %*****************************************************************************************%
    %Calculating the Transition probability (Gamma) and the growth rate of the individual cell%
    %*****************************************************************************************%
    
    R = 5E-6;       %meters
    rho = 1.01;     %g/cm^3
    K_sL = 0.02;    %g/L
    mu_cL = 1;      %hr^-1
    m_cL = 3E-12;   %g
    phi = 6E-5;     %g/cm^2?hr
    C_s = 0.03;     %g/L

    rprime = 2/(R*rho)*phi*C_s.*m./(K_sL+C_s);  %Calculates r' (Equation 12)
    rdubprime = mu_cL.*m;                       %Calculates r" (Equation 13)
    r = rprime - rdubprime;                     %Calculates r from r' and r" (Equation 11)


    L = 2/(epsilon*sqrt(pi)).*exp(-((m-m_cL)./epsilon).^2)./erfc((m-m_cL)./epsilon).*r; %Transition probability (Equation 10)

    figure(5);
    plot(m, L);
    title('Transition Probability')
    xlabel('cell mass')
    ylabel('Gamma(m)')
    
    %*************************************************************%
    %Calculate the Mass distribution denstiy of the daughter cells%
    %*************************************************************%
    
    m_new = input('Enter arbitrary parent mass to determine daughter cell mass probabilities: ');
    p = 30.*m.^2.*(m_new - m).^2./(m_new).^5;

    figure(6)
    plot(m, p);
    title(sprintf('Probability of Daughter Cell Having mass, m, with Parent Having mass, %d', m_new))
    xlabel('mass of daughter cell (grams)')
    ylabel('p(m, m_new)')

%============================== Determine when to stop the fermentation

t = 0:3:200;
t0 = length(t)-10;                              %last time points before the end of the fermentation
tf = length(t)-2;

final_t = t(t0:tf);                             %vectorized for input into polyfit
final_C_e = C_e(t0:tf);                         %vectorized for input into polyfit

extrapolation = polyfit(final_t, final_C_e, 1); %last ten values of the fermentation

stop_ferment = 0.95*extrapolation(2);           %To maximize efficiency, we will suggest stopping the fermentation once it is 95% of the expected ethanol to be produced is achieved

for n = 1:length(C_e)                           %Determines the indices for the time at which the ethanol concentration is close to the 95%
    TOL = 1;                                
    
    if abs(C_e(n)-stop_ferment) < TOL
        break
    end
    
end

t = 0:3:200;

fprintf('To achieve a 95 percent ethanol yield and maximize efficiency, \nstop the fermentation at %d\n hours', t(n))






