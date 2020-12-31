%Greenhouse Temp (No Heating)
%HELP - Line 137 Unit problem?
%useful source - https://www.engineeringtoolbox.com/heat-loss-buildings-d_113.html
%Includes infliration, conduction and radiation heat sources
%Using Double-Poly material data as an equivalent for SolaWrap
%Page 161-162 on "Climate battery details" for material data. (NEED SOURCE)
clear; clc;
%% Parameters 
A_earth = 1300;   %ft^2 
A_walls = (330.31*2) + 395.74 + 632.38; %ft^2 
A_roof = 2987.80;   %ft^2 
V = 33454.7; %ft^3
V_m = V*0.0283168; %m^3
A_r = 259.75;%m^2 Area of roof
t_roof = 85; %roof angle ** 70 is optimal for solar but too steep
t_s = .83; %transparency
%T_desiredw = 60; %F
%T_desireds = 90;% F

%Sun Angle
t_winter = 26.5; %degrees
t_summer = 73.5; %degrees
%Air properties -https://www.engineeringtoolbox.com/standard-atmosphere-d_604.html
Cp = 1.006; %kJ/kg*K
rho = 0.6601; %@6000 kg/m^3
m = V_m*rho; %kg
%Imported Data 
T_a = readtable('Circleville.txt');
%TT = timetable(Date,Station Id,Air Temperature Observed (degF), ...
%Solar Radiation Average (watt/m2),Wind Speed Average (mph))
TT = T_a(1:24,1:6);
T_d = TT{1:24,3};
S_sun = TT{1:24,4};
winds = TT{1:24,5};
%T_e = TT{1:24,6}; *using analytical solution average
%Concrete Properties
U_c = .35; %-depends on thickness/type
%Double-Poly properties
ACH = .5; %ACH due to inflitration
m_dot = (m*ACH)/(60*60); % kg/s
Uwind = [[0,5,10,20,25,30];[.535,.631,.675,.716,.728,.736]]; %U due to wind (Btu/h*ft^2*F)
%pre-allocate
U = zeros(1,24);
Delta_T = zeros(1,24);
Delta_Te = zeros(1,24);
for i = 1:length(winds)
    if winds(i)<5 
        U(i) = Uwind(2,1); 
    elseif winds(i)< 10
        U(i) = Uwind(2,2); 
    elseif winds(i)< 20
        U(i) = Uwind(2,3);
    else
        U(i) = Uwind(2,4);
    end
end
%% Diana's 
%Define Variables
a = 0.055; %m^2*s^-1 for loamy soil 
Tinf = [9, 9, 9, 9,8,9,8,8,7,10,13,16,17,19,20,20,18,14,12,11,9,7,4,2]; % SCAN data for air temperature
B = [0.04,0.03,0.03,0.04,0.03,0.04,0.03,2.05,90.80,279.60,424.20,520.20,547.30,517.70,427.10,286.50,112.60,3.65,0.01,0.01,0.01,0.02,0.04,0.12]; %W*m^-2SCAN Data
a0 = 0.6; % percentage off solar radiation absorbed
k = 1.169 ; % W*(m*K)^-1
e = 0.85; %emissivity of soil

% Convert Temp to Kelvin
C = @f_to_c;
Tinf = C(Tinf); % C
% Solve for the temperature of the sky - Tsky
Tsky = Tinf - 12; % C

% Stefan-Boltzman constant: sigma
sigma = 5.670374419*10^-8; % W*(m^-2)*(K^-4)
deltaR = sigma.*((Tinf+273.15).^4-(Tsky+273.15).^4); % temp in Kelvin

Tinf = Tinf +273.15; %K

% Ti is the temperature the soil reaches as you go deeper
Ti = 43; % degree Farenheit - needs to be in K for equation
Ti = C(Ti) +273.15; % K


% Air Velocity
v = [2.4, 2.5,2,3.3,5.1,1.6,4.6,5.1,2.9,7.5,5,9.6,8.1,7.8,11.6,14,14.4,12.2,12.4,9.8,3.6,3.5,2.1,2.4];
v = v*1609.34/(60*60); % convert to m/s

%calculate heat convection coefficient
hc = 2.8+3.*v;
hr = 4*e*sigma.*Tinf.^3;
h = hc + hr;

% Depth 4in
x_i = 0.01; %0-156in in meters(3.9624)
x_f = 3.96;
t = [1:24]; %hours
t = t*60*60; %convert time to seconds
%% Convert from Farenheit to Celsius
F = @c_to_f;
T_i = analyticalsol(x_i,t,h,a,Tinf,B,a0,k,Ti,e,deltaR);
T_i = F(T_i-273.15);
T_f = analyticalsol(x_f,t,h,a,Tinf,B,a0,k,Ti,e,deltaR);
T_f = F(T_f-273.15);
T_e = (T_i+T_f)./2;
%% Equations
%pre-allocate
q_c = zeros(1,24);
q_i = zeros(1,24);
S = zeros(1,24);
Q_solar = zeros(1,24);
T2 = zeros(1,24);
Q_i_cW = zeros(1,24);
Q_T = zeros(1,24);
for i = 1:length(S_sun)
    
    if i ==1
        Delta_T(i) = 0;
        Delta_Te(i) = 0;
        T2(i) = T_d(i);%F intial temp T_d(1);
    else
  %Conduction
        Delta_Te(i) = T_e(i-1)- T2(i-1); %F 
        Delta_T(i) = T_d(i-1,1)-T2(i-1); %F 
        q_c(i) = A_walls*Delta_T(i)*U(i)+A_roof*Delta_T(i)*U(i)*1.15 +... 
        A_earth*U_c*Delta_Te(i); %Btu/h 
        %Inflitration
        %(https://soa.utexas.edu/heat-loss-due-infiltration-using-air-change-method#:~:text=Air%2DChange%20Method-,OVERVIEW,cracks%20around%20windows%20and%20doors.&text=Infiltration%20is%20caused%20by%20wind,movement%20within%20the%20building%20envelope.)
        q_i(i) = ACH*V*Delta_T(i)*.0018; %Btu/hr     
        %Radiation
        S(i) = S_sun(i)*cosd(t_roof-t_winter)*t_s; %W/m^2 Q due to solar
        Q_solar(i) = (S(i)*A_r)/1000; %kW
        
        %Totals & Btu/hr -> kW
        Q_i_cW(i) = (q_c(i)+q_i(i))*(0.29307107/1000); %Conversion to kW
        Q_T(i) = Q_i_cW(i) + Q_solar(i); %kW
    
        %Solving for T(F)         eq.Q = (m*Cp*(T2-T1)) kW = kg/s*kJ/kg/C *C
        T2(i) = c_to_f((f_to_c(T_d(i-1))+(Q_T(i)/(m_dot*Cp)))); %
    end
end

plot(1:24,T2)
xlabel('Hour')
ylabel('Temperature(F)')
title('Greenhouse Temperature Jan 1st w/o heating')
xlim([1 24])
hold on 
plot(1:24,T_d)
legend('Inside','Outside')
legend
function Celsius = f_to_c(Fahrenheit)
Celsius = (5/9)*(Fahrenheit-32);
end
function Fahrenheit = c_to_f(Celsius)
Fahrenheit = (Celsius/(5/9))+32;
end
function T = analyticalsol(x,t,h,a,Tinf,B,a0,k,Ti,e,deltaR)
T =((-sqrt(pi.*a.*t).*(h.*(Ti-Tinf)-B.*a0+e.*deltaR))./(2.*h.*sqrt(pi.*a.*t)+k)).*(erfc((x./sqrt(4.*a.*t)))) + Ti;
end
%% Convert from Farenheit to Celsius


%%%OLD VARIABLES
%B = .75;%heat recovery efficiency (%) rough figure for geothermal
%R = 1.75; %SolaWrap
%U = 1/R;
%Delta_T=  12-60; %F *summer 98-7 
%Humidity range: 60% - 90%;
%Solar
% S_full = 1000; % W/m^2 full sun
% S_light = 850; % W/m^2 light clouds
% S_cloudy = 150; % W/m^2 thick clouds
% S_sun = [S_full S_light S_cloudy];

