%Greenhouse:Heating
clear; clc;
%%Assumptions
%1. Uniformly mixed air inside the greenhouse. Fair assumptions with
%circulation fans
%2. Contact thermal resistance = 0. No data for chosen materials.
%3. Climate battery/geothermal subsystems have equivalent air velocity 
%Plastic-Earth contact resistance
%% Parameters 
%Greenhouse
T_desired = f_to_c(60); %C
A_earth = ft2_to_m2(1300); %m^2
A_walls = ft2_to_m2((330.31*2) + 395.74 + 632.38); %m^2
A_roof = ft2_to_m2(2987.80); %m^2
V = 33454.7; %ft^3
V = V*0.0283168; %m^3
A_r = 259.75;%m^2 Area of roof
t_roof = 85; %roof angle 
% Geothermal/ClimateBattery System as of 1/6/21
r_E = 6.315:10:288; %radi of dirt
CFM = 214.6*2; %ft^3/min per subsystem (2subsystems)
V_s = 0.1012799;% m^3/s
%Main Tubing (6in)
r1 = i_to_m(6.065); %in->m pipe ID
r2 = i_to_m(6.625); %in->m pipe OD
A_tube =2*pi*r1; %m^2
L_tube = i_to_m((2*252.58)*12); %in->m 2 systems ***doublecheck 
v_tube = V_s/(7*A_tube); %m/s 7 tubes per system
%Manifold Tubing (12in)
% r1_M = i_to_m(11.89);%in->m 
% r2_M = i_to_m(12.75);%in->m
% A_Manifold = 2*pi*r1_M; %m^2
% L_mainfold = i_to_m((61.1*2)*12); %ft
% v_manifold = V_s/(A_Manifold);
%Sun Angle
t_winter = 26.5; %degrees
%Air properties -https://www.engineeringtoolbox.com/standard-atmosphere-d_604.html
Cp = 1006; %J/kg*K
rho = 0.6601; %@6000 kg/m^3
m = V*rho; %kg
U_c = .3/.152; %6in concrete(.152m) W/m2K 
m_dot = (V_s*rho); % kg/s
m_dot_infl = (.5*V*rho)/(60*60);%Estimated to be .5V/hr for most buildings
%Double-Poly properties
t_s = .83; %transparency of SolarWrap/double poly
Uwind = [[0,5,10,20,25,30];[.535,.631,.675,.716,.728,.736]*5.6745]; %U-doubly poly due to wind (Btu/h*ft^2*F)->W/m^2K
%Imported Data
T_a = readtable('Circleville.txt');%(Date,Station Id,Air Temperature Observed (degF), ...
%Solar Radiation Average (watt/m2),Wind Speed Average (mph))
TT = T_a(1:24,1:10);
T_air = TT{1:24,3}';
S_sun = TT{1:24,4}';
winds = TT{1:24,5}';
T_Earth = TT{1:24,10}';
T_EarthC = f_to_c(T_Earth);
T_airC = f_to_c(T_air);
%pre-allocate
h_c = zeros(1,24);
for i = 1:length(winds) %Heat transfer coefficent based on Wind 
    if winds(i)<5 
        h_c(i) = Uwind(2,1); 
    elseif winds(i)< 10
        h_c(i) = Uwind(2,2); 
    elseif winds(i)< 20
        h_c(i) = Uwind(2,3);
    else
        h_c(i) = Uwind(2,4);
    end
end
%Conductivity
%https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4813881/
hair = 12.12-1.16*v_tube+11.6*v_tube^1/2; %W/m^2*K FORCED AIR can be approx this eq. BETTER TO FIND REAL VALUE VIA NUSSELT#
k_PEH = .5; %W/m*K
k_sand_dry = .4; %W/m*K
k_soil = 2.1; %W/m*K
k_clay_moist = 1.7;%W/m*K
%Diffusivity 
diff_s = .91*10^-6;% m^2/s *loam type soil
diff_clay_moist = .5*10^-6;% m^2/s 
diff_sand_dry = .3*10^-6;% m^2/s
%Pre-allocate
Delta_T = zeros(1,24);
Delta_Te  = zeros(1,24);
Q_Conv = zeros(1,24);
Q_Inf = zeros(1,24);
S = zeros(1,24);
Q_Solar = zeros(1,24);
T_out = zeros(1,24);
Q_CB = zeros(1,24);
Q_Heater = zeros(1,24);
T_G = zeros(1,24);

%% Heat Equations 
%put in for loop for every hour
tspan =[0 60*60]; %1hr in seconds
%dp = 2.3*sqrt(diff_s*t^2);%pentration depth defined as 90% T_i into medium
for i = 1:24
    if i == 1
        T_G(1) = T_airC(1);
    else
        T_G0 = T_G(i-1);
        %Convection
        Delta_T(i) = T_G(i)-T_airC(i);%C
        Delta_Te(i) = T_G(i) - T_EarthC(i); %C
        Q_Conv(i) = A_walls*Delta_T(i)*h_c(i)+A_roof*Delta_T(i)*h_c(i)+ A_earth*U_c*Delta_Te(i); %W
        %Q_Conv(i) = Btuhr_to_W(Q_Conv(i)); %W
        
        %Inflitration
        %(https://soa.utexas.edu/heat-loss-due-infiltration-using-air-change-method#:~:text=Air%2DChange%20Method-,OVERVIEW,cracks%20around%20windows%20and%20doors.&text=Infiltration%20is%20caused%20by%20wind,movement%20within%20the%20building%20envelope.)
        %Q_Inf(i) = ACH_infl*V*Delta_T(i)*.0018; %Btu/hr 
        Q_Inf(i) = m_dot_infl*Cp*(T_G(i)-T_air(i)); %W
        %Q_Inf(i) = Btuhr_to_W(Q_Inf(i));
        
        %Radiation
        S(i) = S_sun(i)*cosd(t_roof-t_winter)*t_s; %W/m^2 Q due to solar
        Q_Solar(i) = (S(i)*A_r); %W
        
        %Climate Battery/Geo
        %T_out(i) = T_G(i)-T_EarthC(i)*exp(((-R_T/(m_dot*Cp))*L_tube)+T_EarthC(i));
        %Q_CB(i) = k_soil*(T_EarthC(i)-T_G(i))/sqrt(pi*diff_s*t);
        %[t,L_tube]=ode45(@(t,L_tube) (erf(L_tube/(2*sqrt(diff_s*t)))*(T_G(i)-T_EarthC(i))) + T_EarthC(i),tspan,L0);
        %T_out(i) = (erf(L_tube/(2*sqrt(diff_s*t)))*(T_G(i)-T_EarthC(i))) + T_EarthC(i);
        %Q_CB = m_dot*Cp*(T_out-T_G);
        
        %Heater
        %Q_Heater(i) = Cp*m_dot*(T_desired-T_G(i));
        
        %Greenhouse T
        %T_G(i) = (Q_CB(i) + Q_Heater(i) + Q_Solar(i) - Q_Conv(i) - Q_Inf(i))/(rho*V_s*Cp);
        %[t,T_G] = ode45(@(t,T_G) (m_dot*Cp*(((erf(L_tube/(2*sqrt(diff_s*t)))*(T_G-T_EarthC(i))) + T_EarthC(i))-T_G))...
         %   + (Q_Heater(i) + Q_Solar(i) - Q_Conv(i) - Q_Inf(i))/(rho*V_s*Cp),tspan,T_G0);
        [t,T_G] = ode45(@(t,T_G) (m_dot*Cp*(((erf(L_tube/(2*sqrt(diff_s*t)))*(T_G-T_EarthC(i))) + T_EarthC(i))-T_G))...
            + (Q_Heater(i) + Q_Solar(i) - Q_Conv(i) - Q_Inf(i))/(rho*V_s*Cp),tspan,T_G0);
        hold on
        plot(t,T_G)
        legend
    end 
end
%plot(1:24,T_G)
%% Functions
%Unit conversion
function Celsius = f_to_c(Fahrenheit)
Celsius = (5/9)*(Fahrenheit-32);
end
function Fahrenheit = c_to_f(Celsius)
Fahrenheit = (Celsius/(5/9))+32;
end
function Meters = i_to_m(In)
Meters = In*.0254;
end
function Watts = Btuhr_to_W(Btuhr)
Watts = 0.29307107*Btuhr;
end
function Meters2 = ft2_to_m2(ft2)
Meters2 = ft2*0.092903;
end
%Unused
%Thermal Resistance C/W
% R_conv = 1/hair*A_tube; 
% R_pipe = (log(r2/r1))/(2*pi*k_PEH*L_tube);
% R_earth = (log(r_E/r2))/(2*pi*k_soil*L_tube);
% R_contact =1-(6/100); %Contact resistance ~6% with grout- https://info.ornl.gov/sites/publications/Files/Pub57560.pdf
% R_T = R_contact*(R_conv+R_pipe+R_earth);  