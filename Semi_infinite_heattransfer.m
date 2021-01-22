clear, clc
%Pg.339 of Fundamental of Heat Transfer 7th ed.
%% parameters
%Air properties
V = 33454.7; %ft^3
V = V*0.0283168; %m^3
Cp = 1006; %J/kg*K
rho = 1.225; %0.6601; %@6000 kg/m^3
V_s = 0.1012799;% m^3/s 
m = V*rho; %kg
m_dot = (V_s*rho); % kg/s
diff_air = 20*10^-6; %m^2/s
hair = 12; %turbulent air W/m*K
%Tube properties
L_tube = 153;
diff_PEH = 2.77*10^-7; %m2/s
k_PEH = .5; %W/m*K
r1 = 0.1541;%m inner r
r2 = 0.1683;%m outer r 
A_tube =2*pi*r1; %m^2
t = 3600; %1hr in seconds
SA_tubes  = r2*2*pi*L_tube;
%Dirt properties
diff_soil = .91*10^-6; %m2/s
k_soil = 2.1; %W/m*K
%Initial Condition
T_G = -12; %C Inlet air temp 
T_EarthC = 6; %C Constant surface earth temperature 
%Thermal resistances
R_conv = (1/hair)/A_tube; %K/W
R_pipe = (log(r2/r1))/(2*pi*k_PEH*L_tube);
R_contact =1-(6/100);
R_T = R_contact*(R_conv+R_pipe);

%% EQ's
%METHOD 1 - FLUX into pipe -> Q into pipe -> Temp out factoring in thermal
%resistances  ** NOT WORKING, temps go out of possible range
Flux_Earth = (k_soil*(T_G - T_EarthC))/(2*sqrt(pi*diff_PEH*t));
Q_Earth = Flux_Earth*SA_tubes*R_T;
T_delta = Q_Earth *R_T;
T_out1 = T_G + T_delta

%METHOD 2 - T out based on equation pg.339. 
%**NOT WORKING, does not factor in thermal resistance between
%earth/tubes/air. Way too effective, reaches earth temp within 1.5 meters
T_out2 = (erf(L_tube/(2*sqrt(diff_PEH*t)))*(T_EarthC- T_G)) + T_G