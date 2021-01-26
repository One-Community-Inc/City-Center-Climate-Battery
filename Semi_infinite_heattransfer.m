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
diff_air = 20;%*10^-6; %m^2/s
hair = 12; %turbulent air W/m*K
%Tube properties
L_tube = 1.53;
diff_PEH = 2.77*10^-7; %m2/s
k_PEH = .5; %W/m*K
r1 = 0.1541;%m inner r
r2 = 0.1683;%m outer r 
A_tube =2*pi*r1; %m^2
t = 3600; %1hr in seconds
SA_tubes  = r2*2*pi*L_tube;
%Dirt properties
diff_soil = .91;%10^-6; %m2/s
k_soil = 2.1; %W/m*K
%Initial Condition
T_G = -12+273; %C Inlet air temp 
T_EarthC = 6+273; %C Constant surface earth temperature 
%Thermal resistances
R_conv = (1/hair)/A_tube; %K/W
R_pipe = (log(r2/r1))/(2*pi*k_PEH*L_tube);
R_contact =1-(6/100);
R_T = R_contact*(R_conv+R_pipe);

%% EQ's
%METHOD 1 - FLUX into pipe -> Q into pipe -> Temp out factoring in thermal
%resistances  ** NOT WORKING, temps go out of possible range
Flux_Earth = (k_soil*(T_G - T_EarthC))/(sqrt(pi*diff_PEH*t)); %W/m^2
Q_Earth = Flux_Earth*SA_tubes;%W
T_delta1 = Q_Earth *R_T;%K
T_out1 = T_G + T_delta1-273%K
%METHOD 2 - T out based on equation pg.339. 
%**NOT WORKING, does not factor in thermal resistance between
%earth/tubes/air. Way too effective, reaches earth temp within 1.5 meters
T_out2 = (erf(L_tube/(2*sqrt(diff_PEH*t)))*(T_EarthC- T_G)) + T_G;
T_out2 = T_out2-273
%METHOD 3 - Row of pipes in semi-infinite solid 2.1.19
%pg.75 in Conduction Heat Transfer Solutions
s = 0.3048; %12in ->m
d = 3.9624; %13feet 
Bi_1 = (hair*r1)/k_PEH;
Bi_2 = (k_soil*d)/k_soil; %??? h2? d?
t2 = T_EarthC;
t1 = T_G;
D = d/s;
q = 2*pi*k_soil*(t2-t1)/((1/Bi_1)+log((d/(pi*r1*D))*sinh(2*pi*(D+(D/Bi_2)))));
T_delta3 = q *R_T*SA_tubes;
T_out3 = T_G + T_delta3 -273