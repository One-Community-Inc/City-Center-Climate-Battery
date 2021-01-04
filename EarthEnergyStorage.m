%EarthEnergyStorage
%Underground Pipe Pressure -https://www.engineeringtoolbox.com/underground-pipe-pressure-soil-transport-d_2145.html
%Insulated Pipe Heat Transfer - https://www.engineersedge.com/heat_transfer/heat_loss_insulated_pipe_13865.htm
%%Parameters
%Temps/Conductivity
Tair_in = f_to_c(100);%F->C
Tearth_inf = f_to_c(55);%F->C
v = 0.635;% m/s
hair = 10.45 - v + 10*v^1/2; %W/m*K FORCED AIR can be approx this eq. BETTER TO FIND REAL VALUE VIA NUSSELT#
k_PEH = .5; %W/m*K
k_sand_dry = .20; %W/m*K
k_soil = .835; %W/m*K
k_clay_moist = 1.17;%W/m*K
%Geometry
r1 = i_to_m(6); %in->m pipe ID
r2 = i_to_m(6.315); %pipe OD
r3 = 6.315:10:288; %radi of dirt
A =2*pi*r1;
L = i_to_m(12*12);%203.56 Taken from Climate Battery spreadsheet(SS1)
%pressure on pipes
p_s = 1850; %Density Soil kg/m^3
p_sw = 1100; %Density soil below groundwater level
g = 9.8;%m/s
h = i_to_m(15*12 +4); %in
h_w = i_to_m(87.36*12); %Depth from groundwater level to object ... 
%For BEAVER COUNTY -> https://waterdata.usgs.gov/ut/nwis/current/?type=gw
psoil = p_s*g*(h-hw)+p_sw*g*h_w;
%Resistance C/W
R_conv = 1/hair*A; 
R_pipe = (log(r2/r1))/(2*pi*k_PEH*L);
R_earth = (log(r3/r2))/(2*pi*k_soil*L);
R_total = R_conv+R_pipe+R_earth;
%Heat EQ
Q = (Tair_in-Tearth_inf)./(R_total); %W
T_E = c_to_f(Tair_in-(Q.*R_earth));
%Graph
plot(r3/12,T_E)
title('Temp Around Climate Battery')
xlabel('Distance(ft)')
ylabel('Temperature (F)')
function Celsius = f_to_c(Fahrenheit)
Celsius = (5/9)*(Fahrenheit-32);
end
function Fahrenheit = c_to_f(Celsius)
Fahrenheit = (Celsius/(5/9))+32;
end
function Meters = i_to_m(In)
Meters = In*.0254;
end