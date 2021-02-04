clear, clc
%Pg.339 of Fundamental of Heat Transfer 7th ed.
%% parameters
T_a = readtable('Circleville2019.txt');
TT = T_a(:,1:10);
%T_air = TT{:,3}';
T_air = f_to_c(90);%C
T_Earth = TT{:,10}';
TE_slist = T_Earth(1,2880:6552);% months 5-9 % months 4, april and october excluded
TE_slist(141) = TE_slist(140); %141 is a NAN
TE_wlist = T_Earth(1,[7296:end 1:2160]);%months 11-3
TE_s = f_to_c(sum(TE_slist)/length(TE_slist));%C
TE_w = f_to_c(sum(TE_wlist)/length(TE_wlist));%C
%Tube properties
L_tubes = 153; %m
L_tube = 9.75; %m  13ft
diff_PEH = 2.77*10^-7; %m2/s
k_PEH = .5; %W/m*K
r1 = 0.1541;%m inner r
r2 = 0.1683;%m outer r 
A_tube =2*pi*r1; %m^2
%t = 3600; %1hr in seconds
SA_tubes  = r2*2*pi*L_tubes;%total all tube length
SA_tube = r2*2*pi*L_tube; %1 section of tube
burial_D = 3.9624;% 13 feet
%Air properties
%https://www.engineeringtoolbox.com/convective-heat-transfer-d_430.html
V = 33454.7; %ft^3
V = V*0.0283168; %m^3
Cp = 1006; %J/kg*K
rho = 1.225; %0.6601; %@6000 kg/m^3
V_s = 0.1012799;% m^3/s 
v = V_s/A_tube;
hair = 10.45 - v + 10*v^(1/2);      
%hair = 12; %turbulent air W/m*K
%Dirt properties
diff_soil = .91*10^-6; %m2/s
k_soil = 2.1; %W/m*K 
%Initial Condition
T_Gs = f_to_c(90); %C Inlet air temp 
T_Earths = TE_s; %C earth temperature summer
T_Gw = f_to_c(60); %C
T_Earthw = TE_w; %C earth temperature winter
%Thermal resistances
%R = L/k
R_conv = (1/hair)/A_tube; %K/W
R_pipe = (log(r2/r1))/(2*pi*k_PEH*L_tube);
R_contact =1-(6/100);
R_T = R_contact*(R_conv+R_pipe);

%% EQ's  
%METHOD - Row of pipes in semi-infinite solid 2.1.19
%pg.75 in Conduction Heat Transfer Solutions
s = 0.3048; %12in ->m
Bi_1 = (hair*r1)/k_PEH;
Bi_2 = (k_soil*burial_D*burial_D)/k_soil; % h=k*burialD
D = burial_D/s;
%SUMMER
t2 = T_Earths;
t1 = T_Gs;
q = 2*pi*k_soil*(t2-t1)/((1/Bi_1)+log((burial_D/(pi*r1*D))*sinh(2*pi*(D+(D/Bi_2)))));
T_delta3 = q *R_T*SA_tube*7; %7 or 8 is the number of tubes within each system
T_EarthPS = T_Earths-T_delta3 %Earth Temp Post-Summer
%WINTER
t2w = T_EarthPS;
t1w = T_Gw; %Was T_Gw, TE_w
q_w = 2*pi*k_soil*(t2w-t1w)/((1/Bi_1)+log((burial_D/(pi*r1*D))*sinh(2*pi*(D+(D/Bi_2)))));
%Semi-infinite solid. Surface condition: constant surface heat flux. pg337
%in Fundamental of Heat Transfer
%t = 60;
x = burial_D;
%T_EarthPW = t2w-((2*q_w*((diff_soil*t/pi)^1/2)/k_soil)*exp((-x^2)/(4*diff_soil*t)))-(((q_w*x)/k_soil)*erfc(x/(2*sqrt(diff_soil*t))))
t = 1:7.776e+6;
T_EarthPW = zeros(1,length(t));
q_w = zeros(1,length(t));
flag = 0;
for i=1:length(t)
    if i == 1 
        q_w(i) = 2*pi*k_soil*(t2w-t1w)/((1/Bi_1)+log((burial_D/(pi*r1*D))*sinh(2*pi*(D+(D/Bi_2)))));
        T_EarthPW(i) = t2w-((2*q_w(i)*((diff_soil*t(i)/pi)^1/2)/k_soil)*exp((-x^2)/(4*diff_soil*t(i))))-(((q_w(i)*x)/k_soil)*erfc(x/(2*sqrt(diff_soil*t(i)))));
    else
        T_EarthPW(i)=T_EarthPW(i-1)-((2*q_w(i-1)*((diff_soil*t(i)/pi)^1/2)/k_soil)*exp((-x^2)/(4*diff_soil*t(i))))-(((q_w(i)*x)/k_soil)*erfc(x/(2*sqrt(diff_soil*t(i)))));
        q_w(i) = 2*pi*k_soil*(T_EarthPW(i)-t1w)/((1/Bi_1)+log((burial_D/(pi*r1*D))*sinh(2*pi*(D+(D/Bi_2)))));
        if  T_EarthPW(i) < t1w*1.1 && T_EarthPW(i) > t1w && flag==0
            t_end = t(i)/(60*60);
             
            flag =1;
        end
    end
end

%% Functions
%Unit conversion
function Celsius = f_to_c(Fahrenheit)
    Celsius = (5/9)*(Fahrenheit-32);
end
function Fahrenheit = c_to_f(Celsius)
    Fahrenheit = (Celsius/(5/9))+32;
end

%summer_T = (sum([71.4	81.8	88.9	85.9	78.0])/5)
%winter_T = (sum([52.6	43.5	42.3	46.5	53.4])/5)
%https://wrcc.dri.edu/cgi-bin/cliMAIN.pl?utcirc 
%https://wcc.sc.egov.usda.gov/reportGenerator/view/...
%customSingleStationReport/hourly/start_of_period/2125:UT:SCAN%7Cid=%22%22%7Cname/2019-01-01,2019-01-31:M%7C1,D%7C1/stationId,TOBS::value,SRADV::value,WSPDV::value,STO:-2:value,STO:-4:value,STO:-8:value,STO:-20:value,STO:-40:value?fitToScreen=false
%METHOD 1 - FLUX into pipe -> Q into pipe -> Temp out factoring in thermal
%resistances  ** NOT WORKING, temps go out of possible range
%T_out2 = (erf(L_tube/(2*sqrt(diff_PEH*t)))*(T_EarthC- T_G)) + T_G;
%(T_out2-T_G)/(T_EarthC- T_G) = erf(L_tube/(2*sqrt(diff_PEH*t)))
%TRANSIENT SOLUTION SOLVE TO SEE T_EARTH GO BACK TO ~6C, t=?
%Q_CB3 = m_dot*Cp*((T_out3)-T_Gs)
%T_exp(i) = TE_s - q_w(i)*R_T*SA_tube*7;
%V = 33454.7; %ft^3
%V = V*0.0283168; %m^3
%Cp = 1006; %J/kg*K
%rho = 1.225; %0.6601; %@6000 kg/m^3
%V_s = 0.1012799;% m^3/s 
%m = V*rho; %kg
%m_dot = (V_s*rho); % kg/s
%diff_air = 20*10^-6; %m^2/s