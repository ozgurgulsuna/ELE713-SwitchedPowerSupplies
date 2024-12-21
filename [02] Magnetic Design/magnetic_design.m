close all
clear 
clc
%% Specify project requirements

Vimax = 60;
Vimin = 36;
Vo = 5;
Po = 40;
Vaux = 12;
Paux = 4;

% Choose an estimated efficiency
eff_est = 0.82;

% Choose inductor ripple ratio (ismail wants this)
ripple_ratio = 0.1625;

% Voltage ripple defined by specifications
delta_Vo = 0.1;  % Vpp


Io = Po/Vo;
Io_max = Io + Io*ripple_ratio/2;
Io_min = Io - Io*ripple_ratio/2;

% Auxiliary
Iaux = Paux/Vaux;


%% Duty Range vs Turns Ratio
% We observe the effect of turns ratio on the duty range of our converter
% given the input and output voltage range

% Calculated turns ratio for forward converter
% Vo/Vin = D*N2/N1

N1 = 1;

% We sweep N1/N2 between 1 and 5
N2 = linspace(1,5.9,50);

Dmax_vec = zeros(1,50);
Dmin_vec = zeros(1,50);

% Calculate min and max duty values for all N values
for i=1:50
    d1 = ((N2(i)) * ((Vo+Io*0.05)/Vimin));  % output inductance resistance
    d2 = ((N2(i)) * ((Vo+Io*0.05)/Vimax));
    dmax = max(d1,d2);
    dmin = min(d1,d2);
    Dmax_vec(i) = dmax;
    Dmin_vec(i) = dmin;
end

fig1 = figure;
plot(N2/N1,Dmax_vec);
hold on
plot(N2/N1,Dmin_vec);
% draw a line for the duty cycle limit
plot([0 5.9], [0.39 0.39], 'r--');
% draw a line for n = 2.6
plot([2.6 2.6], [0 0.5], 'r--');
legend
grid minor
title("Duty cycle limits vs N")
legend("D_{max}", "D_{min}")
xlabel("Turns ratio N"), ylabel("Duty cycle")
ylim([0 0.5])
hold off
% exportgraphics(fig1, "../../4-Report/img/DvsN.pdf")

% from the duty ranges we can select 
% Dmax = 0.39  and
% N1/N2 = 2.6

N2  = 1/2.6 ;
Dmax = (Vo+Io*0.05)/Vimin/N2;
%% Transformer Design
% Switch Voltage Limit: 200 V
% Vsw >= Vimax + Vimax*N1/N3 + Vsnub

Vsnub = 10 ;  % approximately
Vsw_max = 120;
N3_min = 1/((Vsw_max - Vimax - Vsnub)/(Vimax*N1));
% Vimax+Vinmax*N1/N3+Vsnub;
% N3_min = (Vsw - Vimax - Vsnub)/Vimax ;

% Dmax should be less than the Dmax level
Dmax_limit = 1/(1+N3_min/N1) ;


% Pick N3 as 1
N3 = 1;
Dmax_limit = 1/(1+N3/N1) ;


%% Determine Lm
clear dmin dmax Dmin_vec Dmax_vec d1 d2 i N
% 
Vin_vec = linspace(60, 36, 50);
D_vec = ((Vo+Io*0.05)*N1./(Vin_vec*N2));
fs = linspace(50e3, 200e3, 50);
Lmvec = zeros(1,50);
% delta_ILm_vec = ones(50, 50);  % Initial guess for delta_ILm, e.g., 0.1 A


%
Lm = 185e-6;
Fsw = 200e3;

% for i = 1:50
%     for j = 1:50
%         Lmvec(j) = Vimin*Dmax/(fs(i)*delta_ILm_vec(j))
%         delta_ILm_vec(i,j) = Vimin*Dmax/(fs(i)*Lmvec(j));
%     end
% end

% plot(fs,delta_ILm_vec)

delta_ILm = Vimin*Dmax/(Fsw*Lm);

% Choose max duty as 0.3900, calculate N
N2 = Vo/(Dmax*Vimin);

% Auxulary winding
N4 = Vaux/(Dmax*Vimin);

% Check switch peak current 
Io_pk_ref = Io_max*N2/N1/0.92;  % eff estimate correct ?
Isw_min = Io_min*N2/N1/0.92;
Isw_peak = delta_ILm+Io_pk_ref;   

% Choose switching frequency
%SOLVE DOUBLE ITERATION TO FIND DELTA IL AND LM 

% Determine the transformer current ripple


% Find required Lm versus Fsw
for i=1:50
    Lmvec(i) = (Vimin*Dmax)/(fs(i)*delta_ILm);
end

fig2= figure;
plot(fs, Lmvec/1e-6);
grid minor
xlabel("Switching Frequency (Hz)"), ylabel("L_m value (H)")
title("Lm for the ripple constraint vs. Fsw")

Ireset = (delta_ILm*(N3/N1)*Dmax)/2*1.83;

Vreset = Vimin;


%% Transformer Design POT-Core
% https://www.tdk-electronics.tdk.com/download/519704/069c210d0363d7b4682d9ff22c2ba503/ferrites-and-accessories-db-130501.pdf
Kw = 0.4;

% Select Core Material (N49)
performance_factor = 22500;


% Area Product Calculation: First Method
J = 5e6; % Current Density (A/m^2)
delta_B = performance_factor/Fsw; % Flux Density Swing (T)
P_trans = Vo*Io + Vaux*Iaux + (Vo*Io + Vaux*Iaux)/eff_est+ Vreset*Ireset; % Transmitted Power (P)

% First Ap is from the lecture notes
% Area Product Calculation: Method-1
display("Method-1")
Ap = ((11.1*Vo*Io/eff_est)/(0.141*delta_B*Fsw))^(1.31)*1e4

% Area Product Calculation: Method-2
display("Method-2")
Ap = ((Vo*Io/eff_est)/(0.014*delta_B*Fsw))^(4/3)*1e4

% Area Product Calculation: Method-3
display("Method-3")
C = 0.71; % Converter Coefficient (Forward Converter)
Ap = P_trans/(C*delta_B*Fsw*J)*1e12


% Area Product Calculation: Method-4
display("Method-4")
Kconv = 0.5;
Ap = Kconv*P_trans/(Kw*Fsw*delta_B*J)*1e12



% Selected RM type Core (RM 8 N49)
Ac = 64e-6;
Aw = 10.8*(17-8.55)/2e6;
Ap_core_rm8 = Ac*Aw*1e12
Al = 2200; % (nH/N^2)

% % Selected RM type Core (RM 7 N49 )
% Ac = 43e-6;
% Aw = 8.4*(14.75-7.25)/2e6;
% Ap_core_rm7 = Ac*Aw*1e12
% Al = 1900; % (nH/N^2)

% Ac = 31.3e-6;
% Aw = 8*(12.6-6.4)/2e6;
% Ap_core_rm6 = Ac*Aw*1e12
% Al = 2400; % (nH/N^2)

% % Selected ETD type Core (E 20/10/6 EF20)
% Ac = 32e-6;
% Aw = 14*(14.1-5.9)/2e6;
% Ap_core = Ac*Aw
 
% % Selected PQ type Core (PQ20/20)
% Ac = 62.1e-6;
% Aw = 14.3*(18-8.8)/2e6;
% Ap_core = Ac*Aw
% Al = 5200 % (nH/N^2)

% RM 10 N49
Ac = 98e-6;
Aw = 12.4*(21.2-10.9)/2e6;
Ap_core = Ac*Aw*1e12
Al = 2900; % (nH/N^2)

% 
% % RM 12
% Ac = 146e-6;
% Aw = 16.8*(25-12.8)/2e6;
% Ap_core = Ac*Aw
% Al = 5300; % (nH/N^2)

% % EFD 20/10/7
% Ac = 31e-6;
% Aw = 7.77*2*(15.4-8.9)/2;
% Ap_core = Ac*Aw*1e12
% Al = 910;


%% Turn Number and Saturation Calculations
% From primary side 
% Vp = N1 dΦ/dt

% Expected operation
Np = (Vimax)*D_vec(1)/(Fsw*delta_B*Ac);
Np = (Vimin)*Dmax/(Fsw*delta_B*Ac)



% At max instance, allow for max B. if turns ratio higher than this, core
% will not saturate.
B_sat = 0.3;

% Min turn numbers, less than that with saturate the core 
Np_min = (Vimax)*0.5/(Fsw*B_sat*Ac)

Np = 8
Nr = 8


Ns_min = Np*N2
Naux_min = Np*N4

Ns = 3
Na = 8

% Np = 8  PRIMARY
% Nr = 8  RESET
% Ns = 3  SECONDARY
% Na = 6  AUX


% Let's back-calculate the core flux density
% Lm = Np dΦ/i in microhenries
Lm_core = Al*Np*Np*1e-9*1e6

% in Tesla
delta_B_core = Lm*delta_ILm/(Np*Ac)

% in Tesla
delta_B_core = Lm_core*delta_ILm/(Np*Ac) 



%% Wire Gauge Calculation
% Determine Cable Length and Type
% https://masteringelectronicsdesign.com/the-rms-value-of-a-trapezoidal-waveform-part-2/
% https://masteringelectronicsdesign.com/how-to-derive-the-rms-value-of-a-triangle-waveform/
% Primary side assume square type voltage

Isw_mean = (Isw_peak+Isw_min)/2;
Isw_rms = Isw_mean*sqrt(Dmax);
Isw_rms = sqrt(Dmax/3*(Isw_peak^2+Isw_peak*Isw_min+Isw_min^2));

Isw_avg = Isw_mean*Dmax

Ipri_rms = Isw_rms;
Isec_rms = sqrt(Dmax/3*(Io_max^2+Io_max*Io_min+Io_min^2))

Ireset_rms = delta_ILm*(N3/N1)*sqrt(Dmax/3);


Iaux_rms = sqrt(Dmax/3*((Iaux*1.2)^2+Iaux^2+(Iaux*0.9)^2))

% Winding Area Distribution


A1 = Ipri_rms*Np
A2 = Isec_rms*Ns
A3 = Ireset_rms*Nr
A4 = Iaux_rms*Na

Atotal = A1 + A2 + A3 + A4

A1 = A1/Atotal
A2 = A2/Atotal
A3 = A3/Atotal
A4 = A4/Atotal



Pri_copper_cs = A1*Kw*Aw/Np*1e6
Sec_copper_cs = A2*Kw*Aw/Ns*1e6
Reset_copper_cs = A3*Kw*Aw/Nr*1e6
Aux_copper_cs = A4*Kw*Aw/Na*1e6



Jp=Ipri_rms/Pri_copper_cs
Js=Isec_rms/Sec_copper_cs
Jr=Ireset_rms/Reset_copper_cs
Ja=Iaux_rms/Aux_copper_cs














%% Litz
% Determine Cable Length and Type

% Set cable diameters:

% https://www.elektrisola.com/en/Litz-Wire/Info
litz_packing_factor = 1.28 ;


% TO DO: winding configuration sec/2 : pri : sec/2

%% Output Filter Calculation

Cf_cross_Lf = (Vo*(1-D_vec(1)))/(Fsw^2*8*delta_Vo) 

Lf = 16e-6 % SRP1265CC-220M reduces to 16 uH @ 10 A
Cf_min = Cf_cross_Lf/Lf

Cf = (5 + 5 + 5)*1e-6 ; % uF

delta_Vo_calculated = ((1-D_vec(1)))/(8*Lf*Cf*Fsw^2)*Vo;

delta_Il_max = ((Ns/Np*Vimax-Vo)*(Dmax)/(Fsw*Lf))
delta_Il = ((Ns/Np*Vimax-Vo)*(D_vec(1))/(Fsw*Lf))
Il_peak = Po/Vo+delta_Il/2

Caux_cross_Laux = (Vaux*(1-D_vec(1)))/(Fsw^2*8*Vaux*0.2);
Laux = 2.2e-6;
Caux_min = Caux_cross_Laux/Laux;
Caux = (5 + 5)*1e-6 ; % uF

delta_Iaux_max = ((Vaux)*(1-Dmax)/(Fsw*Laux))  %% simulation check




%% Power Loss

% Core Loss
% https://elnamagnetics.com/wp-content/uploads/library/Ferroxcube-Materials/3C94_Material_Specification.pdf

% Coefficients (with 95% confidence bounds):
p1 =   8.892e-05;
p2 =     0.01501;
p3 =      -0.775;
p4 =       17.63;

P_density_center = p1.*(B_field_center.*1000.*ripple_ratio).^3 + p2.*(B_field_center.*1000.*ripple_ratio).^2 + p3.*(B_field_center.*1000.*ripple_ratio) + p4; % mT -> kW per m^3
P_center = P_density_center.*(C*D*F*2)*1000 ;% core center volume

P_density_sides = p1.*(B_field_sides.*1000.*ripple_ratio).^3 + p2.*(B_field_sides.*1000.*ripple_ratio).^2 + p3.*(B_field_sides.*1000.*ripple_ratio) + p4; % mT -> kW per m^3
P_sides = P_density_sides.*(((A-B)*F*2*D)+((B-C)*(D-E)*F*2))*1000 ;  % core center volume

P_core = (P_center + P_sides)*2 ; % an error margin of 200 % is given

% Copper Loss for litz wire
% AC losses proximity and skin effect can be ignored
% Length of the wires is estimated from the middle radius

end_winding_pri = 0.265;
end_winding_sec = 0.080;

lenght_pri = primary_cable_count*(2*(C+0.005+primary_cable_diameter/2)+2*(F+0.005+primary_cable_diameter/2));
lenght_sec = secondary_cable_count*(((A+B)/2*2)+(5e-3+F)*2) ;

R_meas_pri = 18e-3;
R_pri = (cu_resistivity*lenght_pri)/(primary_parallel*pi*(primary_cable_diameter/2)^2*primary_parallel/litz_packing_factor);
R_pri = R_meas_pri*lenght_pri/(lenght_pri+end_winding_pri);

R_meas_sec = 102e-3;
R_sec = (cu_resistivity*lenght_sec)/(secondary_parallel*pi*(secondary_cable_diameter/2)^2*secondary_parallel/litz_packing_factor);
R_sec = (R_meas_sec*lenght_sec/secondary_parallel)/(lenght_sec/secondary_parallel+end_winding_sec)

P_cu_pri = R_pri*5.3*5.3;
P_cu_sec = R_sec*1.45*1.45;

P_cu_total = P_cu_pri + P_cu_sec ;

% https://downloads.hindawi.com/archive/2012/635715.pdf
L_leak = (mu0*N_pri^2*((0.00016)*2+((B-C)/2))*(F*E+D*(C+2*(((B-C)/2)))))/(3*2^2*E^2);