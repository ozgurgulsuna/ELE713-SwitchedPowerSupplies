close all
clear 
clc
%% Specify project requirements

Vimax = 60;
Vimin = 36;
Vo = 5;
Po = 40;
Vaux = 12;
Iaux = 0.3;

% Choose an estimated efficiency
eff_est = 0.82;
% Choose inductor ripple ratio (ismail wants this)
ripple_ratio = 0.10;

Io = Po/Vo;
Io_max = Io + Io*ripple_ratio/2;
Io_min = Io - Io*ripple_ratio/2;


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
legend
grid minor
title("Duty cycle limits vs N")
legend("D_{max}", "D_{min}")
xlabel("Turns ratio N"), ylabel("Duty cycle")
ylim([0 0.5])
hold off
% exportgraphics(fig1, "../../4-Report/img/DvsN.pdf")

% from the duty ranges we can select 
% Dmax = 0.3611  and
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
delta_ILm_vec = ones(50, 50);  % Initial guess for delta_ILm, e.g., 0.1 A


%
Lm = 185-6;
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


% Check switch peak current 
Io_pk_ref = Io_max*N2/N1/0.92;  % eff estimate correct ?
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

Ireset = (delta_ILm*(N3/N1)*D_vec(end))/2*1.83;
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

Ap = ((11.1*Vo*Io/eff_est)/(0.141*delta_B*Fsw))^(1.31)*1e4

Ap = ((Vo*Io/eff_est)/(0.014*delta_B*Fsw))^(4/3)*1e4

C = 0.71; % Converter Coefficient (Forward Converter)
Ap = P_trans/(C*delta_B*Fsw*J)*1e12

% Area Product Calculation: Second Method
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
Np = (Vimax)*D_vec(1)/(Fsw*delta_B*Ac)
Np = (Vimin)*D_vec(end)/(Fsw*delta_B*Ac)

% At max instance, allow for max B. if turns ratio higher than this, core
% will not saturate.
B_sat = 0.3;
Np_min = (Vimax)*0.5/(Fsw*B_sat*Ac)

% The results is yielded as 10.41, Pick N1 = 10;

Ns = Np*N2

Np = 8
Ns = 3

% Let's back-calculate the core flux density
% Lm = Np dΦ/i

Lm_core = Al*Np*Np*1e-9

delta_B_core = Lm*delta_ILm/(Np*Ac)
delta_B_core = Lm_core*delta_ILm/(Np*Ac) 


Lm_core = Al*Np^2*1e-9*1e6

Np = Lm_core*1e-6*delta_ILm/(Ac*delta_B_core)



Vo/(Fsw*delta_B*Ac)



%% and Wire Gauge Calculation




%% Choose Gap Length and Turn Number
% https://www.ferroxcube.com/upload/media/product/file/Pr_ds/E30_15_7.pdf



%% AWG
% Determine Cable Length and Type




%% Litz
% Determine Cable Length and Type

% Set cable diameters:

% https://www.elektrisola.com/en/Litz-Wire/Info
litz_packing_factor = 1.28 ;


% TO DO: winding configuration sec/2 : pri : sec/2

%% Determine Max B


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