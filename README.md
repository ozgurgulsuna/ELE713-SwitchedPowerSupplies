# INTRODUCTION

This report presents the design and implementation of a forward converter that operates with an input voltage range of 36-60 V DC, providing a primary output of 5 V at 8 A along with a secondary auxiliary output of 12 V at 4 W. The forward converter is a type of DC-to-DC converter that transfers power during the forward operation of a transformer, making it suitable for applications requiring isolated and regulated outputs.

To achieve precise control over the output voltage, the design employs an analog controller integrated with a feedforward loop. The primary side of the converter utilizes a single switch, while the transformer flux is reset using an additional winding and a diode. This approach ensures efficient operation and prevents magnetic saturation of the transformer core.

In low-voltage, high-current applications, traditional diode-based rectification becomes less efficient due to the voltage drop across the diodes, typically around 0.5 V. To address this inefficiency, synchronous rectification is implemented using MOSFETs at the output. When conducting, a MOSFET with a 10-milliohm resistance produces a voltage drop of only 0.1 V at 10 A, significantly improving overall efficiency.

The printed circuit board (PCB) for this design is manufactured using a photolithography process. This process includes applying a photoresistive UV coating and etching the copper layer using a solution of hydrogen peroxide and hydrochloric acid. The circuit's functionality is thoroughly tested under various load conditions, including overload scenarios and auxiliary load operation, to validate its performance and reliability.

This report further discusses the chosen topology, magnetic design, and key considerations in achieving an efficient, compact, and robust forward converter design.

| Specification        |          |
|--------------|-------|
| Input Voltage                   | 36-60 Vdc (Nominal 48Vdc)            |
| Main Output Voltage             | 5V (min. 4.9V, max. 5.1V)            |
| Output Voltage Ripple           | 100mV peak-to-peak                   |
| Main Output Nominal Current     | 8A                                   |
| Continuous Output Power         | 40W (Main Output)                    |
| Overload Current                | 9A                                   |
| Auxiliary Output Voltage        | 12V                                  |
| Auxiliary Output Nominal Current| 0.3A                                 |
| Auxiliary Output Power          | 4W                                   |
| Line Regulation                 | $\pm 0.2$%                           |
| Load Regulation                 | $\pm 0.5$%                           |
| Transient Peak Deviation        | 10% of Vout (10% load to 100% @ 48 V) |
| Transient Response Recovery     | 200 $\mu sec$ (10% load to 100% @ 48 V) |
| Overall Efficiency @ Full Load  |$\ge$ 80%                             |


---

# TOPOLOGY
The forward converter is a widely used topology in DC-to-DC conversion, especially for applications requiring electrical isolation and multiple outputs. This topology operates by transferring energy through a transformer during the on-state of the primary-side switch, with energy delivery controlled by the duty cycle of the switch. It is a buck-type converter implemented with the output LC filter structure.

![Forward Converter](https://github.com/user-attachments/assets/e35db0a6-66df-4666-a457-e994187df0b1)
*Schematic of Forward Converter*

The duty ratio, defined as the proportion of time the primary switch is on during a switching cycle, is a critical parameter in determining the output voltage of the forward converter. The relationship between the input voltage, duty ratio, and transformer turns ratio governs the output voltage. A feedback mechanism ensures the duty ratio is adjusted dynamically to maintain a stable output under varying input and load conditions. In a forward converter, the maximum duty cycle is required to be less than 50% to ensure the transformer reset.

$${V_o \over V_i} = {Ns\over Np}*D$$

At the first cycle when the transistor T is ON, the transformer transfers power to the secondary side, and reset winding is reverse biased. In the second cycle, reset winding resets the transformer and supplies the magnetizing current back to the supply, and the load is supplied by the inductor and capacitor. Reset winding is very important to reset to transformer to prevent core saturation. In our design, we need to add two more windings to the transformer, one is to supply the primary side control circuitry, and another is to create the required 12V auxiliary supply.

The forward converter design includes a bias supply for the primary side, which provides the necessary power to the control circuitry and auxiliary components. This supply is derived from the bias winding connected to the primary side. The circuit initially starts with a voltage divider and a zener diode. As the bias winding experiences a change in flux, it supplies the rest of the circuit with the required voltage. The output of the bias winding is configured as a rectifier of a forward converter and is regulated with a linear regulator constructed from an NPN transistor. This sets the bias supply to approximately 12 volts.

![bias](https://github.com/user-attachments/assets/b54f5da7-5755-426d-8909-4d48eae5e08a)
*Schematic of the Bias Voltage Generation Subcircuit.*


The secondary side of the converter includes an auxiliary output in addition to the primary output. This output is set to 12 volts and is rated for 0.3 amperes, resulting in a 4 W auxiliary output. The auxiliary output is supplied from another auxiliary winding with a forward rectifier output. This output is also regulated using a 12-volt linear regulator and is referenced to the ground of the secondary winding.

![aux](https://github.com/user-attachments/assets/d044b516-b183-4307-97ac-a812cf31d17e)
*Schematic of the Auxiliary Voltage Generation Subcircuit.*

---

# MAGNETIC DESIGN

Magnetic design is a critical aspect of the forward converter, as it directly impacts efficiency, size, and reliability. The transformer in the forward converter serves two primary purposes: providing electrical isolation and stepping the voltage up or down based on the turn ratio.

A forward converter utilizes a transformer operating in forward mode, meaning that energy is transferred directly through the transformer core during the on-state of the primary switch. Unlike some other topologies, energy is not stored in the core during one cycle and transferred in the next; instead, the transfer is immediate.

The magnetic design is done in a MATLAB script that calculates key parameters for a forward converter, focusing on component sizing and operational constraints. It starts by defining input/output voltage ranges, power requirements, efficiency assumptions, and ripple constraints. The output current and its ripple are calculated to establish operating boundaries.

A significant part of the script focuses on determining the duty cycle range as a function of the transformer turns ratio $$(\N_1/N_2\)$$ and selecting an optimal ratio to ensure the converter operates within specified limits. Transformer design calculations include determining winding turns for primary, secondary, and auxiliary windings while ensuring compliance with voltage and current constraints, such as switch voltage limits. The calculated turn ratios are given in the table below.

| Winding      | Turns |
|--------------|-------|
| Primary      | 10    |
| Reset        | 10    |
| Secondary    | 4     |
| Bias         | 7     |
| Auxiliary    | 7     |

*Selected Turns Ratios.*

The script calculates the magnetizing inductance $$\(L_m\)$$ and evaluates the impact of switching frequency on inductor performance. It also designs the output filter components to meet ripple voltage constraints.

Core selection is supported by area product $$\(A_p\)$$ calculations using various methods, ensuring the chosen core can handle the required flux and thermal limits. Wire gauge and winding distribution are calculated to balance current density and copper usage.

Finally, the script includes a loss analysis, estimating core and copper losses, and evaluating leakage inductance to ensure the design meets efficiency and performance goals. This approach focuses on key aspects of forward converter design while avoiding unnecessary complexity.

The code can be found in [magnetic design folder](%5B02%5D%20Magnetic%20Design/magnetic_design.m).

## TRANSFORMER REALIZATION ##

After evaluating the available cores, calculations are made, and a few candidate cores are chosen based on area product criteria. As a final decision, E30/15/7 3C94 is chosen (by Ferroxcube). There were smaller options; however, this core is favored considering the ease of implementation. Our switching frequency is chosen as 200kHz. At 200kHz, which is moderately high, we used $1 mm^2$ litz wires to avoid the skin effect and proximity effect. For the primary we have 10 turns with no parallel wires, at the secondary, we have 4 turns with 2 parallel wires. The auxiliary winding is winded using 0.2 $mm^2$ litz cable to have lower losses. Reset winding and primary bias circuitry are winded with thin copper cables so that they won't carry high currents. The order of windings is secondary-auxiliary-reset-bias-primary.


<div style="text-align: center;">
  <img src="https://github.com/user-attachments/assets/303ca2ac-c7ce-4901-8131-cbb1feaca27f" alt="PCB Preprocess" width="500" height="750" />
  <p><em>Transformer Windings Drawing</em></p>
</div>


---
# CONTROLLER DESIGN
For the controller, a specific analog integrated circuit, LT1952-1 produced by Analog Devices, is used to implement the controller. LT 1952-1 is a single-switch synchronous forward controller for forward controller designs within the range of 25-500 W. This controller is favored because it has the option to implement synchronous rectification. and high efficiencies can be achieved. LT1952 is the main controller; however, there are 2 more integrated circuits for the essential working of the converter. TL431 is used as a feedback compensator for controlling the output voltage. Moreover, LTC3900 is used for synchronous rectification which takes the SYNC signal from the LT1952-1 with the help of pulse transformer and switches two MOSFETs for improving the efficiency (with the risk of short-circuiting the secondary side MOSFETs).

---
# PRODUCTION AND ASSEMBLED PCB
The printed circuit board (PCB) for this design is manufactured using a photolithography process. This process includes applying a photoresistive UV coating and etching the copper layer using a solution of hydrogen peroxide and hydrochloric acid. The circuit's functionality is thoroughly tested under various load conditions, including overload scenarios and auxiliary load operation, to validate its performance and reliability.

<div style="text-align: center;">
  <img src="https://github.com/user-attachments/assets/dce49cd1-3a43-4067-a34e-a12fab11281c" alt="PCB Preprocess" width="500" />
  <p><em>Final Results of the etched PCB.</em></p>
</div>

<div style="text-align: center;">
  <img src="https://github.com/user-attachments/assets/f288cc9c-4eab-49d1-bcc6-069adfa5558a" alt="PCB Preprocess" width="450" />
  <img src="https://github.com/user-attachments/assets/e04cd08b-fa40-4647-a380-e396ca59438b" alt="PCB Preprocess" width="500" />
  <p><em>Assembled PCB</em></p>
	        
</div>

## RESULTS ##
After the design is completed and assembled some measurements are taken. At first, some problems occurred at the current sense pin of the controller, which was due to a grounding problem. We solved it by fixing it. 
![cs_fault](https://github.com/user-attachments/assets/d1a69e69-13e3-4eed-968b-fc26dccb8cf2)
*Converter Output in the case of fault.*

After fixing that, tests were continued without synchronous switches which provided an efficiency of around 75%. We activated the synchronous rectification to increase efficiency and obtain better results. At this part, we had some problems with the synchronization and failed the synchronize both switches. However, synchronizing one switch had no disadvantage but a slight improvement so we obtained a  semi-synchronous converter. Forward converter without auxiliary is quite efficient; however, the auxiliary linear regulator creates quite a loss. The results are as follows:

![5v_rated](https://github.com/user-attachments/assets/05da8ae4-b400-4627-9c6d-f33e8ce721b7)
*5V Rated Operation*
![Ripple_v1](https://github.com/user-attachments/assets/afe493c5-0d2a-45a8-b23a-3706a55110a7)
*Voltage Ripple Measurement*

Efficiency Calculations(@48V):

Only Forward Converter

${Efficiency} = {Pout \over Pin}$ =  $`{(5.15 × 7.21) \over (48.5 × 1.04)}`$ = 0.8602

Forward Converter with 4W auxiliary load

${Efficiency} = {Pout \over Pin}$ =  $`{(5.15 × 7.21 + 4.3) \over (48.5 × 1.04)}`$ = 0.8214

I would like to point out that @60V efficiency gets worse due to auxiliary efficiency. 

---

# CONCLUSION

In conclusion, there are several isolated DC-DC converter topologies with all their advantages
and disadvantages. In this project, forward converter topology is implemented to convert 36-
60V to 5V with the addition of 12V auxiliary voltage. Moreover, the project covers magnetic design, analog design, and feedback design
for isolated DC/DC converters. The design showed that the topology proved to be working;
however, there is still room for improvement in terms of PCB design, feedback design, and
magnetic design. 

---
# APPENDIX

```matlab
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
Lm = 200e-6;
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
performance_factor = 20000;


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

% % RM 10 N49
% Ac = 98e-6;
% Aw = 12.4*(21.2-10.9)/2e6;
% Ap_core = Ac*Aw*1e12
% Al = 2900; % (nH/N^2)

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

% % ETD 34/17/11
% Ac = 97.1e-6;
% Aw = 11.8*2*(25.6-11.1)/2e6;
% Ap_core_ETD = Ac*Aw*1e12
% Al = 2600;

% % P18/11 SER
% Ac = 43.3e-6;
% Aw = (14.9-7.60)/2*7.2/1e6;
% Ap_core_P18 = Ac*Aw*1e12
% Al = 260;

% % EF 25 6 11 N97
% Ac = 67.2e-6;
% Aw = (21.7-9.4)/2*6.2/1e6;
% Ap_core_EF25 = Ac*Aw*1e12
% Al = 4100;

% E30/15/7 3C94
Ac = 60e-6;
Aw = (19.5-7.2)/2*9.7*2/1e6;
Ap_core_E30 = Ac*Aw*1e12
Al = 1900;

% % E20/10/6
% Ac = 32e-6;
% Aw = (14.1-5.9)/2*7*2/1e6;
% Ap_core_E20 = Ac*Aw*1e12
% Al = 1350;

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

Np = 10
Nr = 10


Ns_min = Np*N2
Naux_min = Np*N4

Ns = 4
Na = 10

% Np = 8  PRIMARY
% Nr = 8  RESET
% Ns = 3  SECONDARY
% Na = 6  AUX


% Let's back-calculate the core flux density
% Lm = Np dΦ/i in microhenries
Lm_core = Al*Np*Np*1e-9

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
```

