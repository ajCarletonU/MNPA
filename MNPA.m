%ELEC 4700
%PA7 MNA Building
%Alina Jacobson  (4 March 2022)

clc; 
close all; 
clear all;

set(0, 'DefaultFigureWindowStyle', 'docked')
global G C b; %define global variables


%equations KCL
%node1: 
        %0 = IE - (V2-V1)/R1 - d/dt (V1-V2)C
        %IE = (V2-V1)/R1 - jwC (V2-V1)
        %and V1 = Vin

%node2:
        %0=(V2-V1)/R1 - d/dt (V2-V1)C + V2/R2 + IL
        %0=(V2-V1)/R1 - jwC (V2-V1) + V2/R2 + 1/jwL (V2-V3)

%node3:
        %0=V3/R3 -IL
        %0=V3/RS - 1/jwL (V2-V3)
        
%node4:
        % V4= a*I3 
        % V4= a* V3/R3
        
%node5:
        %0=(V4-V5)/R4 - V5/RO / v5
        %0 = V4(Ro / R4+Ro


%Components
%withR1 =1,C=0.25,R2 =2,L=0.2,R3 =10,α=100,R4 =0.1andRO =1000.
r1=1;
r2=2;
r3=10;
r4=0.1;
ro=1000;

c1=0.25;
l=0.2;
a=100;


%conductance
g1 =1/r1;       %node 1 to 2
g2 =1/r2;       %node 2 to ground
g3 =1/r3;       %node 3 to ground
g4 =1/r4;       %node 4 to 5
go =1/ro;       %node 5 to ground

g22=g1+g2;      % at node 2
g55=g4+go;      % at node 5

aa=-a/r3;

%current
IE = 1;
IL=1;
Iv =1;

% # of nodes = 5
% # of inudance = 1
% # of independent voltage sources = 1
% # of voltage control voltage source VCVS = 1
% MNA size = 8

%the G matrix:
%   [v1  v2  v3  v4  v5  IE IL Iv]
G = [g1 -g1  0   0   0   IE 0  0         %v1
    -g1  g22 0   0   0   0  IL 0         %v2
     0   0   g3  0   0   0 -IL 0         %v3
     0   0   0   g4 -g4  0  0  Iv        %v4
     0   0   0  -g4  g55 0  0  0         %v5
     IE  0   0   0   0   0  0  0         %ie
     0   IL -IL  0   0   0  0  0         %il
     0   0   aa  Iv  0   0  0  0];       %iv

%the C matrix:  
%   [v1  v2 v3 v4 v5 IE IL IV]
C = [c1 -c1 0  0  0  0  0  0        %v1
    -c1  c1 0  0  0  0  0  0        %v2
     0   0  0  0  0  0  0  0        %v3
     0   0  0  0  0  0  0  0        %v4
     0   0  0  0  0  0  0  0        %v5
     0   0  0  0  0  0  0  0        %ie
     0   0  0  0  0  0 -l  0        %il
     0   0  0  0  0  0  0  0];      %iv
 
 
%source column vector b:
b = [0;  0; 0; 0; 0; 1; 0; 0];  

%--------------------------------------------------------------------------
%(b) For the DC case sweep the input voltage V1 from -10V to 10V and plot VO and the voltage at V3.
%--------------------------------------------------------------------------
vmin = -10;
vmax = 10;
num=1000;

v1_in = linspace(vmin,vmax,num);    %input voltage V1 from -10V to 10V
V_3 = zeros(size(v1_in));
V_O = zeros(size(v1_in));
s=0;                                %at 0 for DC 

% %Loop thru until system reaches end of size
for j = 1:size(v1_in,2)
    A=(G+s.*C);
    B=(b.*v1_in(j));
    X = A\B;                      % calc  AX=b  =>  x = A\b    A=G+sC
    V_3(j) = X(3);                % v3 at node 3    
    V_O(j) = X(5);                % VO at node 5
end


%plot output DC sweep voltages
figure('Name',"Part B - DC Sweep")
plot(v1_in,V_3,'k','LineWidth',3);
grid;
hold on
plot(v1_in,V_O,'r','LineWidth',3);
title('VO and the voltage at V3 (RANGE: -10 V to 10 V)');
xlabel('Input Voltage (Vin)');
ylabel('Voltage (V)');
legend('Voltage at V3', 'Output voltage (VO)');


%--------------------------------------------------------------------------
%(c) For the AC case plot VO as a function of ω also plot the gain VO/V1 in dB.
%--------------------------------------------------------------------------
fmin= 1;      %hz
fmax=100;     %hz
itr=10000;        

F = linspace(fmin, fmax, itr);
Vout = zeros(size(F));
for j=1:itr
    w = 2*pi*F(j); %solve omega
    s = 1i*w;       %solve s
    A = G+s.*C;      
    X = A\b;        % solve AX=b ->  x = A\b
    
    Vout(j) = X(5);  % The voltage at output node = 5
                     % is collected in an array "Vout(n)" for every frequency.
    
end

%plot output AC sweep voltages
figure('Name',"Part C - AC Sweep")
plot(F, abs(Vout),'LineWidth',3);
grid;
title('AC sweep');
xlabel('w (rad/s)');
ylabel('|V_{out}|  (Volts)');
legend('Vout at Node 5');

figure('Name',"Part C - AC Sweep Gain ")
vOutLog = 20*log(abs(Vout));            % convert vout value to log
semilogx(F, vOutLog,'LineWidth',3);       %plot(x,y) but make x axis log scale
grid;
title('AC sweep in Vo/V1 dB vs Frequency (hz) log scale');
xlabel('frequency (hz)');
ylabel('Vout/V1  (dB)');


%--------------------------------------------------------------------------
%(d) For the AC case plot the gain as function of random perturbations on C 
% using a normal distribution with std = .05 at ω = π. Do a histogram of the gain.
%--------------------------------------------------------------------------

std=0.05;
c_norm_dist = std.*randn(1000,1) +c1;
V_out = zeros(size(c_norm_dist));

for j = 1:size(c_norm_dist)
    
    %4 elements of c in matrix
    C(1,1) =  c_norm_dist(j);
    C(1,2) = -c_norm_dist(j);
    C(2,1) = -c_norm_dist(j);
    C(2,2) =  c_norm_dist(j);
    
    
    w = pi;
    s = 1i*w;       %solve s
    A = G+s.*C;      
    X = A\b;        % solve AX=b ->  x = A\b
    V_out(j) = X(5);  % The voltage at output node = 5
    
end

figure('Name',"Part D - histogram of C")
histogram(c_norm_dist);
title('C Dist ');
xlabel('C (F)');
ylabel('number');

figure('Name',"Part D - histogram of the gain")
histogram(abs(V_out));
title('V_out ');
xlabel(' Vout','FontSize',12);
ylabel('number');