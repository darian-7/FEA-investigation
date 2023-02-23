clear all
close all
%% Question 2a

% problem variables

sig_max = 40; % MPa
sig_min = 0;
t = 10e-3; % thickness
D = 200e-3; % outer diameter
R = 90e-3; % inner diameter
K_IC = 90; % MPa / m^0.5
C = 1e-12; % m/cycle per ~(MPa / m^0.5)^m
m = 4;
af = 10e-3; % crack length at failure due to leakage --> af = thickness
ao = 0.1e-3; % arbitrary starting crack length

sig_axial = (sig_max*((2*R)/2))/(2*t);

N = [];
n1 = 1;

% Crack Growth Rate Iterations

while ao <= 0.01

    N(n1,1) = ao*1000; % initial half crack length in mm

    fun = @(ao) 1./((C).*((0.728 + (0.373*(ao/0.01)^2) - (0.029*(ao/0.01)^4)).*sig_axial.*sqrt(pi.*ao)).^4);

    N(n1,2) = integral(fun,ao,af,"ArrayValued",true);

    ao = ao + 1e-4; % increase crack length by 0.1mm
    
    n1 = n1 + 1; 

end

figure
plot(N(:,2),N(:,1))
ylabel('Half Crack Length [mm]')
xlabel('Number of Cycles until Failure')

% Number of cycles from initial to failure

paris = [];
n2 = 1;
a = 1e-3; % initial crack length
a1 = 1;

while a <= af

    Y = 0.728 + (0.373*(a/0.01)^2) - (0.029*(a/0.01)^4);

    K = Y*sig_axial*sqrt(pi*a);

    paris(n2,1) = a*1000; % crack radius

    paris(n2,2) = 1 - fn_pod(2018717,a*1000); % prob of not detection cracks above this radius

    paris(n2,3) = fn_pod(2018717,a*1000); % prob of detection
    
    paris(n2,4) = n2;

    dadN(n2,1) = a*1000;

    dadN(n2,2) = C*(K^m); % cack growth rate at that crack radius

    a = a + dadN(n2,2);

    n2 = n2 + 1;


end

figure
plot(paris(:,1),paris(:,3))
xlabel('Crack Radius "a" [mm]')
ylabel('Probability of Detection Beyond "a"')

figure
plot(paris(:,1),paris(:,2))
xlabel('Crack Length "a" [mm]')
ylabel('Probability of Undetection Beyond "a"')

%% Failure Variables

NC = length(dadN);
xaxis1(:,1) = 1:NC;

failure = [];

cycles = 13000:1:16000; % cycles per inspection
cycles = transpose(cycles);

insp_int = [];
a2 = 1;
failure = [];
IIT = {};

%% Finding Probability of undetection at inspection number

for j = 1:length(cycles)

    a2 = 1;

    cycle = cycles(j);

    for i = 1:length(paris)
        
        insp_int(i,1) = xaxis1(i,1)/cycle;
    
        if rem(insp_int(i,1),1) == 0
    
            IIT{j,1}(a2,1) = insp_int(i,1); % inspection number
            IIT{j,1}(a2,2) = paris(i,2); % prob. of undetection
            IIT{j,1}(a2,3) = paris(i,1); % crack radius @ inspection

            a2 = a2 + 1;
    
        end

    end

end

%% Finding probability of failure for each cycle variable

for i = 1:length(IIT)

    for k = 2:length(IIT{i,1})

        IIT{i,1}(1,4) = IIT{i,1}(1,2);
    
        IIT{i,1}(k,4) = IIT{i,1}(k,2) * IIT{i,1}((k-1),4);

    end

    failure(i,1) = cycles(i); % cycles per iteration
    failure(i,2) = IIT{i,1}(end,4); % probability of failure

end

% Graph to show all variations of inspection
figure
scatter(failure(:,1),failure(:,2), 7)
hold on
yline(0.01)
xlabel('No. of Cycles per Inspection')
ylabel('Probability of Failure')
hold off

% Graph based on inspection for 8000:12000 cycles
figure
scatter(IIT{316,1}(:,1),IIT{316,1}(:,4), 7)
hold on
scatter(IIT{317,1}(:,1),IIT{317,1}(:,4), 7)
yline(0.01)
xlabel('Inspection Number')
ylabel('Probability of Failure @ Inspection')
hold off

% Graph based on inspection for 13000:16000 cycles
figure
scatter(IIT{361,1}(:,1),IIT{361,1}(:,4), 7)
hold on
scatter(IIT{362,1}(:,1),IIT{362,1}(:,4), 7)
yline(0.01)
xlabel('Inspection Number')
ylabel('Probability of Failure @ Inspection')
hold off

% Graph of Probability of non-detection vs number of cycles
figure
plot(paris(:,4), paris(:,2))
xlabel('Cycle Number')
ylabel('Probability of Non-detection')