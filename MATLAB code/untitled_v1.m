clear all
close all
%% Question 1b

E = 70*10^9;
v = 0.3;

ex = fn_strain(1944637,0);
    


ey = fn_strain(1944637,90);

Ek = E/(1-v^2);

sxx = Ek*(ex+v*ey);
syy = Ek*(ey+v*ex);

%% Question 1c

spp = -400000000:1000000:400000000; % range of stress values to test for applies stress

% k1 = (sxx+syy+spp(i))/2;

% k2 = ((sxx/2)-(syy/2)+spp(i)/2)^2;

% k3 = (3/16)*(syy-sxx)^2;

sig1 = zeros(size(spp)); % principal stress 1 for all stress values

sig2 = zeros(size(spp)); % principal stress 2 for all stress values

for i = 1:801
    
    sig1(i) = (sxx+syy+spp(i))/2 + sqrt((((sxx/2)-(syy/2)+spp(i))/2)^2 + ((3/16)*(syy-sxx)^2));

    sig2(i) = (sxx+syy+spp(i))/2 - sqrt((((sxx/2)-(syy/2)+spp(i))/2)^2 + ((3/16)*(syy-sxx)^2));

end

sig1 = transpose(sig1);
sig2 = transpose(sig2);

% for loop to find tresca values

tresca1 = zeros(size(spp));
tresca2 = zeros(size(spp));
tresca3 = zeros(size(spp));
trescamax = zeros(size(spp));

for t1 = 1:801
    
    tresca1(t1) = abs(sig1(t1)-sig2(t1));
    tresca2(t1) = abs(sig2(t1));
    tresca3(t1) = abs(-sig1(t1));
    tresca_all = [tresca1(t1) tresca2(t1) tresca3(t1)];
    trescamax(t1) = max(tresca_all);

end

trescamax = transpose(trescamax);

% for loop to find von mises values

mises = zeros(size(spp));

r2 = 1/sqrt(2);

for t2 = 1:801

    mises(t2) = r2*sqrt(((sig1(t2)-sig2(t2))^2)+((sig2(t2))^2)+((-sig1(t2))^2));

end

mises = transpose(mises);

% converging von mises and tresca

convergetresca = [];
convergemises = [];

for c1 = 1:801

    if trescamax(c1) <= 400*10^6

        convergetresca(c1,1) = spp(c1);
        convergetresca(c1,2) = trescamax(c1);
        
    else
        
        continue

    end

end

for c2 = 1:801

    if mises(c2) <= 400*10^6

        convergemises(c2,1) = spp(c2);
        convergemises(c2,2) = mises(c2);
        
    else
        
        continue

    end

end

% safety factors for both von mises and tresca

SFV = [];
SFT = [];
sy = 400*10^6;

for sft = 1:length(convergetresca)
    
%     if convergetresca(sft,2) ~= 0

        SFT(sft) = sy/convergetresca(sft,2);
    
%     else
%         
%         continue
%         
%     end

end

for sfv = 1:length(convergemises)
    
%     if convergemises(sfv,2) ~= 0

        SFV(sfv) = sy/convergemises(sfv,2);
        
%     else
%         
%         continue
% 
%     end

end

SFV = transpose(SFV);
SFT = transpose(SFT);

%% Crack Orientation Problem

a = 0.0025;

K = zeros(size(289161));

design_stress = zeros(size(289161));

orientation = [];

sigma_p = [];

j = 1;
b = 1;

C = {};

Sxx = (1/4)*(3*sxx + syy);
Syy = (1/4)*(sxx + 3*syy);
Sxy = (sqrt(3)/4)*(syy-sxx);

%% LEFM & EPFM

for p = -400*10^6:10^6:400*10^6

    for theta = 0:(pi/180):pi

        final_sig = [ Sxx + p, Sxy; Sxy, Syy ];

        A = [ cos(theta), sin(theta); -sin(theta), cos(theta)];

        AT = transpose(A);

        t_final_sig = A*final_sig*AT;

        design_stress(j) = t_final_sig(2,2) - t_final_sig(1,2); % subtracted parallel stress from perpendicular

        orientation(j) = theta;

        sigma_p(j) = p;

        K(j) = (t_final_sig(2,2) - t_final_sig(1,2))*sqrt(pi*a);

        Keff(j) = K(j).*sqrt(1-0.5*((t_final_sig(2,2) - t_final_sig(1,2))/sy));

        j = j+1;

    end

end


% full_orientation=repmat(orientation,26,1);
% 
% A = [];
% 
% u = 1;
% 
% for m = 1:length(design_stress)
%     
%     for a = 0.0025:0.0001:0.005
% 
%         A(u) = a;
%     
%         K(u) = design_stress(m)*sqrt(pi*a);
% 
%         u = u + 1;
%     
%     end
% 
% end


design_stress = transpose(design_stress);
orientation = transpose(orientation);
sigma_p = transpose(sigma_p);

K = transpose(K);
Keff= transpose(Keff);

FINALMATRIX = [sigma_p orientation K];
FINALMATRIX_EFF = [sigma_p orientation Keff];

% Creating cell array

row = 181;
grp = size(FINALMATRIX,1)/row;

C = mat2cell(FINALMATRIX,row*ones(1,grp),size(FINALMATRIX,2));

CC = mat2cell(FINALMATRIX_EFF,row*ones(1,grp),size(FINALMATRIX_EFF,2));

for w = 1:grp
    C{w};
end

for w = 1:grp
    CC{w};
end

% finding max fracture toughness & orientaion for each applied stress

for q = 1:length(C)

    [M, I] = max(C{q,1}(:,3));
    C1{q,3} = M;
    C1{q,2} = C{q,1}(I,2);
    C1{q,1} = C{q,1}(I,1);

end

for q = 1:length(CC)

    [MM, II] = max(CC{q,1}(:,3));
    CC1{q,3} = MM;
    CC1{q,2} = CC{q,1}(II,2);
    CC1{q,1} = CC{q,1}(II,1);

end

mesh_k = [];

for a1 = 1:801

    mesh_k1 = [C{a1,1}(:,3)];

    mesh_k = horzcat(mesh_k, mesh_k1);    

end

% mesh_angle = 0:(pi/180):2*pi;
mesh_angle = 0:(pi/180):pi;
mesh_stress = -400e6:1e6:400e6;
[X,Y] = meshgrid(mesh_stress, mesh_angle);

fracture_polar = [];
f1 = 1;
f2 = 1;
j=1;

for i = 1:181
    
    for j = 1:801

        if mesh_k(i,j) > 2e7
    
            truth(i,j) = 1;  % number of fractures for that value of theta

        else

            truth(i,j) = 0;
    
        end

    end

end


for f2 = 1:181

    fracture_polar(f2,1) =  mesh_angle(1,f2);
    
    fracture_polar(f2,2) = sum(truth(f2,:));
    
end

figure
plot(fracture_polar(:,1),fracture_polar(:,2))
xlabel('Crack Orientation')
ylabel('Number of Failures from Varying Applied Stresses')

figure
cont = contourf(X,Y,mesh_k);
hold on
co = colorbar;
ylabel('Crack Orientation [radians]')
xlabel('Applied Stress [MPa]')
co.Label.String = 'Stress Intensity Factor - K [MPam^0.5]';
hold off

C2_LEFM = cell2mat(C1); % all applied stress values matched with K value and crack angle
C2_EPFM = cell2mat(CC1); % for epfm values
C2L1 = [];
C2L2 = [];
p2 = 1;
p3 = 1;

for p1 = 1:801

    if C2_LEFM(p1,2) < 3

        C2L1(p2,1) = C2_LEFM(p1,1);
        C2L1(p2,2) = C2_LEFM(p1,2);
        C2L1(p2,3) = C2_LEFM(p1,3);

        p2 = p2 + 1;

    elseif C2_LEFM(p1,2) > 3

        C2L2(p3,1) = C2_LEFM(p1,1);
        C2L2(p3,2) = C2_LEFM(p1,2);
        C2L2(p3,3) = C2_LEFM(p1,3);

        p3 = p3 + 1;

    end
    
end

% graph to show
figure
plot(C2L1(:,1),C2L1(:,2),'k')
xlabel('Applied Stress [MPa]')
ylabel('Crack Orientation [radians]')

% graph to show most dangerous crack orientation
figure
polarplot(pi/6 + fracture_polar(:,1),mesh_k1)
pax = gca;
pax.ThetaAxisUnits = 'degrees';
thetaticks(0:45:360)
thetatickformat('degrees')

% graph for different k values for max tresca stress
figure
polarplot(pi/6 + C{783,1}(:,2),C{783,1}(:,3))
pax = gca;
pax.ThetaAxisUnits = 'degrees';
thetaticks(0:45:360)
thetatickformat('degrees')

% graph for all max K values for all applied stress
figure
plot(C2_LEFM(:,1),C2_LEFM(:,3))
xlabel('Applied Stress [MPa]')
ylabel('Fracture Toughness [MPam^(0.5)]')
title('LEFM')

% graph for all max K values for all applied stress
figure
plot(C2_EPFM(:,1),C2_EPFM(:,3))
xlabel('Applied Stress [MPa]')
ylabel('Fracture Toughness [MPam^(0.5)]')
title('EPFM')

% tresca & von mises
figure
plot(spp,mises)
hold on
plot(spp,trescamax)
xlabel('Applied Stress [MPa]')
ylabel('Estimated Yield Stress')
legend('Von Mises','Tresca')
hold off

% safety factor plots
figure
plot(spp, SFV) % von mises safety factor
hold on
plot(convergetresca(:,1),SFT) % tresca
xlabel('Applied Stress [MPa]')
ylabel('Safety Factor')
legend('Von Mises','Tresca')
hold off

ydiag = -4e7:1e7:4e7;
xdiag = -4e7:1e7:4e7;

% plot of K (LEFM) vs Keff (EPFM)
figure
plot(K,Keff)
hold on
plot(xdiag,ydiag,'-')
xlabel('LEFM Stress Intensity Factor')
ylabel('EPFM Stress Intensity Factor')