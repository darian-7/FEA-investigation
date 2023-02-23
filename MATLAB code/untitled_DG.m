clc
clear
%material properties and given values
E=70*10^3;
v=0.3;
yield_stress=400;
K_c=20;
theta_a=0;
theta_b=45;
theta_c=90;
theta_p=30;
strain_a=fn_strain(2049155,theta_a);
strain_b=fn_strain(2049155,theta_b);
strain_c=fn_strain(2049155,theta_c);
strain_p=fn_strain(2049155,theta_p);

%calculate stress
sig_xx=(E/(1-v^2))*(strain_a+v*strain_c);
sig_yy=(E/(1-v^2))*(strain_c+v*strain_a);
sig_xy=strain_b*(E/(1-v));

%obtain principal stress 
%p1=((stress_xx+stress_yy)/2)+sqrt(((stress_xx-stress_yy)/2)^2+stress_xy^2);
%p2=((stress_xx+stress_yy)/2)-sqrt(((stress_xx-stress_yy)/2)^2+stress_xy^2);

%original stress tensor
sig_ij=[sig_xx sig_xy; sig_xy sig_yy];

%rotate stress tensor
Q=[cos(theta_p) sin(theta_p);-sin(theta_p) cos(theta_p)];
sig_30=Q*sig_ij*transpose(Q);

%%

a = 0.0025;

for i = 0:1:180
        k = i+1;
        %define angle phi 
        phi = (i*pi)/180;
        
        R{k} = [(cos(phi))^2, (sin(phi))^2, 2*cos(phi)*sin(phi);
            (sin(phi))^2, (cos(phi))^2, -2*cos(phi)*sin(phi);
            -cos(phi)*sin(phi), cos(phi)*sin(phi), (cos(phi))^2 - (sin(phi))^2];
        
        omega(k) = i;
end


for i = -400:1:400
    j = i + 401;
    sigmap(j) = i;
    
    % define stress tensor of extra stress in rotated system 
    stressptensor{j} = [sigmap(j), 0; 0, 0];
    
    % total stress tensor
    stot{j} = Q + sig_ij;
    stotvec{j} = [stot{j}(1,1); stot{j}(2,2); stot{j}(1,2)];
    
    % new stress tensor for crack rotation
    for k = 1:1:181
        
        svec{j,k} = R{k}*stotvec{j};
        syy = svec{j,k}(2);
        
        if syy > 0
            tensilesyy(j,k) = syy;
        else
            tensilesyy(j,k) = 0;
        end
      
        KIC(j,k) = tensilesyy(j,k)*sqrt(pi*a);
        
        Keff(j,k) = (tensilesyy(j,k)*sqrt(pi*a))./(sqrt(1-0.5*((tensilesyy(j,k))/(yield_stress)).^2));
        KICnorm(j,k) = 0.05*tensilesyy(j,k)*sqrt(pi*a);
        
        SFKIC(j,k) = 1/KICnorm(j,k);
          
    end
    % maximum KICnorm for varying sp. Gives max KICnorm for varying stress
    % and lets me know at what angle
    [KICnormaxs(j,2), KICnormaxs(j,1)] = max(KICnorm(j,:));
    
end

SFKIC(isinf(SFKIC)|isnan(SFKIC)) = 0;

%maximum tensile stress for varying angles
SYmaxangle = max(tensilesyy,[],1);
%maximum tensile stress for varying sp
SYmaxstress = max(tensilesyy,[],2);

% plot stress sp and angle that cause the worst effects
figure

plot(KICnormaxs(:,1),sigmap)

% Figure 1 - plot max stress with angle
figure(1)
hold on
box on
grid on

plot(omega,SYmaxangle,'-','Color',[0, 0.4470, 0.7410],'LineWidth',1.5)

plot(omega,0*omega + 400,'Color',[0.9290 0.6940 0.1250],'LineWidth',1.2)
xlimits=[0,180];
xlim(xlimits)
ylimits=[0,500];
ylim(ylimits)

xlabel('angle in deg')
ylabel('max tensile stress')

set(gca,'XColor','k', 'YColor','k','FontSize',12,'XTick',0:30:180,'YTick',0:100:500)

%%

p=-400:1:400; %range of sigma p

%find principal stresses
for j=1:length(p);
    sig_p1(j)=((sig_30(1,1)+sig_30(2,2)+p(j))/2)+sqrt(((sig_30(1,1)+p(j)-sig_30(2,2))/2)^2+sig_30(1,2)^2);
    sig_p2(j)=((sig_30(1,1)+sig_30(2,2)+p(j))/2)-sqrt(((sig_30(1,1)+p(j)-sig_30(2,2))/2)^2+sig_30(1,2)^2);
end

%von mises' stress
for j=1:length(p);
VM(j)=(1/sqrt(2))*sqrt(((sig_p1(j)-sig_p2(j))^2)+((sig_p2(j))^2)+((sig_p1(j))^2))
end

%tresca stress
for j=1:length(p);
    av1=abs(sig_p1(j)-sig_p2(j));
    av2=abs(sig_p2(j));
    av3=abs(-sig_p1(j));
   tres=[av1 av2 av3];
   tresca(j)=max(tres);
end

SF_VM=[];
%safety factor
for j=1:length(p)
SF_VM(j)=yield_stress/VM(j);
SF_Tres(j)=yield_stress/tresca(j);
end
maxtres=max(SF_Tres);
maxvm=max(SF_VM);


%% crack analysis

%rotate stress tensor
sig_new = [];

a=0.0025;
sqr=sqrt(pi*a);

%all orientation for every sigma p
for p = -400:1:400
    i = p + 401;
    
    sig_p{i} = [sig_30(1,1) + p, sig_30(1,2); sig_30(2,1), sig_30(2,2)];
    
    for t = 1:1:181
         rotate{t} = [cosd(t-1), sind(t-1); -sind(t-1), cosd(t-1)];
         rotatetensor{i,t} = rotate{t}*sig_p{i}*transpose(rotate{t});
%          max_K(i,t) = rotatetensor{i,t}(2,2);
    end 
end


% max k for every sigma p at every orientation
for i = 1:801
    for j = 1:181

        max_K(i,j) = rotatetensor{i,j}(2,2)*sqr;

    end
end

%max k for every sigma p at 30 deg
k_30=max_K(:,31);


%% Plotting

% 3D plot
figure      
mesh(max_K)
xlabel('Angle(degrees)')
ylabel('sigma p (MPa)')
zlabel ('K')

% 2D contour plot
sig_p_plot = [-400*10^6:10^6:400*10^6];
orien_angle = [0:1:180];
[X,Y] = meshgrid(orien_angle,sig_p_plot);
figure
c = contourf(X,Y,max_K);
hold on
a = colorbar
xlabel('Crack orientation (degrees)');
ylabel('Additional stress "P" (MPa)');
a.Label.String = 'Stress Intensity Factor (MPa m^0.5)';
hold off

% Tresca and Von Mises plot
tresca_SF = transpose(SF_Tres);
vonmises_SF = transpose(SF_VM);
figure
y1 = tresca_SF;
plot(sig_p_plot,y1);
hold on
y2 = vonmises_SF;
plot(sig_p_plot,y2);
hold off
xlabel('Additional stress "P" (MPa)');
ylabel('Safety factor');
yline(2);
legend('Tresca stress','Von Mises stress')


