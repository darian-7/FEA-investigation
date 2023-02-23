clc
clear
%material properties and given values
E=70*10^9;
v=0.3;
yield_stress=400*10^6;
K_c=20;
theta_a=0;
theta_b=45;
theta_c=90;
theta_p=30;
strain_a=fn_strain(1944637,theta_a);
strain_b=fn_strain(1944637,theta_b);
strain_c=fn_strain(1944637,theta_c);
strain_p=fn_strain(1944637,theta_p);

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


p=-400*10^6:10^6:400*10^6; %range of p

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

%safety factor
for j=1:length(p)
SF_VM(j)=yield_stress/VM(j);
SF_Tres(j)=yield_stress/tresca(j);
end
maxtres=max(SF_Tres);
maxvm=max(SF_VM);


%%crack analysis

%rotate stress tensor
sig_new = [];
a=0.0025;
sqr=sqrt(pi*a);
i = 1;

for sigmap = -400*10^6:10^6:400*10^6
for theta = 0:pi/180:2*pi
        sigmatensor = [sig_ij(1,1) + sigmap, sig_ij(1,2); sig_ij(2,1) sig_ij(2,2)];
        rotate = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        transformedtensor = rotate*sigmatensor*transpose(rotate);
        sig_new(i,1) = sigmap;
        sig_new(i,2) = theta;
        sig_new(i,3) = transformedtensor(2,2);
        sig_new(i,4) = transformedtensor(2,2)*sqr;
        i = i+1;
     end
end 
% 
% kval={};
% p_value= -400*10^6;
% %p_value=transpose(p_value);
% k=1;
% j=1;
% 
% while p_value <= 400*10^6
%     for i=1:length(sig_new)
%         if sig_new(i,1)==p_value
%             kval{k}(j,1)=p_value;
%             kval{k}(j,2)=sig_new(i,2);
%             kval{k}(j,3)=sig_new(i,4);
%         else
%             continue
%         end
%         j=j+1;
%     end
%     j=1;
%     p_value=p_value + 10^6;
%     k=k+1;
% end
% 
% %max_kval=cellfun(max{i,3},max_k);
% 
% max_k=[];
% i=1;
% for i=1:length(kval)
%     max_k(i,1)=max(kval{i}(:,1)); % sigma_p
%     [max_k(i,2),I]=max(kval{i}(:,3)); % K
%     max_k(i,3) = kval{i}(I,2); %angle in rad
%     i=i+1;
% end
% 
% x=max_k(:,1);
% y=max_k(:,3)*180/pi;
% z=max_k(:,2);
% contourf(x,y,z)
