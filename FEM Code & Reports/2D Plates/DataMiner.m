clear all
%This section gathers data for the stress con factor

if 1==0
LoadType=menu('Choose the type of load that will be applied',{'Tension','Bending'});
LoadString=inputdlg('Load','Input the wanted load');

save nii.txt LoadString LoadType
timerval=[];
MaxStresses=[];
save MaxStress.txt MaxStresses timerval
for nee=5:2:11
  tic
  load nii.txt
  n=nee
  save nii.txt LoadString LoadType n
  PiercedLamina
  PreProcessor
  Solver
  StressFactor
  load StressFact.txt
  load MaxStress.txt
  MaxStresses=[MaxStresses;Max_Stress];
  timerval=[timerval;toc];
  save MaxStress.txt MaxStresses timerval
endfor

endif
%This section plots the data
load MaxStress.txt
load nii.txt



if LoadType==1
Force=str2num(LoadString{:});
NormalStress=Force/[(0.1-0.03)*0.004];
StressConFactor=MaxStresses/NormalStress;
[StressConFactor transpose(5:2:37)]

Kt=3 - 3.13*(0.03/0.1) + 3.66*(0.03/0.1)^2 - 1.53*(0.03/0.1)^3;
figure
hold on
plot([5:2:37],StressConFactor,'b')
plot([5 40],[Kt Kt],'linestyle','--','color','r','linewidth',1)
plot([5:2:37],StressConFactor,'^','color','b')
grid on
axis([5 40 2 2.5])
ylabel('Kt','fontsize',15)
xlabel('Number of Nodes in a Quadrant','fontsize',15)
title('Stress Concentration Factor in Tension','fontsize',15)
legend({'Finite Element Approach' 'Analytical Approach'},'fontsize',15)

figure
semilogy([5:2:37],timerval,'b')
semilogy([5:2:37],timerval,'b','marker','o')
grid on

ylabel('t (s)','fontsize',30)
xlabel('Number of Nodes in a Quadrant','fontsize',30)
title('Computation Time Dependant on Number of Nodes','fontsize',30)



elseif LoadType==2
Torque=str2num(LoadString{:});
NormalStress=(6*Torque*0.03)/((0.1^3-0.03^3)*0.004);
StressConFactor=MaxStresses/NormalStress;

[StressConFactor transpose(5:2:37)]
Kta=2;
figure
hold on
plot([5:2:37],StressConFactor,'b')
plot([5 40],[Kta Kta],'linestyle','--','color','r','linewidth',1)
plot([5:2:37],StressConFactor,'^','color','b')
grid on
axis([5 40 1 2.2])
ylabel('Kta','fontsize',15)
xlabel('Number of Nodes in a Quadrant','fontsize',15)
title('Stress Concentration Factor in Bending','fontsize',15)
legend({'Finite Element Approach' 'Analytical Approach'},'fontsize',15)

figure
semilogy([5:2:37],timerval,'b')
semilogy([5:2:37],timerval,'b','marker','o')
grid on

ylabel('t (s)','fontsize',30)
xlabel('Number of Nodes in a Quadrant','fontsize',30)
title('Computation Time Dependant on Number of Nodes','fontsize',30)



endif


