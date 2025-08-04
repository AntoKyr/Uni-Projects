
clear all
M=[1.25:0.25:5];
thita_v=deg2rad([0:0.2:70]);
gamma=1.4;

%Solutions
ColorM=jet(length(M));
figure
hold on
colormap('Jet')
counter=0;
delta=0;
thita_s=deg2rad([1:0.5:90]);
for M1=M
  %Wedge Shock
  counter=counter+1;
  d=atan(2*cot(thita_s).*(M1^2*sin(thita_s).^2-1)./(M1^2*(gamma+cos(2*thita_s))+2));
  M2=sqrt(((sin(thita_s)*M1).^2+2/(gamma-1))./((2*gamma/(gamma-1))*(sin(thita_s)*M1).^2-1))./sin(thita_s-d);
  %Cone Shock
  if 1==0
    AngleRange=find(d>0);
  else
    AngleRange=1:length(d);
  endif
  Vn=(2./((gamma-1)*M2.^2)+1).^(-1/2);
  Vthb=-sin(thita_s(AngleRange)-d(AngleRange)).*Vn(AngleRange);
  Vrb=cos(thita_s(AngleRange)-d(AngleRange)).*Vn(AngleRange);
  WODE=waitbar(0,'Solving ODEs');
  for i=1:length(Vthb)
   waitbar(i/length(Vthb),WODE);
   thitas=linspace(thita_s(i),0,10^5);
   [thitas Vs]=ode15s("conedif",thitas,[Vrb(i) Vthb(i)],odeset('RelTol',30));
   Vth=Vs(:,2);
   [Vthc,spot]=min(abs(Vth));
   thita_c(i)=thitas(spot);
  endfor
  close (WODE)

  %Post-Processing Cone solutions
  AngleRange=find(d>0);
  %This loop trims off extreme values which arise due to the very high instability of the equations
  for i=1:3
    %Fitting data into plynomial
    Pc=polyfit(thita_s(AngleRange),thita_c(AngleRange),10);
    %Calculating std
    StnDev=std(thita_c(AngleRange)-polyval(Pc,thita_s(AngleRange)));
    %Removing extreme values
    AngleRange=AngleRange(find(abs(thita_c(AngleRange)-polyval(Pc,thita_s(AngleRange)))<i*StnDev));
  endfor

  %Cone plot
  subplot(1,11,6:11)
  hold on
  %plot(rad2deg(thita_c(AngleRange)),rad2deg(thita_s(AngleRange)),'color',ColorM(counter,:))
  plot(rad2deg([polyval(Pc,thita_s(1:end-1)) 0]),rad2deg(thita_s),'color',ColorM(counter,:))
  %Wedge plot
  subplot(1,11,1:5)
  hold on
  plot(rad2deg(d),rad2deg(thita_s),'color',ColorM(counter,:))
  hold off
endfor

subplot(1,11,1:5)
axis([0 45 0 90])
title('θ-β-M Diagram of 2D Wedge Shock')
xlabel('Deflection angle θ, degrees')
ylabel('Shock wave angle β, degrees')
grid on

subplot(1,11,6:11)
axis([0 60 0 90])
title('θ-β-M Diagram of Cone Shock')
xlabel('Cone half-angle θc, degrees')
ylabel('Shock wave angle β, degrees')
grid on

%ColorBar
FringNums=linspace(1.25,5,6);
for i=1:6
  FringeString{i}=num2str(FringNums(i));
endfor
FringeBar=colorbar('yticklabel',FringeString,'ytick',linspace(0,1,6));
title(FringeBar,'Mach Number')






