%PreviewPlotter
%POWERED BY OCTAVE 
%The preview plotter takes the data from the PreProccessor or the simple model problem, 
%and plots previews of the structures. It also provides the ability to apply external 
%forces and single point constraints by mouse on the plot.

load PPOutput6337.txt
BarNodes;BeamNodes;Beams;Bars;BarPCs;BeamPCs;BarForces;BeamForces;

close all
figure
%Plotting Bar Structure 
subplot(1,2,1);
plot3(transpose(BarNodes(:,1)),transpose(BarNodes(:,2)),transpose(BarNodes(:,3)),'.','color','k')
hold on
Amax=max(Bars(:,4));
for i=1:rows(Bars)
  Node1=Bars(i,1);
  Node2=Bars(i,2);
  x=[BarNodes(Node1,1);BarNodes(Node2,1)];
  y=[BarNodes(Node1,2);BarNodes(Node2,2)];
  z=[BarNodes(Node1,3);BarNodes(Node2,3)];
  greyscale=0.8*(ones(1,3)-Bars(i,4)/Amax);
  line(x,y,z,'color',greyscale,'linewidth',4*Bars(i,4)/Amax)
endfor
for i=1:rows(BarForces)
  Node1=BarForces(i,1);
  quiver3(BarNodes(Node1,1),BarNodes(Node1,2),BarNodes(Node1,3),BarForces(i,2)/4000,BarForces(i,3)/4000+rand(1)/100,BarForces(i,4)/4000,'filled','color','r','maxheadsize',0.2);
  text(BarNodes(Node1,1)+1+BarForces(i,2)/4000,BarNodes(Node1,2)+1+BarForces(i,3)/4000,BarNodes(Node1,3)+1+BarForces(i,4)/4000,['F = ' num2str(norm(BarForces(i,2:end))) 'N'],'color','r')
endfor
for i=1:length(BarPCs)
  BarPC=BarPCs{i};
  plot3(BarNodes(BarPC(1,1:end-1),1),BarNodes(BarPC(1,1:end-1),2),BarNodes(BarPC(1,1:end-1),3),'x','color','b','markersize',10)
  plot3(BarNodes(BarPC(1,1:end-1),1),BarNodes(BarPC(1,1:end-1),2),BarNodes(BarPC(1,1:end-1),3),'o','color','b','markersize',10)
endfor
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
axis('equal')
title('Bar Structure Preview')
grid on
hold off

%Plotting Beam Structure 
subplot(1,2,2)
plot3(transpose(BeamNodes(:,1)),transpose(BeamNodes(:,2)),transpose(BeamNodes(:,3)),'.','color','k')
hold on
Amax=max(Beams(:,10));
for i=1:rows(Beams)
  Node1=Beams(i,1);
  Node2=Beams(i,2);
  x=[BeamNodes(Node1,1);BeamNodes(Node2,1)];
  y=[BeamNodes(Node1,2);BeamNodes(Node2,2)];
  z=[BeamNodes(Node1,3);BeamNodes(Node2,3)];
  greyscale=0.8*(ones(1,3)-Beams(i,10)/Amax);
  line(x,y,z,'color',greyscale,'linewidth',4*Beams(i,10)/Amax)
endfor
for i=1:rows(BeamForces)
  Node1=BeamForces(i,1);
  quiver3(BeamNodes(Node1,1),BeamNodes(Node1,2),BeamNodes(Node1,3),BeamForces(i,2)/4000,BeamForces(i,3)/4000+rand(1)/100,BeamForces(i,4)/4000,'filled','color','r','maxheadsize',0.2);
  if norm(BeamForces(i,2:4))~=0
    text(BeamNodes(Node1,1)+1+BeamForces(i,2)/4000,BeamNodes(Node1,2)+1+BeamForces(i,3)/4000,BeamNodes(Node1,3)+1+BeamForces(i,4)/4000,['F = ' num2str(norm(BeamForces(i,2:4))) 'N'],'color','r')
  endif
  quiver3(BeamNodes(Node1,1),BeamNodes(Node1,2),BeamNodes(Node1,3),BeamForces(i,5)/400,BeamForces(i,6)/400+rand(1)/10^4,BeamForces(i,7)/400,'filled','color','m','maxheadsize',0.2);
  if norm(BeamForces(i,5:7))~=0
    text(BeamNodes(Node1,1)+1+BeamForces(i,2)/400,BeamNodes(Node1,2)+1+BeamForces(i,3)/400,BeamNodes(Node1,3)+1+BeamForces(i,4)/400,['M = ' num2str(norm(BeamForces(i,5:7))) 'Nm'],'color','m')
  endif
endfor
for i=1:length(BeamPCs)
  BeamPC=BeamPCs{i};
  plot3(BeamNodes(BeamPC(1,1:end-1),1),BeamNodes(BeamPC(1,1:end-1),2),BeamNodes(BeamPC(1,1:end-1),3),'x','color','b','markersize',10)
  plot3(BeamNodes(BeamPC(1,1:end-1),1),BeamNodes(BeamPC(1,1:end-1),2),BeamNodes(BeamPC(1,1:end-1),3),'o','color','b','markersize',10)
endfor
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
axis('equal')
title('Beam Structure Preview')
grid on
hold off

clear all