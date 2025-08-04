%PostProcessor
%Antonios-Ioakeim Kyriakopoulos 6337
%Coded on Octave
%Ôhe post processor uses the outputs of the solver, to help visualize and better understand the behaviour
%of the structures. 

load SolverOutput6337.txt
BarNodes;BeamNodes;Beams;Bars;BarPCs;BeamPCs;BarForces;BeamForces;BarDisplacements;BeamDisplacements;BarStress;

ColorM=jet(1001);
%Bar Structure Visualization
%Plotting Bar Structure Stresses
figure
hold on
colormap('jet')
MaxStress=max([abs(max(BarStress)) abs(min(BarStress))]);
Amax=max(Bars(:,4));
Loading=waitbar(0,'Plotting Element Stresses');
div=rows(Bars);
for i=1:rows(Bars)
  waitbar(i/div,Loading)
  Node1=Bars(i,1);
  Node2=Bars(i,2);
  A=Bars(i,4);
  x=[BarNodes(Node1,1);BarNodes(Node2,1)];
  y=[BarNodes(Node1,2);BarNodes(Node2,2)];
  z=[BarNodes(Node1,3);BarNodes(Node2,3)];
  line(x,y,z,'linewidth',3*A/Amax,'color',ColorM(round(500*BarStress(i)/MaxStress+500)+1,:))
endfor
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
axis('equal')
title('Bar Structure Stresses')
grid on
FringeBar=colorbar('ytick',linspace(0,1,7),'yticklabel',linspace(-MaxStress/10^6,MaxStress/10^6,7));
title(FringeBar,'Mpa')
hold off
close(Loading)

%Plotting Deformed Bar Structure
figure
hold on
colormap('jet')
ExagFactor=100;
plot3(transpose(BarNodes(:,1)),transpose(BarNodes(:,2)),transpose(BarNodes(:,3)),'o','color',[0.5 0.5 0.5])
BarDisplacementNorm=hypot(hypot(BarDisplacements(:,1),BarDisplacements(:,2)),BarDisplacements(:,3));
MaxDisplacement=max(BarDisplacementNorm);
Loading=waitbar(0,'Plotting Deformed Structure');
Discret=10;
div=rows(Bars)*Discret;
for i=1:rows(Bars)
  Node1=Bars(i,1);Node2=Bars(i,2);A=Bars(i,4);
  x=[BarNodes(Node1,1)+BarDisplacements(Node1,1)*ExagFactor;BarNodes(Node2,1)+BarDisplacements(Node2,1)*ExagFactor];
  y=[BarNodes(Node1,2)+BarDisplacements(Node1,2)*ExagFactor;BarNodes(Node2,2)+BarDisplacements(Node2,2)*ExagFactor];
  z=[BarNodes(Node1,3)+BarDisplacements(Node1,3)*ExagFactor;BarNodes(Node2,3)+BarDisplacements(Node2,3)*ExagFactor];
  xx=linspace(x(1),x(2),Discret+1);yy=linspace(y(1),y(2),Discret+1);zz=linspace(z(1),z(2),Discret+1);
  InterpDisplacement=round(1000*(([Discret:-1:1]*BarDisplacementNorm(Node1)+[1:Discret]*BarDisplacementNorm(Node2))/(Discret+1))/MaxDisplacement)+1;
  for j=1:Discret
    waitbar(((i-1)*Discret+j)/div,Loading)
    line(xx(j:j+1),yy(j:j+1),zz(j:j+1),'linewidth',2*A/Amax,'color',ColorM(InterpDisplacement(j),:))  
  endfor
  x=[BarNodes(Node1,1);BarNodes(Node2,1)];
  y=[BarNodes(Node1,2);BarNodes(Node2,2)];
  z=[BarNodes(Node1,3);BarNodes(Node2,3)];
  line(x,y,z,'linestyle','--','color',[0.5 0.5 0.5])
endfor  
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
axis('equal')
title(['Bar Structure Displacements Exagerated x' num2str(ExagFactor)])
grid on
FringeBar=colorbar('ytick',linspace(0,1,7),'yticklabel',linspace(0,MaxDisplacement*10^3,7));
title(FringeBar,'mm')
hold off
close(Loading)  
  
%Plotting Reaction forces and final structure view
figure
plot3(transpose(BarNodes(:,1)),transpose(BarNodes(:,2)),transpose(BarNodes(:,3)),'.','color','k')
hold on
Amax=max(Bars(:,4));
Loading=waitbar(0,'Plotting Final Structure View');
div=rows(Bars)+rows(BarForces)+length(BarPCs);
for i=1:rows(Bars)
  waitbar(i/div,Loading)
  Node1=Bars(i,1);
  Node2=Bars(i,2);
  x=[BarNodes(Node1,1);BarNodes(Node2,1)];
  y=[BarNodes(Node1,2);BarNodes(Node2,2)];
  z=[BarNodes(Node1,3);BarNodes(Node2,3)];
  greyscale=0.8*(ones(1,3)-Bars(i,4)/Amax);
  line(x,y,z,'color',greyscale,'linestyle','--')
endfor
MaxForce=max(hypot(hypot(BarForces(:,2),BarForces(:,3)),BarForces(:,4)));
for i=1:rows(BarForces)
  waitbar((i+rows(Bars))/div,Loading)
  xcord=BarNodes(BarForces(i,1),1);
  ycord=BarNodes(BarForces(i,1),2);
  zcord=BarNodes(BarForces(i,1),3);
  if (norm(BarForces(i,2:4))>10^-9*MaxForce)&&(norm(BarForces(i,2:4))>10^-9)
    quiver3(xcord,ycord,zcord,BarForces(i,2)/8000,BarForces(i,3)/8000+rand(1)/100,BarForces(i,4)/8000,'filled','color','r','maxheadsize',0.15);
    text(xcord+BarForces(i,2)/8000,ycord+BarForces(i,3)/8000,zcord+BarForces(i,4)/8000,['F = ' num2str(norm(BarForces(i,2:end))) 'N'],'color','r')
  endif
endfor
for i=1:length(BarPCs)
  waitbar((i+rows(Bars)+rows(BarForces))/div,Loading)
  BarPC=BarPCs{i};
  xcord=BarNodes(BarPC(1,1:end-1),1);
  ycord=BarNodes(BarPC(1,1:end-1),2);
  zcord=BarNodes(BarPC(1,1:end-1),3);
  plot3(xcord,ycord,zcord,'x','color','b','markersize',10)
  plot3(xcord,ycord,zcord,'o','color','b','markersize',10)
endfor
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
axis('equal')
title('Bar Structure External Forces')
grid on
hold off 
close(Loading)
  
  
%Beam Structure Visualization
%Plotting Beam Structure Stresses
%Plotting Deformed Beam Structure
figure
hold on
colormap('jet')
ExagFactor=100;
plot3(transpose(BeamNodes(:,1)),transpose(BeamNodes(:,2)),transpose(BeamNodes(:,3)),'o','color',[0.5 0.5 0.5])
BeamDisplacementNorm=hypot(hypot(BeamDisplacements(:,1),BeamDisplacements(:,2)),BeamDisplacements(:,3));
MaxDisplacement=max(BeamDisplacementNorm);
Loading=waitbar(0,'Plotting Deformed Beam Structure');
Discret=10;
div=rows(Beams)*Discret;
for i=1:rows(Beams)
  Node1=Beams(i,1);Node2=Beams(i,2);A=Beams(i,10);
  x=[BeamNodes(Node1,1)+BeamDisplacements(Node1,1)*ExagFactor;BeamNodes(Node2,1)+BeamDisplacements(Node2,1)*ExagFactor];
  y=[BeamNodes(Node1,2)+BeamDisplacements(Node1,2)*ExagFactor;BeamNodes(Node2,2)+BeamDisplacements(Node2,2)*ExagFactor];
  z=[BeamNodes(Node1,3)+BeamDisplacements(Node1,3)*ExagFactor;BeamNodes(Node2,3)+BeamDisplacements(Node2,3)*ExagFactor];
  xx=linspace(x(1),x(2),Discret+1);yy=linspace(y(1),y(2),Discret+1);zz=linspace(z(1),z(2),Discret+1);
  InterpDisplacement=round(1000*(([Discret:-1:1]*BeamDisplacementNorm(Node1)+[1:Discret]*BeamDisplacementNorm(Node2))/(Discret+1))/MaxDisplacement)+1;
  for j=1:Discret
    waitbar(((i-1)*Discret+j)/div,Loading)
    line(xx(j:j+1),yy(j:j+1),zz(j:j+1),'linewidth',2*A/Amax,'color',ColorM(InterpDisplacement(j),:))  
  endfor
  x=[BeamNodes(Node1,1);BeamNodes(Node2,1)];
  y=[BeamNodes(Node1,2);BeamNodes(Node2,2)];
  z=[BeamNodes(Node1,3);BeamNodes(Node2,3)];
  line(x,y,z,'linestyle','--','color',[0.5 0.5 0.5])
endfor  
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
axis('equal')
title(['Beam Structure Displacements Exagerated x' num2str(ExagFactor)])
grid on
FringeBar=colorbar('ytick',linspace(0,1,7),'yticklabel',linspace(0,MaxDisplacement*10^3,7));
title(FringeBar,'mm')
hold off
close(Loading)  
  
%Plotting Reaction forces and final structure view
figure
plot3(transpose(BeamNodes(:,1)),transpose(BeamNodes(:,2)),transpose(BeamNodes(:,3)),'.','color','k')
hold on
Amax=max(Beams(:,10));
Loading=waitbar(0,'Plotting Final Structure View');
div=rows(Beams)+rows(BeamForces)+length(BeamPCs);
for i=1:rows(Beams)
  waitbar(i/div,Loading)
  Node1=Beams(i,1);
  Node2=Beams(i,2);
  x=[BeamNodes(Node1,1);BeamNodes(Node2,1)];
  y=[BeamNodes(Node1,2);BeamNodes(Node2,2)];
  z=[BeamNodes(Node1,3);BeamNodes(Node2,3)];
  greyscale=0.8*(ones(1,3)-Beams(i,10)/Amax);
  line(x,y,z,'color',greyscale,'linestyle','--')
endfor
MaxForce=max(hypot(hypot(BeamForces(:,2),BeamForces(:,3)),BeamForces(:,4)));
MaxMoment=max(hypot(hypot(BeamForces(:,5),BeamForces(:,6)),BeamForces(:,7)));
for i=1:rows(BeamForces)
  waitbar((i+rows(Beams))/div,Loading)
  Node1=BeamForces(i,1);
  xcord=BeamNodes(Node1,1);
  ycord=BeamNodes(Node1,2);
  zcord=BeamNodes(Node1,3);
  if (norm(BeamForces(i,2:4))>10^-9*MaxForce)&&(norm(BeamForces(i,2:4))>10^-9)
    quiver3(xcord,ycord,zcord,BeamForces(i,2)/8000,BeamForces(i,3)/8000+rand(1)/10^5,BeamForces(i,4)/8000,'filled','color','r','maxheadsize',0.15);
    text(BeamNodes(Node1,1)+BeamForces(i,2)/8000,BeamNodes(Node1,2)+BeamForces(i,3)/8000,BeamNodes(Node1,3)+BeamForces(i,4)/8000,['F = ' num2str(norm(BeamForces(i,2:4))) 'N'],'color','r')
  endif
  if (norm(BeamForces(i,5:7))>10^-9*MaxMoment)&&(norm(BeamForces(i,5:7))>10^-9)
    quiver3(BeamNodes(Node1,1),BeamNodes(Node1,2),BeamNodes(Node1,3),BeamForces(i,5)/8000,BeamForces(i,6)/8000+rand(1)/10^4,BeamForces(i,7)/8000,'filled','color','m','maxheadsize',0.15);
    text(BeamNodes(Node1,1)+BeamForces(i,5)/8000,BeamNodes(Node1,2)+BeamForces(i,6)/8000,BeamNodes(Node1,3)+BeamForces(i,7)/8000,['M = ' num2str(norm(BeamForces(i,5:7))) 'Nm'],'color','m')
  endif
endfor
for i=1:length(BeamPCs)
  waitbar((i+rows(Beams)+rows(BeamForces))/div,Loading)
  BeamPC=BeamPCs{i};
  xcord=BeamNodes(BeamPC(1,1:end-1),1);
  ycord=BeamNodes(BeamPC(1,1:end-1),2);
  zcord=BeamNodes(BeamPC(1,1:end-1),3);
  plot3(xcord,ycord,zcord,'x','color','b','markersize',10)
  plot3(xcord,ycord,zcord,'o','color','b','markersize',10)
endfor
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
axis('equal')
title('Beam Structure Forces and Review')
grid on
hold off 
close(Loading)

clear all