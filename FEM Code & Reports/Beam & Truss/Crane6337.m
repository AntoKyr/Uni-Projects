%Antonios-Ioakeim Kyriakopoulos 6337
%Coded on Octave
%In the preprocessor 2 cranes are "built", one structured by densely placed bar elements and one with
%thick beam elements. Both structures weigh the same. First are placed the nodes, the connecting points
%of the elements, then the elements, the constraints of the nodes that will stay in place and finally
%the external forces applied to the structures.
%Edit AEM
AEM=6337;
%Creating data
psifia = [];
while AEM>=1
  AEM=AEM/10;
  psifia=[mod(AEM,1)*10 psifia];
  AEM=floor(AEM);
endwhile
AEM1=psifia(end-3);AEM2=psifia(end-2);AEM3=psifia(end-1);AEM4=psifia(end);
F=20000;
A=1.2*(1+AEM1/10+AEM2/100);
B=1+AEM2/10+AEM3/100;
fi=deg2rad(60+AEM2+AEM4/10);
thita=deg2rad(45+AEM1+AEM3/10);
n=5+mod(AEM1,AEM4); 
L=1.5*(1+AEM4/10+AEM3/100);
A0=6*(0.5+(10*AEM4+AEM3)/100)*1/(100^2);
E=210*10^9;
G=80*10^9;

%Placing Nodes
%Base Nodes
BNodes=[0 L/2 B;0 0 B;0 -L/2 B;A L/2 0;A -L/2 0];
%Tower Nodes
TNodes=[];
for i=0:n
  TNodes=[TNodes;A+L/2 L/2 L+L*i;A+L/2 -L/2 L+L*i;A-L/2 L/2 L+L*i;A-L/2 -L/2 L+L*i];
endfor
%Towerhead nodes
HNodes=[A+L/2 -L/2 n*L+2*L;A+L/2 L/2 n*L+2*L;A+L*(3/2) 0 n*L+L*(3/2)];
%All nodes
A1Nodes=[BNodes;TNodes;HNodes];

%Rotations
fip=(pi/2-fi);
%z rotation
A2Nodes=transpose([cos(thita) -sin(thita) 0;sin(thita) cos(thita) 0;0 0 1]*transpose(A1Nodes));
%rotation around the base joint
%Select only moving nodes
MNodes=A2Nodes(6:end,:);
%Change the coordination system and rotate
Joint_Displacement=[cos(thita)*A;sin(thita)*A;0]*ones(1,rows(MNodes));
Back_z_Rotation=[cos(-thita) -sin(-thita) 0;sin(-thita) cos(-thita) 0;0 0 1];
Joint_Rotation=[cos(fip) 0 sin(fip);0 1 0;-sin(fip) 0 cos(fip)];
z_Rerotation=[cos(thita) -sin(thita) 0;sin(thita) cos(thita) 0;0 0 1];
%Perform Rotation
MNodes=transpose(z_Rerotation*(Joint_Rotation*(Back_z_Rotation*(transpose(MNodes)-Joint_Displacement)))+Joint_Displacement);
Nodes=[A2Nodes(1:5,:);MNodes];

%Placing Elements
%Thin and Thick Element cross sections
Athin=A0/2;
Athick=3*A0/2;

%Wire Elements
WBars=[1 16 E Athin;3 17 E Athin;2 8+n*4 E Athin;2 9+n*4 E Athin];

%Thick Elements
%Base
BThickBars=[4 5 E Athick;4 8 E Athick;4 6 E Athick;5 9 E Athick;5 7 E Athick;
            6 7 E Athick;7 9 E Athick;6 8 E Athick;8 9 E Athick];
%Head_
HThickBars=[n*4+6 n*4+11 E Athick;n*4+7 n*4+10 E Athick;n*4+8 n*4+11 E Athick;n*4+9 n*4+10 E Athick;n*4+10 n*4+11 E Athick
            n*4+6 n*4+12 E Athick;n*4+7 n*4+12 E Athick;n*4+10 n*4+12 E Athick;n*4+11 n*4+12 E Athick];
%Tower 
TThickBars=[];
for i=1:n
  TThickBars=[TThickBars;4*i+6 4*i+7 E Athick;4*i+7 4*i+9 E Athick;4*i+6 4*i+8 E Athick;4*i+8 4*i+9 E Athick;
                         4*i+2 4*i+6 E Athick;4*i+3 4*i+7 E Athick;4*i+4 4*i+8 E Athick;4*i+5 4*i+9 E Athick];
endfor
%All Thick Elements
ThickBars=[BThickBars;TThickBars;HThickBars];

%Thin Elements
%Base Elements
BThinBars=[5 8 E Athin;5 6 E Athin;4 9 E Athin;4 7 E Athin;6 9 E Athin;7 8 E Athin];
%Head Elements
HThinBars=[n*4+8 n*4+10 E Athin;n*4+9 n*4+11 E Athin;n*4+6 n*4+10 E Athin;n*4+7 n*4+11 E Athin];
%Tower Elements
TThinBars=[];
for i=1:n
  TThinBars=[TThinBars;i*4+3 i*4+6 E Athin;i*4+4 i*4+6 E Athin;
                       i*4+2 i*4+7 E Athin;i*4+5 i*4+7 E Athin;
                       i*4+3 i*4+9 E Athin;i*4+4 i*4+9 E Athin;
                       i*4+5 i*4+8 E Athin;i*4+2 i*4+8 E Athin;
                       i*4+7 i*4+8 E Athin;i*4+6 i*4+9 E Athin]; 
endfor
%All Thin Elements 
ThinBars=[BThinBars;TThinBars;HThinBars];
%All Bar ELements
Bars=[ThinBars;ThickBars;WBars];

%V1=0;
for i=1:rows(Bars)
 % V1=V1+Bars(i,4)*norm(Nodes(Bars(i,1))-Nodes(Bars(i,2)));
endfor
%Beam Structure Elements
%Beam structure cross sectional area calculation
Lthin=0;Lthick=0;
for i=1:rows(ThinBars)
  Node1=ThinBars(i,1);
  Node2=ThinBars(i,2);
  Lthin=Lthin+norm(Nodes(Node2,:)-Nodes(Node1,:));
endfor
for i=1:rows(ThickBars)
  Node1=ThickBars(i,1);
  Node2=ThickBars(i,2);
  Lthick=Lthick+norm(Nodes(Node2,:)-Nodes(Node1,:));
endfor
Abeam=Athick+Athin*Lthin/Lthick;
%Beam Elements
I=Abeam^2/(4*pi);
%Due to the cylindrical cross section the vector of rotation is unimportant
WBeams=[WBars(:,1:2) ones(rows(WBars),1)*[rand(1,3) E 0 0 0 Athin]];
Beams=[ThickBars(:,1:2) ones(rows(ThickBars),1)*[rand(1,3) E G I I Abeam];WBeams];
          
%Point Constraints
%SPCs
for i=1:4
  BarPCs{3*i-2}=[i 0;1 0;0 0;0 0];
  BarPCs{3*i-1}=[i 0;0 0;1 0;0 0];
  BarPCs{3*i}=[i 0;0 0;0 0;1 0];
endfor
BarPCs{13}=[5 0;0 0;0 0;1 0];
%MPCs
BarPCs{14}=[5 0;cos(thita) 0;sin(thita) 0;0 0];

%SPCs
for i=1:3
  BeamPCs{6*i-5}=[i 0;1 0;0 0;0 0;0 0;0 0;0 0];
  BeamPCs{6*i-4}=[i 0;0 0;1 0;0 0;0 0;0 0;0 0];
  BeamPCs{6*i-3}=[i 0;0 0;0 0;1 0;0 0;0 0;0 0];
  BeamPCs{6*i-2}=[i 0;0 0;0 0;0 0;1 0;0 0;0 0];
  BeamPCs{6*i-1}=[i 0;0 0;0 0;0 0;0 0;1 0;0 0];
  BeamPCs{6*i}=[i 0;0 0;0 0;0 0;0 0;0 0;1 0];
endfor
BeamPCs{19}=[4 0;1 0;0 0;0 0;0 0;0 0;0 0];
BeamPCs{20}=[4 0;0 0;1 0;0 0;0 0;0 0;0 0];
BeamPCs{21}=[4 0;0 0;0 0;1 0;0 0;0 0;0 0];
BeamPCs{22}=[5 0;0 0;0 0;1 0;0 0;0 0;0 0];
%MPCs
BeamPCs{23}=[5 0;cos(thita) 0;sin(thita) 0;0 0;0 0;0 0;0 0];


%External Forces
%Bar Forces [Node fx fy fz]
BarForces=[4*n+12 0 0 -F];

%Beam Forces [Node x y z Mx My Mz]
BeamForces=[4*n+12 0 0 -F 0 0 0];

%Splitting Nodes
BarNodes=Nodes;
BeamNodes=Nodes;

%Outputs
%Nodes, Bars, Beams, External Forces, MPCs and SPCs 
save PPOutput6337.txt BarNodes BeamNodes Beams Bars BarPCs BeamPCs BarForces BeamForces 

clear all
