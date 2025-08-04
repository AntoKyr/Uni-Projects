%Antonios-Ioakeim Kyriakopoulos 6337
%Coded on Octave
%In the preprocessor 2 pylons are "built", one structured by densely placed bar elements and one with
%sparcer beam elements. First are placed the nodes, the connecting points
%of the elements, then the elements, the constraints of the nodes and finally the external forces
%applied to the structures.

%Properties
E=210*10^9;
G=80*10^9;
A0=5.64*10^-4;
L=4.2;
n=3; %the number of "floors" on the pylon is a variable parameter. Just for fun.

%Adding Nodes
%Base
BaseNodes=[-L -L 0;-L L 0;L L 0;L -L 0];
%Tower 
TowerNodes=[repmat([BaseNodes(:,1:2)/2],n,1) kron(transpose(1:n),ones(4,1))*L];
%Head
HeadNodes=[[BaseNodes(:,1:2)/2 ones(4,1)*(n+0.75)*L];[1.5*[-L 0;0 L;L 0;0 -L] ones(4,1)*n*L]];
%All Nodes
Nodes=[BaseNodes;TowerNodes;HeadNodes];
BarNodes=Nodes;
BeamNodes=Nodes;

%Adding Elements
%Bar Elements
%Base
BaseThickBars=[transpose(1:4) transpose(5:8)];
BaseThinBars=[1 6;1 8;2 7;2 5;3 8;3 6;4 7;4 5];
%Tower
ThickHBars=transpose([5:8+n*4;6:9+n*4]);
ThickHBars(4:4:end,2)=ThickHBars(4:4:end,2)-4;
TowerThickBars=[[transpose(5:4+4*n) transpose(9:8+4*n)];ThickHBars];
ThinHBars=[transpose([5:4:5+n*4;7:4:7+n*4]);transpose([6:4:6+n*4;8:4:8+n*4])];
ThinV1Bars=transpose([5:4+n*4;10:9+n*4]);
ThinV1Bars(4:4:end,2)=ThinV1Bars(4:4:end,2)-4; 
ThinV2Bars=transpose([5:4+4*n;8:7+4*n]);
ThinV2Bars(1:4:end,2)=ThinV2Bars(1:4:end,2)+4;
TowerThinBars=[ThinHBars;ThinV1Bars;ThinV2Bars];
%Head
numb=transpose(1:3);
HeadBars=[numb+4*n numb+8+4*n;numb+1+4*n numb+8+4*n;numb+4+4*n numb+8+4*n;numb+5+4*n numb+8+4*n;1+4*n 12+4*n;4+4*n 12+4*n;5+4*n 12+4*n;8+4*n 12+4*n];

ThickBars=[BaseThickBars;TowerThickBars;HeadBars];
ThinBars=[BaseThinBars;TowerThinBars];
%Adding Properties to Elements
%Bars
Bars=[[ThickBars ones(rows(ThickBars),1)*[E 1.5*A0]];[ThinBars ones(rows(ThinBars),1)*[E 0.5*A0]]];
%Beams
Ithick=((A0*1.5)^2)/(4*pi);
Beams=[ThickBars ones(rows(ThickBars),1)*[rand(1,3) E G Ithick Ithick 1.5*A0]];

%The 2 next lines exist just for fun. They are the diagonal elements missing from the beam structure. They only exist for comparison with the bar structure. 
Ithin=(A0*0.5)^2/(4*pi);
Beams=[Beams;[ThinBars ones(rows(ThinBars),1)*[rand(1,3) E G Ithin Ithin 0.5*A0]]];

%Adding BCs
%Bars
%Constraining Nodes 
%Each node is constrained in all 3 axes. 
for i=1:4
  BarPCs{3*i-2}=[i 0;1 0;0 0;0 0];
  BarPCs{3*i-1}=[i 0;0 0;1 0;0 0];
  BarPCs{3*i}=[i 0;0 0;0 0;1 0];
endfor
%Beams
for i=1:4
  BeamPCs{6*i-5}=[i 0;1 0;0 0;0 0;0 0;0 0;0 0];
  BeamPCs{6*i-4}=[i 0;0 0;1 0;0 0;0 0;0 0;0 0];
  BeamPCs{6*i-3}=[i 0;0 0;0 0;1 0;0 0;0 0;0 0];
  BeamPCs{6*i-2}=[i 0;0 0;0 0;0 0;1 0;0 0;0 0];
  BeamPCs{6*i-1}=[i 0;0 0;0 0;0 0;0 0;1 0;0 0];
  BeamPCs{6*i}=[i 0;0 0;0 0;0 0;0 0;0 0;1 0];
endfor       

%Adding External Forces
%Bars
BarForces=[11+4*n 0 8000 -15000;10+4*n 8000 0 -9000;9+4*n 0 8000 -11000;12+4*n 8000 0 -7000];
%Beams
BeamForces=[BarForces zeros(4,3)];

%Outputs
%Nodes, Bars, Beams, External Forces, MPCs and SPCs 
save PPOutput6337.txt BarNodes BeamNodes Beams Bars BarPCs BeamPCs BarForces BeamForces 

clear all