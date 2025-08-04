%Stress Factor
%Coded on Octave
load SolverOutput6337.txt
Elements;Nodes;Forces;PCs;Displacements;Deformations;Stresses;PrincStresses;thita;

XStresses=Stresses{2};

Node1=find(round(Nodes(:,1)*10^8)==0);

[Waste Node2]=min(abs(Nodes(Node1,2)));
Node1=Node1(Node2);

[Element1 Waste]=find(Elements==Node1);
Element1=Element1(1);

Max_Stress=abs(XStresses(Element1));

save StressFact.txt Max_Stress

clear all



