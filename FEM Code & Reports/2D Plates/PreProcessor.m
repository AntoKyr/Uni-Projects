%CODED ON OCTAVE
%UI
clear all

load Mesh6337.txt

Elements;Nodes;LoadLines;ConstraintLine;Split4Elements;

load nii.txt
LoadType;LoadString;


  for i=1:rows(Elements)
    Node1=Elements(i,1);Node2=Elements(i,2);Node3=Elements(i,3);Node4=Elements(i,4);s=Elements(i,5);E=Elements(i,6);v=Elements(i,7);
    if hypot(Nodes(Node1,1)-Nodes(Node3,1),Nodes(Node1,2)-Nodes(Node3,2))>hypot(Nodes(Node2,1)-Nodes(Node4,1),Nodes(Node2,2)-Nodes(Node4,2))
      ElementsTrig([2*i-1 2*i],1:6)=[Node1 Node2 Node4 s E v;Node2 Node3 Node4 s E v];
    else
      ElementsTrig([2*i-1 2*i],1:6)=[Node1 Node3 Node4 s E v;Node2 Node3 Node1 s E v];
    endif
  endfor
  Elements=ElementsTrig;
  Split4Elements=[Split4Elements*2-1 Split4Elements*2];

%Forces and Constraints
for i=1:rows(LoadLines)
  ForceNodes(:,i)=find(round(Nodes(:,1)*10^7)==round(LoadLines(i,2)*10^7));
endfor
ForceNodes=transpose(ForceNodes);
ForceNodes(:,end+1)=LoadLines(:,1);
FNN=columns(ForceNodes)-1;
BCNodes=find(round(Nodes(:,1)*10^7)==round(ConstraintLine*10^7));
for i=1:rows(BCNodes)
  PCs{i}=[BCNodes(i) 0;1 0;0 0];
endfor

%Adding Constraints in symmetrical fashion
if mod(length(BCNodes),2)==0
  PCs{end+1}=[BCNodes(floor(length(BCNodes)/2)) BCNodes(floor(length(BCNodes)/2)+1) 0;0 0 0;1 1 0];
else
  PCs{end+1}=[BCNodes(ceil(length(BCNodes)/2)) 0;0 0;1 0];
endif

Forces=[];

if LoadType==1
  for i=1:rows(ForceNodes)
    T=ForceNodes(i,end)*str2num(LoadString{:});
    Fc=T/((FNN-1)*2);
    F=[transpose(ForceNodes(i,1:end-1)) ones(FNN,1)*[Fc*2 0]];
    F([1 end],2)=Fc*[1;1];
    Forces=[Forces;F];
  endfor
  
elseif LoadType==2
  for i=1:rows(ForceNodes)
    Mn=str2num(LoadString{:})*ForceNodes(i,end);
    ForceLine=ForceNodes(i,1:end-1);
    L=abs(max(Nodes(ForceLine,2))-min(Nodes(ForceLine,2)));
    Fc=3*Mn/(2*(L/2)^3);
    F=[transpose(ForceLine) Fc*Nodes(ForceLine,2)*L/(FNN-1) zeros(FNN,1)];
    F(find(Nodes(F(:,1),2)==max(Nodes(F(:,1),2))),2)=Fc*(max(Nodes(F(:,1),2))-L/(FNN-1)/2)*L/(FNN-1)/2;
    F(find(Nodes(F(:,1),2)==min(Nodes(F(:,1),2))),2)=Fc*(min(Nodes(F(:,1),2))+L/(FNN-1)/2)*L/(FNN-1)/2;
    Forces=[Forces;F];
  endfor

endif

save PPOutput6337.txt Elements Nodes Forces PCs Split4Elements

clear all


