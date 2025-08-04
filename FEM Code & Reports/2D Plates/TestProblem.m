%OCTAVE
clear all
E=210*10^9;
v=0.3;
Data=inputdlg({'x Nodes' 'y Nodes' 'Length(x) in m' 'Heigh(y) in m' 'Width(t) in m'},'Input data for the rectangular steel lamina. (E = 210*10^9 v=0.3)');
Nx=str2num(Data{1});
Ny=str2num(Data{2});
L=str2num(Data{3});
d=str2num(Data{4});
s=str2num(Data{5});
%Building Mesh
xStep=L/(Nx-1);
yStep=d/(Ny-1);
xvalues=repmat([0:xStep:L],Ny,1);
yvalues=repmat(transpose([-d/2:yStep:d/2]),1,Nx);
Nodacius(:,1:2:2*Nx)=xvalues;
Nodacius(:,2:2:2*Nx)=yvalues;
%Packing Nodes and elements
Nodes=NodePacker(Nodacius);
Elements=[];
for j=1:Nx-1
  for i=1:Ny-1
    Elements(end+1,:)=[j*Ny+i+1 j*Ny+i (j-1)*Ny+i (j-1)*Ny+i+1];
  endfor 
endfor


Elements=[fliplr(Elements) ones(rows(Elements),1)*[s E v]];
LoadLines=[1 L];
ConstraintLine=0;
Split4Elements=[];

save Mesh6337.txt Elements Nodes LoadLines ConstraintLine Split4Elements

plotyman=questdlg('Do you wish to view the generated mesh?','Plot Mesh');
if strcmp(plotyman,'Yes')==1
  figure
  hold on
  for i=1:rows(Nodacius)
    plot(Nodacius(i,1:2:end),Nodacius(i,2:2:end),'k')
  endfor
  for i=1:columns(Nodacius)/2
    plot(Nodacius(:,2*i-1),Nodacius(:,2*i),'k')
  endfor
  axis('equal')
  grid on
  title('Generated Mesh')
endif

clear all
