%Octave
%Mesh Generator for pierced lamina
clear all
%Dimensions
w=0.170;
a=0.100;
r=0.015;
s=0.004;
E=210*10^9;
v=0.3;

load nii.txt
n;


%Square Points
Sqr=transpose([linspace(-a/2,a/2,n) (a/2)*ones(1,n-2) linspace(a/2,-a/2,n) -(a/2)*ones(1,n-2);(a/2)*ones(1,n) linspace(a/2-a/n,-a/2+a/n,n-2) -(a/2)*ones(1,n) linspace(-a/2+a/n,a/2-a/n,n-2)]);
%Circle Points
Crcl=transpose([r*cos(linspace(3*pi/4,-5*pi/4+2*pi/(4*n-4),4*n-4));r*sin(linspace(3*pi/4,-5*pi/4+2*pi/(4*n-4),4*n-4))]);
%Internal Nodes
m=n;
k=n;
FractS=([0:m].*[[0:m]+k*ones(1,m+1)])/max([0:m].*[[0:m]+k*ones(1,m+1)]);
for j=1:4*n-4
Radiance(:,2*j-1)=FractS*(Sqr(j,1)-Crcl(j,1))+Crcl(j,1);
Radiance(:,2*j)=FractS*(Sqr(j,2)-Crcl(j,2))+Crcl(j,2);
endfor
%Calculating Border Nodes
Radiance=Radiance(1:ceil(rows(Radiance)/pi),:);
Radiance=[Radiance Radiance(:,1:2)];
SqrCrclSide=2*Radiance(end,2);
xNodeNum=2*floor((n-1)*w/SqrCrclSide/2)+mod(n,2);
xNodeSpace=w/(xNodeNum-1);
yNodeNum=2*floor((n-1)*a/SqrCrclSide/2)+mod(n,2);
yNodeSpace=a/(yNodeNum-1);

%Placing Cross Limbs
xMidBorderNodes(1:2:2*n)=linspace(-(n-1)*xNodeSpace/2,(n-1)*xNodeSpace/2,n);
xMidBorderNodes(2:2:2*n)=(a/2)*ones(1,n);
TopLimb=transpose(linspace(Radiance(end,1:2*n),xMidBorderNodes,(yNodeNum-n)/2+1));
BotLimb=transpose(linspace(Radiance(end,4*n-3:6*n-4),-xMidBorderNodes,(yNodeNum-n)/2+1));
yMidBorderNodes(2:2:2*n)=linspace((n-1)*yNodeSpace/2,-(n-1)*yNodeSpace/2,n);
yMidBorderNodes(1:2:2*n)=(w/2)*ones(1,n);
RightLimb=transpose(linspace(Radiance(end,2*n-1:4*n-2),yMidBorderNodes,(xNodeNum-n)/2+1));
LeftLimb=transpose(linspace(Radiance(end,[6*n-5:8*n-6]),-yMidBorderNodes,(xNodeNum-n)/2+1));

%Placing Corner Parts
yCornerBorderNodes(1:(yNodeNum-length(yMidBorderNodes)/2)/2+1,2)=transpose(linspace(yMidBorderNodes(end),-a/2,(yNodeNum-length(yMidBorderNodes)/2)/2+1));
yCornerBorderNodes(1:(yNodeNum-length(yMidBorderNodes)/2)/2+1,1)=-w/2*ones((yNodeNum-length(yMidBorderNodes)/2)/2+1,1);
xCornerCords=linspace(yCornerBorderNodes(:,1),BotLimb(:,end-1),(xNodeNum-length(xMidBorderNodes)/2)/2+1);
yCornerCords=linspace(yCornerBorderNodes(:,2),BotLimb(:,end),(xNodeNum-length(xMidBorderNodes)/2)/2+1);

TopLeftCorner(:,1:2:2*columns(xCornerCords))=xCornerCords;TopRightCorner(:,2:2:2*columns(xCornerCords))=-yCornerCords;
BotLeftCorner(:,1:2:2*columns(xCornerCords))=xCornerCords;BotRightCorner(:,2:2:2*columns(xCornerCords))=yCornerCords;
TopRightCorner(:,1:2:2*columns(xCornerCords))=-xCornerCords;TopLeftCorner(:,2:2:2*columns(xCornerCords))=-yCornerCords;
BotRightCorner(:,1:2:2*columns(xCornerCords))=-xCornerCords;BotLeftCorner(:,2:2:2*columns(xCornerCords))=yCornerCords;
%These are flipped, to not end up having negative elements with the loop below
TopRightCorner=flipud(TopRightCorner);
BotLeftCorner=flipud(BotLeftCorner);

%Packing all generated nodes into a n*2 matrix.
Nodes=[NodePacker(Radiance);NodePacker([BotLeftCorner BotLimb BotRightCorner]);NodePacker([RightLimb LeftLimb]);NodePacker([TopLeftCorner TopLimb TopRightCorner])];
%Creating Elements
Elements=[];
n=0;
Parts={Radiance BotLeftCorner BotLimb BotRightCorner RightLimb LeftLimb TopLeftCorner TopLimb TopRightCorner};
for k=1:length(Parts)
  part=Parts{k};R=rows(part);C=columns(part)/2;
  for j=1:C-1
    for i=1:R-1
      Elements(end+1,:)=[i+1+R*(j-1) i+R*(j-1) i+R*j i+1+R*j]+n;
    endfor
  endfor
  n=n+C*R;
endfor

%Merging Nodes and adjusting Elements to new nodes
[Nodes Waste Transform]=unique(round(Nodes*10^7),'rows');
Nodes=Nodes/10^7;
Elements=Transform(Elements);

%Adding Properties to Elements
Elements=[Elements ones(rows(Elements),1)*[s E v]];

%More Data
Split4Elements=[1:(rows(Radiance)-1)*(columns(Radiance)/2-1)];
LoadLines=[-1 -w/2;1 w/2];
ConstraintLine=0;

save Mesh6337.txt Elements Nodes LoadLines ConstraintLine Split4Elements

clear all
