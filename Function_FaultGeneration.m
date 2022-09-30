function [FaultElementCenter,FaultElementLength,FaultAngle,FaultAngleRad,FaultRRLL,InitialNormalStress,InitialShearStress,VlStress,FaultElementCount,SegmentCountPerFault,FaultLength_Overall, FaultNumberForElement,FaultLength_Bulk,FaultCenter_Bulk]...
    = Function_FaultGeneration( FaultCount,FaultAngle1,FaultAngle2,FaultLengthMin,FaultLengthMax,DomainMinX,DomainMaxX, DomainMinY,DomainMaxY ,Sigma1,Sigma3,LoadingRate )
% Generate Random Faults We need 
% Apply NormalStresses

DiscretizeLength=1000;
FaultLengthMin_BG=50;
FaultLengthMax_BG=3000;
MainFaultPortion=0.2;
FaultElementCount=FaultCount;
FaultCenterMinX_BG=DomainMinX;
FaultCenterMaxX_BG=DomainMaxX;
FaultCenterMinY_Fault=12000;
FaultCenterMaxY_Fault=28000;
Xdistribution=50;



RandNumber = slicesample(1,FaultCount,'pdf',@exppdf,'thin',5,'burnin',1000);
RandNumber=RandNumber/max(RandNumber);
FaultLength_Bulk=10.^(RandNumber*(log10(FaultLengthMax)-log10(FaultLengthMin))+log10(FaultLengthMin));
% % 
% %Fully Random
% FaultCenter_Bulk=[rand(FaultElementCount,1)*(FaultCenterMaxX-FaultCenterMinX)+FaultCenterMinX,rand(FaultElementCount,1)*(FaultCenterMaxY-FaultCenterMinY)+FaultCenterMinY];

% Random Directional
FaultCenter_RandNumber=rand(FaultElementCount,1);
[Center_Repulsive,L_Repulsive]=Function_RepulsiveDistribution(FaultCount)
FaultCenter_Y=Center_Repulsive(:,2)*(FaultCenterMaxY_Fault-FaultCenterMinY_Fault)/L_Repulsive+FaultCenterMinY_Fault;
FaultCenter_Bulk=[(FaultCenter_Y-(FaultCenterMaxY_Fault+FaultCenterMinY_Fault)/2)*tan(deg2rad(FaultAngle1))+(FaultCenterMaxY_Fault+FaultCenterMinY_Fault)/2+normrnd(0,Xdistribution,[FaultElementCount,1]),FaultCenter_Y];

% BackGround
% % repulsive

RandNumber_BG = slicesample(1,FaultCount,'pdf',@exppdf,'thin',5,'burnin',1000);
RandNumber_BG=RandNumber_BG/max(RandNumber_BG);
FaultLength_Bulk_BG=10.^(RandNumber_BG*(log10(FaultLengthMax_BG)-log10(FaultLengthMin_BG))+log10(FaultLengthMin_BG));

[Center_Repulsive_BG,L_Repulsive_BG]=Function_RepulsiveDistribution(FaultCount)
FaultCenter_Bulk_BG=Center_Repulsive_BG*(FaultCenterMaxX_BG-FaultCenterMinX_BG)/L_Repulsive_BG+FaultCenterMinX_BG;

FaultCenter_Bulk(FaultCount*MainFaultPortion+1:FaultCount,:)=FaultCenter_Bulk_BG(FaultCount*MainFaultPortion+1:FaultCount,:);
FaultLength_Bulk(FaultCount*MainFaultPortion+1:FaultCount)=FaultLength_Bulk_BG(FaultCount*MainFaultPortion+1:FaultCount);
% FaultLength_Bulk(1)=FaultLengthMax;

% FaultLength_Bulk(1)=100;
% FaultCenter_Bulk(1,:)=[DomainMaxX/1.5,DomainMaxY/2.5];
% FaultLength_Bulk(2)=200;
% FaultCenter_Bulk(2,:)=[DomainMaxX/3,DomainMaxY/3];
FaultSetRand_Bulk=2*rand(FaultCount,1)-1; %+1 for Positive angle and left latteral
% FaultSetRand_Bulk(1)=-1;
for i=1:FaultCount
    if FaultSetRand_Bulk(i)>0.8
        FaultAngle_Bulk(i,1)=FaultAngle1;
        FaultAngleRad_Bulk(i,1)=deg2rad(FaultAngle_Bulk(i,1));
        FaultSetRand_Bulk(i,1)=1;
        
    else
        FaultAngle_Bulk(i,1)=FaultAngle2;
        FaultAngleRad_Bulk(i,1)=deg2rad(FaultAngle_Bulk(i,1));
        FaultSetRand_Bulk(i,1)=-1;
    end
    if FaultAngle_Bulk(i,1)>0
        FaultRRLL_Bulk(i,1)=-1;
    else
        FaultRRLL_Bulk(i,1)=1;
    end
        InitialNormalStress_Bulk(i,1)=((Sigma1+Sigma3)/2+(Sigma1-Sigma3)*cos(2*(pi/2-FaultAngleRad_Bulk(i)))/2);
        InitialShearStress_Bulk(i,1)=(abs((Sigma1-Sigma3)*sin(2*(pi/2-FaultAngleRad_Bulk(i)))/2));
        FaultIndex_Bulk(i,1)=i;
end

FaultX1_Bulk=FaultCenter_Bulk(:,1)+FaultLength_Bulk(:)/2.*sin(FaultAngleRad_Bulk(:));
FaultX2_Bulk=FaultCenter_Bulk(:,1)-FaultLength_Bulk(:)/2.*sin(FaultAngleRad_Bulk(:));
FaultY1_Bulk=FaultCenter_Bulk(:,2)+FaultLength_Bulk(:)/2.*cos(FaultAngleRad_Bulk(:));
FaultY2_Bulk=FaultCenter_Bulk(:,2)-FaultLength_Bulk(:)/2.*cos(FaultAngleRad_Bulk(:));
Vl_Bulk=(LoadingRate.*cos(FaultAngleRad_Bulk));


for i=1:FaultCount %how many segments in each bulk?
    if FaultLength_Bulk(i) > DiscretizeLength % Define Target
%     if i==1 % Define Target
        SegmentCountPerFault(i)=ceil(FaultLength_Bulk(i)/DiscretizeLength);
    else 
        SegmentCountPerFault(i)=1;
    end
end

%Adjust for first large one
ElementIndex=0;
for FaultIndex=1:FaultCount %how many segments in each bulk?
    SegmentIndex=0;
    for j=1:SegmentCountPerFault(FaultIndex)
        ElementIndex=ElementIndex+1;
        SegmentIndex=SegmentIndex+1;
        if SegmentCountPerFault(FaultIndex)==1
            FaultElementCenter(ElementIndex,1)=FaultCenter_Bulk(FaultIndex,1);
            FaultElementCenter(ElementIndex,2)=FaultCenter_Bulk(FaultIndex,2);
            FaultElementLength(ElementIndex,1)=FaultLength_Bulk(FaultIndex);
            FaultLength_Overall(ElementIndex,1)=FaultLength_Bulk(FaultIndex);
            FaultSetRand(ElementIndex,1)=FaultSetRand_Bulk(FaultIndex);
            FaultAngle(ElementIndex,1)=FaultAngle_Bulk(FaultIndex);
            FaultAngleRad(ElementIndex,1)=FaultAngleRad_Bulk(FaultIndex);
            FaultX1(ElementIndex,1)=FaultX1_Bulk(FaultIndex);
            FaultX2(ElementIndex,1)=FaultX2_Bulk(FaultIndex);
            FaultY1(ElementIndex,1)=FaultY1_Bulk(FaultIndex);
            FaultY2(ElementIndex,1)=FaultY2_Bulk(FaultIndex);
            FaultRRLL(ElementIndex,1)=FaultRRLL_Bulk(FaultIndex);
            InitialNormalStress(ElementIndex,1)=InitialNormalStress_Bulk(FaultIndex);
            InitialShearStress(ElementIndex,1)=InitialShearStress_Bulk(FaultIndex);
            VlStress(ElementIndex)=Vl_Bulk(FaultIndex,1);
            FaultNumberForElement(ElementIndex,1)=FaultIndex;
        else
%             FaultLength(ElementIndex,1)=DiscretizeLength;
            FaultElementLength(ElementIndex,1)=FaultLength_Bulk(FaultIndex)/SegmentCountPerFault(FaultIndex);
            FaultAngle(ElementIndex,1)=FaultAngle_Bulk(FaultIndex);
            FaultAngleRad(ElementIndex,1)=FaultAngleRad_Bulk(FaultIndex);
%             if j==FaultSegmentCount(FaultIndex)
%             FaultX1(ElementIndex,1)=FaultX1_Bulk(FaultIndex)-(SegmentIndex-1)*FaultLength(ElementIndex).*sin(FaultAngleRad(ElementIndex));
%             FaultY1(ElementIndex,1)=FaultY1_Bulk(FaultIndex)-(SegmentIndex-1)*FaultLength(ElementIndex).*cos(FaultAngleRad(ElementIndex));
%             FaultX2(ElementIndex,1)=FaultX2_Bulk(FaultIndex);
%             FaultY2(ElementIndex,1)=FaultY2_Bulk(FaultIndex);
%             else
            FaultX1(ElementIndex,1)=FaultX1_Bulk(FaultIndex)-(SegmentIndex-1)*FaultElementLength(ElementIndex).*sin(FaultAngleRad(ElementIndex));
            FaultY1(ElementIndex,1)=FaultY1_Bulk(FaultIndex)-(SegmentIndex-1)*FaultElementLength(ElementIndex).*cos(FaultAngleRad(ElementIndex));
            FaultX2(ElementIndex,1)=FaultX1_Bulk(FaultIndex)-SegmentIndex*FaultElementLength(ElementIndex).*sin(FaultAngleRad(ElementIndex));
            FaultY2(ElementIndex,1)=FaultY1_Bulk(FaultIndex)-SegmentIndex*FaultElementLength(ElementIndex).*cos(FaultAngleRad(ElementIndex));
%             end
            FaultElementCenter(ElementIndex,1)=(FaultX1(ElementIndex,1)+FaultX2(ElementIndex,1))/2;
            FaultElementCenter(ElementIndex,2)=(FaultY1(ElementIndex,1)+FaultY2(ElementIndex,1))/2;
            FaultLength_Overall(ElementIndex,1)=FaultLength_Bulk(FaultIndex);
            
            FaultSetRand(ElementIndex,1)=FaultSetRand_Bulk(FaultIndex);
            FaultRRLL(ElementIndex,1)=FaultRRLL_Bulk(FaultIndex);
            InitialNormalStress(ElementIndex,1)=InitialNormalStress_Bulk(FaultIndex);
            InitialShearStress(ElementIndex,1)=InitialShearStress_Bulk(FaultIndex);
            VlStress(ElementIndex,1)=Vl_Bulk(FaultIndex);
            FaultNumberForElement(ElementIndex,1)=FaultIndex;
        
    end
        
    end
end
% 
% FaultCenter=FaultCenter_Bulk
% FaultLength=FaultLength_Bulk
% FaultSetRand=FaultSetRand_Bulk
% FaultAngle=FaultAngle_Bulk
% FaultAngleRad=FaultAngleRad_Bulk
% FaultX1=FaultX1_Bulk
% FaultX2=FaultX2_Bulk
% FaultY1=FaultY1_Bulk
% FaultY2=FaultY2_Bulk
% FaultRRLL=FaultRRLL_Bulk
% InitialNormalStress=InitialNormalStress_Bulk
% InitialShearStress=InitialShearStress_Bulk,
% Vl=Vl_Bulk
FaultElementCount=ElementIndex;



end

