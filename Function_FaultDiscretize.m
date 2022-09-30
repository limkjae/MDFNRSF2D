function [FaultElementCenter,FaultElementLength,FaultAngle,FaultAngleRad,FaultRRLL,InitialNormalStress,InitialShearStress,VlStress,FaultElementCount,SegmentCountPerFault,FaultLength_Overall, FaultNumberForElement,FaultLength_Bulk,FaultCenter_Bulk]...
    = Function_FaultDiscretize( FaultCount,FaultCenter_Bulk, FaultLength_Bulk, FaultAngle_Bulk,FaultRRLL_Bulk,InitialNormalStress_Bulk, InitialShearStress_Bulk, VlStress_Bulk,DiscretizeLength)
% Generate Random Faults We need 
% Apply NormalStresses

FaultAngleRad_Bulk=deg2rad(FaultAngle_Bulk);
% DiscretizeLength=500;

FaultX1_Bulk=FaultCenter_Bulk(:,1)+FaultLength_Bulk(:)/2.*sin(FaultAngleRad_Bulk(:));
FaultX2_Bulk=FaultCenter_Bulk(:,1)-FaultLength_Bulk(:)/2.*sin(FaultAngleRad_Bulk(:));
FaultY1_Bulk=FaultCenter_Bulk(:,2)+FaultLength_Bulk(:)/2.*cos(FaultAngleRad_Bulk(:));
FaultY2_Bulk=FaultCenter_Bulk(:,2)-FaultLength_Bulk(:)/2.*cos(FaultAngleRad_Bulk(:));

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
            FaultAngle(ElementIndex,1)=FaultAngle_Bulk(FaultIndex);
            FaultAngleRad(ElementIndex,1)=FaultAngleRad_Bulk(FaultIndex);
            FaultRRLL(ElementIndex,1)=FaultRRLL_Bulk(FaultIndex);
            InitialNormalStress(ElementIndex,1)=InitialNormalStress_Bulk(FaultIndex);
            InitialShearStress(ElementIndex,1)=InitialShearStress_Bulk(FaultIndex);
            VlStress(ElementIndex,1)=VlStress_Bulk(FaultIndex,1);
            FaultNumberForElement(ElementIndex,1)=FaultIndex;
        else
%             FaultLength(ElementIndex,1)=DiscretizeLength;
            FaultElementLength(ElementIndex,1)=FaultLength_Bulk(FaultIndex)/SegmentCountPerFault(FaultIndex);
            FaultAngle(ElementIndex,1)=FaultAngle_Bulk(FaultIndex);
            FaultAngleRad(ElementIndex,1)=FaultAngleRad_Bulk(FaultIndex);
            
            FaultX1(ElementIndex,1)=FaultX1_Bulk(FaultIndex)-(SegmentIndex-1)*FaultElementLength(ElementIndex).*sin(FaultAngleRad(ElementIndex));
            FaultY1(ElementIndex,1)=FaultY1_Bulk(FaultIndex)-(SegmentIndex-1)*FaultElementLength(ElementIndex).*cos(FaultAngleRad(ElementIndex));
            FaultX2(ElementIndex,1)=FaultX1_Bulk(FaultIndex)-SegmentIndex*FaultElementLength(ElementIndex).*sin(FaultAngleRad(ElementIndex));
            FaultY2(ElementIndex,1)=FaultY1_Bulk(FaultIndex)-SegmentIndex*FaultElementLength(ElementIndex).*cos(FaultAngleRad(ElementIndex));
            
            FaultElementCenter(ElementIndex,1)=(FaultX1(ElementIndex,1)+FaultX2(ElementIndex,1))/2;
            FaultElementCenter(ElementIndex,2)=(FaultY1(ElementIndex,1)+FaultY2(ElementIndex,1))/2;
            FaultLength_Overall(ElementIndex,1)=FaultLength_Bulk(FaultIndex);
            
            FaultRRLL(ElementIndex,1)=FaultRRLL_Bulk(FaultIndex);
            InitialNormalStress(ElementIndex,1)=InitialNormalStress_Bulk(FaultIndex);
            InitialShearStress(ElementIndex,1)=InitialShearStress_Bulk(FaultIndex);
            VlStress(ElementIndex,1)=VlStress_Bulk(FaultIndex);
            FaultNumberForElement(ElementIndex,1)=FaultIndex;
        
    end
        
    end
end

FaultElementCount=ElementIndex;



end

