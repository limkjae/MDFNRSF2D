function [FaultElementCount_Return, FaultElementCenter_Return, FaultElementLength_Return,...
    FaultAngle_Return, FaultAngleRad_Return, FaultRRLL_Return, InitialShearStress_Return,...
    InitialNormalStress_Return, FaultLength_Overall_Return, FaultNumberForElement_Return]...
    = Function_Check_TooClose(FaultElementCenter,FaultElementLength,FaultAngle,FaultAngleRad,...
    FaultRRLL,InitialNormalStress,InitialShearStress,VlStress,FaultElementCount,SegmentCountPerFault,...
    FaultLength_Overall, FaultNumberForElement,FaultLength_Bulk,FaultCenter_Bulk,...
    StiffnessMatrixShear, StiffnessMatrixNormal)
% This function removes the faults with too strong interaction 

TooClose=zeros(FaultElementCount,1);
UnstableCount=0;
for i=1:FaultElementCount
    K_Self(i,1)=StiffnessMatrixShear(i,i);
    for j=1:FaultElementCount
        if i==j
            ElasticLoadingShear(i,j)=0;
        else
            ElasticLoadingShear(i,j)=StiffnessMatrixShear(i,j)/StiffnessMatrixShear(i,i);
            ElasticLoadingNormal(i,j)=StiffnessMatrixNormal(i,j)/StiffnessMatrixShear(i,i);
            
        end
    end
end

for i=1:FaultElementCount
    for j=1:FaultElementCount
        if abs(ElasticLoadingShear(i,j))>0.8
            if TooClose(i)==0 & TooClose(j)==0
                if FaultLength_Overall(i)<FaultLength_Overall(j)
                    UnstableCount=UnstableCount+1
                    TooClose(i)=1;   
                else
                    UnstableCount=UnstableCount+1
                    TooClose(j)=1;
                end
            end
        end
        if abs(ElasticLoadingNormal(i,j))>0.8
            if TooClose(i)==0 & TooClose(j)==0
                if FaultLength_Overall(i)<FaultLength_Overall(j)
                    UnstableCount=UnstableCount+1
                    TooClose(i)=1;
                else
                    UnstableCount=UnstableCount+1
                    TooClose(j)=1;
                end
            end
        end
    end
end

FaultElementIdx=0;
for i=1:length(TooClose)
    if TooClose(i)==0
        FaultElementIdx=FaultElementIdx+1;
        FaultElementCenter_Reduced(FaultElementIdx,:)=FaultElementCenter(i,:);
        FaultElementLength_Reduced(FaultElementIdx,1)=FaultElementLength(i);
        FaultAngle_Reduced(FaultElementIdx,1)=FaultAngle(i);
        FaultAngleRad_Reduced(FaultElementIdx,1)=FaultAngleRad(i);
        FaultRRLL_Reduced(FaultElementIdx,1)=FaultRRLL(i);
        InitialShearStress_Reduced(FaultElementIdx,1)=InitialShearStress(i);
        InitialNormalStress_Reduced(FaultElementIdx,1)=InitialNormalStress(i);
        FaultLength_Overall_Reduced(FaultElementIdx,1)=FaultLength_Overall(i);
        FaultNumberForElement_Reduced(FaultElementIdx,1)=FaultNumberForElement(i);
    end
end
FaultElementCount_Return=FaultElementIdx;
FaultElementCenter_Return=FaultElementCenter_Reduced;
FaultElementLength_Return=FaultElementLength_Reduced;
FaultAngle_Return=FaultAngle_Reduced;
FaultAngleRad_Return=FaultAngleRad_Reduced;
FaultRRLL_Return=FaultRRLL_Reduced;
InitialShearStress_Return=InitialShearStress_Reduced;
InitialNormalStress_Return=InitialNormalStress_Reduced;
FaultLength_Overall_Return=FaultLength_Overall_Reduced;
FaultNumberForElement_Return=FaultNumberForElement_Reduced;

% Until Here



end

