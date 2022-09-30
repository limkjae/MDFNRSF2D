

function [ StiffnessMatrixShear,StiffnessMatrixNormal ]...
    = Function_StiffnessMatrix( FaultElementCount,FaultElementCenter,FaultLength,FaultAngle,ShearModulus,PoissonRatio, FaultRRLL)
%UNTITLED4 이 함수의 요약 설명 위치
%   자세한 설명 위치


for i=1:FaultElementCount
    for j=1:FaultElementCount
        
        FaultCenterThis=FaultElementCenter(j,:); %m [b,c]
        FaultHalfLengthThis=FaultLength(j)/2; %m
        FaultAngleThis=FaultAngle(j); %Degree from North
        
        ReceiverCenter=FaultElementCenter(i,:);
        ReceiverAngle=FaultAngle(i);

        if FaultAngleThis>0
            FaultShearDisp=-1;
        else
            FaultShearDisp=-1;
        end
        
        nFaultAngle=FaultAngleThis+90;
        lFaultAngle=FaultAngleThis;
        nFault=cos(deg2rad(nFaultAngle));
        lFault=cos(deg2rad(lFaultAngle));
        
        nReceiverAngle=ReceiverAngle+90;
        lReceiverAngle=ReceiverAngle;
        nReceiver=cos(deg2rad(nReceiverAngle));
        lReceiver=cos(deg2rad(lReceiverAngle));
        
        
        Xi=nFault*(ReceiverCenter(1)-FaultCenterThis(1))-lFault*(ReceiverCenter(2)-FaultCenterThis(2)); %twocurve X
        Zeta=lFault*(ReceiverCenter(1)-FaultCenterThis(1))+nFault*(ReceiverCenter(2)-FaultCenterThis(2));  %onecurve Z
        
        F3=((nFault^2-lFault^2)*Zeta-2*nFault*lFault*(Xi+FaultHalfLengthThis))/((Xi+FaultHalfLengthThis)^2+Zeta^2) ...
            - ((nFault^2-lFault^2)*Zeta-2*nFault*lFault*(Xi-FaultHalfLengthThis))/((Xi-FaultHalfLengthThis)^2+Zeta^2);
        F4=-(2*nFault*lFault*Zeta+(nFault^2-lFault^2)*(Xi+FaultHalfLengthThis))/((Xi+FaultHalfLengthThis)^2+Zeta^2) ...
            +(2*nFault*lFault*Zeta+(nFault^2-lFault^2)*(Xi-FaultHalfLengthThis))/((Xi-FaultHalfLengthThis)^2+Zeta^2);
        F5=(nFault*(nFault^2-3*lFault^2)*((Xi+FaultHalfLengthThis)^2-Zeta^2)+2*lFault*(3*nFault^2-lFault^2)*(Xi+FaultHalfLengthThis)*Zeta)/((Xi+FaultHalfLengthThis)^2+Zeta^2)^2 ...
            -(nFault*(nFault^2-3*lFault^2)*((Xi-FaultHalfLengthThis)^2-Zeta^2)+2*lFault*(3*nFault^2-lFault^2)*(Xi-FaultHalfLengthThis)*Zeta)/((Xi-FaultHalfLengthThis)^2+Zeta^2)^2;
        F6=(2*nFault*(nFault^2-3*lFault^2)*(Xi+FaultHalfLengthThis)*Zeta-lFault*(3*nFault^2-lFault^2)*((Xi+FaultHalfLengthThis)^2-Zeta^2))/((Xi+FaultHalfLengthThis)^2+Zeta^2)^2 ...
            -(2*nFault*(nFault^2-3*lFault^2)*(Xi-FaultHalfLengthThis)*Zeta-lFault*(3*nFault^2-lFault^2)*((Xi-FaultHalfLengthThis)^2-Zeta^2))/((Xi-FaultHalfLengthThis)^2+Zeta^2)^2;
        
        StressXX=FaultShearDisp*ShearModulus/2/pi/(1-PoissonRatio)*(2*nFault^2*F3 -2*nFault*lFault*F4 +Zeta*(nFault*F5 -lFault*F6));
        StressZZ=-FaultShearDisp*ShearModulus/2/pi/(1-PoissonRatio)*(2*lFault^2*F3 +2*nFault*lFault*F4 + Zeta*(nFault*F5 -lFault*F6));
        StressXZ=(FaultShearDisp*ShearModulus/2/pi/(1-PoissonRatio)*(F4+Zeta*(lFault*F5+nFault*F6)));
        if abs(StressXX)>1e10; StressXX=0; end
        if abs(StressZZ)>1e10; StressZZ=0; end
        if abs(StressXZ)>1e10; StressXZ=0; end
        
%         ResultSXX(i,j)=StressXX;
%         ResultSZZ(i,j)=StressZZ;
%         ResultSXZ(i,j)=StressXZ;
        
        StiffnessMatrixShear(i,j)=FaultRRLL(i)*FaultRRLL(j)*(nReceiver*lReceiver*(StressXX-StressZZ)+(nReceiver^2-lReceiver^2)*StressXZ);
        StiffnessMatrixNormal(i,j)=FaultRRLL(j)*(lReceiver^2*StressXX+nReceiver^2*StressZZ+2*nReceiver*lReceiver*StressXZ);
        
        
%         ResultReceiverShear(i,j)=ReceiverRRLL*(nReceiver*lReceiver*(StressXX-StressZZ)+(nReceiver^2-lReceiver^2)*StressXZ);        
%         ResultReceiverNormal(i,j)=-(lReceiver^2*StressXX+nReceiver^2*StressZZ+2*nReceiver*lReceiver*StressXZ);

        
        
    end
end

end

