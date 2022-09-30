%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MDFNRSF2D Matlab version
%%% by Kyungjae (KJ) Im. 
%%% limkjae10@gmail.com, kjim@caltech.edu
%%%
%%% For more information -
%%% Cascading Foreshocks and Aftershocks in a Discrete Fault Network
%%% Kyungjae Im and Jean-Philippe Avouac
%%% Geophysical Journal International
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all


% % % Input parameters
ShearModulus=20e9;
PoissonRatio=0.2;
RockDensity=2400; %kg/m3

TotalStep=200000;
RecordStep=20; 
TotalRecord=TotalStep/RecordStep;
Dt=1
ConvergenceCrit=1e-9;
SwitchV=1e-6;
RSF_a=0.003;
RSF_b=0.006;
RSF_Dc=100e-6
LoadingRate=0;


%FaultGeneration
FaultCount=4;
FaultCenter_Bulk=...
    [6000,6000;...
    2000,5000;...
    10000,2000;...
    8000,10000];
FaultLength_Bulk=[8000;1000;1000;1000];
FaultAngle_Bulk=[10;30;30;30];
FaultRRLL_Bulk=[-1;-1;-1;-1]; %+1 is right lateral
InitialNormalStress_Bulk=[1;1;1;1]*15e6;
InitialShearStress_Bulk=InitialNormalStress_Bulk*0.6;%[8660254.03784439;8660254.03784439;8660254.03784439;8660254.03784439]
VlStress_Bulk=[1;1;1;1]*LoadingRate; % Stressing Rate

% Discretize
DiscretizeLength=500;
[FaultElementCenter,FaultElementLength,FaultAngle,FaultAngleRad,FaultRRLL,InitialNormalStress,...
    InitialShearStress,VlStress,FaultElementCount,SegmentCountPerFault,FaultLength_Overall,...
    FaultNumberForElement,FaultLength_Bulk,FaultCenter_Bulk]...
    = Function_FaultDiscretize( FaultCount,FaultCenter_Bulk, FaultLength_Bulk, FaultAngle_Bulk,...
    FaultRRLL_Bulk,InitialNormalStress_Bulk, InitialShearStress_Bulk, VlStress_Bulk,DiscretizeLength);


% Generate Stiffness Matrix
[StiffnessMatrixShear, StiffnessMatrixNormal]...
    = Function_StiffnessMatrix( FaultElementCount,FaultElementCenter,FaultElementLength,FaultAngle,ShearModulus,PoissonRatio, FaultRRLL);

% Remove faults with too strong interaction
[FaultElementCount, FaultElementCenter, FaultElementLength,...
    FaultAngle, FaultAngleRad, FaultRRLL, InitialShearStress,...
    InitialNormalStress, FaultLength_Overall, FaultNumberForElement]...
    = Function_Check_TooClose(FaultElementCenter,FaultElementLength,FaultAngle,FaultAngleRad,...
    FaultRRLL,InitialNormalStress,InitialShearStress,VlStress,FaultElementCount,SegmentCountPerFault,...
    FaultLength_Overall, FaultNumberForElement,FaultLength_Bulk,FaultCenter_Bulk,...
    StiffnessMatrixShear, StiffnessMatrixNormal);

% Recalculate the Stiffness
[StiffnessMatrixShear, StiffnessMatrixNormal] = ...
    Function_StiffnessMatrix( FaultElementCount,FaultElementCenter,FaultElementLength,FaultAngle,ShearModulus,PoissonRatio, FaultRRLL);

% Initial Values
TMax=1e1;
TMin=2.0000e+9; %1e9
TRand=rand(FaultCount,1)*(TMax-TMin)+TMin;
H=0.006/100e-6-1e8/15e6;
Vini_Fault=0.003./H.*(1./TRand);
for i=1:FaultElementCount
    Vini(i,1)=Vini_Fault(FaultNumberForElement(i));
end
ThetaI=ones(FaultElementCount,1)*1e10;



V0=1e-9;
a=RSF_a*ones(FaultElementCount,1);
b=RSF_b*ones(FaultElementCount,1);
Dc=RSF_Dc*ones(FaultElementCount,1);


ElasticLoadingShear=0;
ElasticLoadingNormal=0;
Omega=0;
Mass=0;
K_Self=0;
Mass=RockDensity*FaultLength_Overall/(1-PoissonRatio)/pi/pi;
for i=1:FaultElementCount
    Omega(i,1)=sqrt(StiffnessMatrixShear(i,i)/Mass(i));
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






Far_Load_Disp_Initial=(StiffnessMatrixShear\InitialShearStress); % initial load point

FrictionI=InitialShearStress./InitialNormalStress;
Friction0=FrictionI-a.*log(Vini./V0)-b.*log(ThetaI.*V0./Dc); % Initial friction


% Fault Plot

FaultX1=FaultElementCenter(:,1)+FaultElementLength(:)/2.*sin(FaultAngleRad(:));
FaultX2=FaultElementCenter(:,1)-FaultElementLength(:)/2.*sin(FaultAngleRad(:));
FaultY1=FaultElementCenter(:,2)+FaultElementLength(:)/2.*cos(FaultAngleRad(:));
FaultY2=FaultElementCenter(:,2)-FaultElementLength(:)/2.*cos(FaultAngleRad(:));
% VlStress=(LoadingRate.*cos(FaultAngleRad));
Vl=VlStress./K_Self;

figure(1)
hold on
for i=1:FaultElementCount
    line([FaultX1(i),FaultX2(i)],[FaultY1(i),FaultY2(i)],'Color','red')
end
% xlim([DomainMinY,DomainMaxY])
% ylim([DomainMinY,DomainMaxY])
drawnow



% Set zeros
History_Time=zeros(TotalRecord,1);
History_Disp=zeros(TotalRecord,FaultElementCount);
History_V=zeros(TotalRecord,FaultElementCount);
History_Friction=zeros(TotalRecord,FaultElementCount);
History_Theta=zeros(TotalRecord,FaultElementCount);
History_NormShearStress=zeros(TotalRecord,FaultElementCount);
History_NormalStress=zeros(TotalRecord,FaultElementCount);
History_Accel=zeros(TotalRecord,FaultElementCount);
History_Dt=zeros(TotalRecord,1);
Accel=0;
Far_Load_Disp_Old=Far_Load_Disp_Initial;%Disp_Initial;
TOld=0;
Step=0;
Vmax=Vl;
T=0;

Friction=FrictionI;
DispOld=zeros(FaultElementCount,1);%Disp_Initial;
Disp=zeros(FaultElementCount,1);%Disp_Initial;
Instability=zeros(FaultElementCount,1);%Disp_Initial;
AccelOld=zeros(FaultElementCount,1);
SolverSwitch=zeros(FaultElementCount,1);
Theta=ThetaI;
ThetaOld=Theta;
VOld=Vini;
FrictionOld=Friction;
VmaxOld=0;
EffNormalStress_Old=InitialNormalStress;

SlowOrFast=0;
for i=1:TotalStep
    % Dt Control
    DtOld=Dt;
    VDiff=10; % Arbitrary for initiation
    Iteration=0; % Number of iteration
    V=VOld;
    
    [DtRef] = Function_Dt_Reference(V, Accel,LoadingRate,SlowOrFast);
    %     Dt=DtRef;
    if Dt<DtRef/1.5
        Dt=Dt*1.5;
    else
        Dt=DtRef;
    end
    
    if max(Instability)==2
        Dt=Dt*100;
    end
    
    
    FaultIdx=0;         
    Terminate=0;
    while Terminate==0
        
        DtOld=Dt;
        Disp_i=DispOld;
        Friction_i=FrictionOld;
        V_i=VOld;
        Theta_i=zeros(FaultElementCount,1);
        Dt_All=ones(FaultElementCount,1)*Dt;
        DV=10.^(log10(VOld)-(log10(VOld-AccelOld*DtOld)-log10(VOld))*Dt/DtOld);
        
        EffNormalStress_i=InitialNormalStress+StiffnessMatrixNormal*DispOld;
        for iii=1:FaultElementCount
            if EffNormalStress_i(iii)<1e6
                EffNormalStress_i(iii)=1e6;
            end
        end
        
        Far_Load_Disp=Far_Load_Disp_Old+Vl.*Dt;
        Elastic_Load_Disp=ElasticLoadingShear*(Far_Load_Disp-DispOld);
        Elastic_Load_Vel=ElasticLoadingShear*(Vl-VOld);
        Total_Loading_Disp=Far_Load_Disp+Elastic_Load_Disp;
        
        
        %%%
%         parfor FaultIdx=1:FaultCount
        for FaultIdx=1:FaultElementCount
            
            if SolverSwitch(FaultIdx)==1
                 [V(FaultIdx),Friction(FaultIdx),Disp(FaultIdx),Theta(FaultIdx),EffNormalStress(FaultIdx),Dt_All(FaultIdx),InstabilityThistime(FaultIdx)]...
                    = Function_OneStepSolve_SmallTimeStep_D(FaultIdx,ConvergenceCrit,DispOld(FaultIdx),FrictionOld(FaultIdx),ThetaOld(FaultIdx),VOld(FaultIdx),...
                    Friction0(FaultIdx),a(FaultIdx),b(FaultIdx),Dc(FaultIdx),V0,K_Self(FaultIdx),...
                    Dt,Omega(FaultIdx),Vl(FaultIdx),AccelOld(FaultIdx),...
                    Mass(FaultIdx),ShearModulus, EffNormalStress_Old(FaultIdx),Instability(FaultIdx),...
                    Disp_i(FaultIdx), Friction_i(FaultIdx), V_i(FaultIdx), Theta_i(FaultIdx), DV(FaultIdx),  Elastic_Load_Vel(FaultIdx), EffNormalStress_i(FaultIdx),Total_Loading_Disp(FaultIdx));
            else 
                [V(FaultIdx),Friction(FaultIdx),Disp(FaultIdx),Theta(FaultIdx),EffNormalStress(FaultIdx),Dt_All(FaultIdx),InstabilityThistime(FaultIdx)]...
                    = Function_OneStepSolve_LargeTimeStep_D(FaultIdx,ConvergenceCrit,DispOld(FaultIdx),FrictionOld(FaultIdx),ThetaOld(FaultIdx),VOld(FaultIdx),...
                    Friction0(FaultIdx),a(FaultIdx),b(FaultIdx),Dc(FaultIdx),V0,K_Self(FaultIdx),...
                    Dt,Omega,Vl(FaultIdx),AccelOld(FaultIdx),...
                    Mass(FaultIdx),ShearModulus, EffNormalStress_Old(FaultIdx),Instability(FaultIdx),...
                    Disp_i(FaultIdx), Friction_i(FaultIdx), V_i(FaultIdx), Theta_i(FaultIdx), DV(FaultIdx),  Elastic_Load_Vel(FaultIdx), EffNormalStress_i(FaultIdx),Total_Loading_Disp(FaultIdx));
            end
            if InstabilityThistime(FaultIdx)>0
                SolverSwitch(FaultIdx)=1;
            end
        end
        %
        if min(Dt_All)<max(Dt_All)
            if max(InstabilityThistime)==3
                Dt=1e-5;
            else
            Dt= min(Dt_All);
            end
        else
            Terminate=1;
        end
        
        
    end
    Dt=min(Dt_All);

    T=TOld+Dt;
    
    if rem(i,RecordStep)==0
        Step=Step+1;
        [i/TotalStep,log10(max(V)),log10(Dt), SlowOrFast] 
        
        History_Time(Step,1)=T;
        History_Dt(Step,1)=Dt;
        History_Disp(Step,:)=Disp;
        History_V(Step,:)=V;
        History_Friction(Step,:)=Friction;
        History_Theta(Step,:)=Theta;
        History_NormalStress(Step,:)=EffNormalStress;
%         History_Pressure(Step,:)=Pressure;
        History_Accel(Step,:)=Accel;
        History_Far_Load_Disp(Step,:)=Far_Load_Disp;
    end
    Vmax=max(V);

    
    VmaxOld=max(V);
    Accel=(V-VOld)/Dt;
    AccelOld=Accel;
    DispOld=Disp;
    ThetaOld=Theta;
    VOld=V;
    FrictionOld=Friction;
    Far_Load_Disp_Old= Far_Load_Disp;
    EffNormalStress_Old=EffNormalStress;
    TOld=T;
    Instability=InstabilityThistime;
    
    for idx=1:FaultElementCount
        if V(idx)>SwitchV
            SolverSwitch(idx)=1;
        elseif abs(Accel(idx))>0.001
            SolverSwitch(idx)=1;
        else
            SolverSwitch(idx)=0;
        end
    end
    
    if max(SolverSwitch)>0
        SlowOrFast=1;
    else
        SlowOrFast=0;
    end
    
%     if Vmax>SwitchV
%         SlowOrFast=1;
%     else
%         SlowOrFast=0;
%     end
end

save Result

% figure(2)
% surf(StiffnessMatrixShear)

figure(2)
hold on
plot(History_Time,log10(real(History_V(:,:))))
% 
% figure(3)
% plot(History_Time,History_Friction(:,:))
% 
% figure(4)
% hold on
% plot(History_Time,History_NormalStress/1e6)
% 
% figure(10)
% hold on
% plot(History_Time,(History_NormalStress.*History_Friction)/1e6)
% 
% figure(5)
% plot(History_Time,(History_Disp(:,:)))
% 
% figure(6)
% hold on
% plot(log10(real(History_V(:,:))))
% 
% figure(7)
% hold on
% plot(History_Time,History_Dt)
% 
% figure(8)
% plot(History_Time,log10(History_Dt))
% 
% figure(9)
% plot(History_Time,History_Theta(:,:))