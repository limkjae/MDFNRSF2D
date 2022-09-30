function [Return_V,Return_Friction,Return_Disp,Return_Theta,Return_EffNormalStress,Return_Dt,Return_InstabilityThistime]...
    = Function_OneStepSolve_LargeTimeStep_D(FaultIdx,ConvergenceCrit,DispOld,FrictionOld,ThetaOld,VOld,...
    Friction0,a,b,Dc,V0,K_Self,...
    Dt,Omega,Vl,AccelOld,...
    Mass,ShearModulus,EffNormalStress_Old,Instability,...
    Disp, Friction, V, Theta, DV,  Elastic_Load_Vel, EffNormalStress,Total_Loading_Disp)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
FaultIdx;

    Iteration=0;
    InstabilityThistime=0;
    maxDiff=10;
    ConvergenceCrit_Iter=ConvergenceCrit;

    while maxDiff>ConvergenceCrit_Iter
        Iteration=Iteration+1;
        VTest=V; % Velocity tested in this NR iteration
        
        VOldIter=V;
        
        % Original for NR
        Accel=(V-VOld)/Dt;
        %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
        Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
        Friction=FrictionOld + (K_Self/EffNormalStress*(Vl+Elastic_Load_Vel-V)...
            -(EffNormalStress-EffNormalStress_Old)/Dt*FrictionOld/EffNormalStress...
            -ShearModulus/6000/EffNormalStress*Accel...
            -Mass/EffNormalStress*(Accel-AccelOld)/Dt)*Dt;
        FV=V-V0*exp((Friction-Friction0-b*log(V0*Theta/Dc))/a);
        
        % Deviated value for NR
        V=V+DV;
        
        Accel=(V-VOld)/Dt;
        %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
        Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
        Friction=FrictionOld + (K_Self/EffNormalStress*(Vl+Elastic_Load_Vel-V)...
            -(EffNormalStress-EffNormalStress_Old)/Dt*FrictionOld/EffNormalStress...
            -ShearModulus/6000/EffNormalStress*Accel...
            -Mass/EffNormalStress*(Accel-AccelOld)/Dt)*Dt;
        FV_DeV=V-V0*exp((Friction-Friction0-b*log(V0*Theta/Dc))/a);
        
        V=VTest-FV*DV/(FV_DeV-FV);
        Disp=DispOld+V*Dt;
        Accel=(V-VOld)/Dt;
        
        
        
        
        if V<0
            if Instability==1
                V=1e-20;
                Disp=DispOld+(VOld+V)/2*Dt; % Update disp
                Theta=ThetaOld-V.*ThetaOld./Dc.*log(V.*ThetaOld./Dc)*Dt;
                Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
                InstabilityThistime=1;
                maxDiff=0;
            else
                
                if  Dt>1e-7
                    Dt=Dt/2;
                    V=VOld;
                    Theta=ThetaOld;
                    Friction=FrictionOld;
                    Disp=DispOld;
                    FaultIdx=0;
                    maxDiff=1;
                else
                    V=VOld;%10.^(log10(VOld)-(log10(VOld-Accel*DtOld)-log10(VOld))*Dt/DtOld);
                    Instability=1;
                    InstabilityThistime=3;
                    Disp=DispOld+(VOld+V)/2*Dt; % Update disp
                    Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
                    Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
                    maxDiff=0;
                end
            end
        elseif isreal(V)==0
            if Instability==1
                V=1e-20;
                Disp=DispOld+(VOld+V)/2*Dt; % Update disp
                %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
                Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
                Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
                InstabilityThistime=1;
                maxDiff=0;
            else
                if  Dt>1e-7
                    Dt=Dt/2;
                    V=VOld;
                    Theta=ThetaOld;
                    Friction=FrictionOld;
                    Disp=DispOld;
                    FaultIdx=0;
                    maxDiff=0;
                else
                    V=VOld;
                    Disp=DispOld+(VOld+V)/2*Dt; % Update disp
                    %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
                    Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
                    Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
                    Instability=1;
                    InstabilityThistime=1;
                    maxDiff=0;
                end
            end
            
        elseif isnan(V)==1
            if Instability==1
                V=1e-20;
                Disp=DispOld+(VOld+V)/2*Dt; % Update disp
                %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
                Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
                Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
                InstabilityThistime=1;
                maxDiff=0;
            else
                if  Dt>1e-7
                    Dt=Dt/2;
                    V=VOld;
                    Theta=ThetaOld;
                    Friction=FrictionOld;
                    Disp=DispOld;
                    FaultIdx=0;
                    maxDiff=0;
                else
                    V=VOld;
                    Disp=DispOld+(VOld+V)/2*Dt; % Update disp
                    %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
                    Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
                    Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
                    Instability=1;
                    InstabilityThistime=1;
                    maxDiff=0;
                end
            end
            
        else
            
            Disp=DispOld+(VOld+V)/2*Dt; % Update disp
            %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
            Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
            Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
            %             Instability=0;
            maxDiff=abs(VTest./V-1);
        end
        
        
        if rem(Iteration,500)==0;
            ConvergenceCrit_Iter=ConvergenceCrit_Iter*2; % Only used when convergence is failed
            %                 Dt=Dt/2;
        end
        %             maxDiff=max(abs(VTest./V-1));
        
        
    end
    Return_V=V;
    Return_Friction=Friction;
    Return_Disp=Disp;
    Return_Theta=Theta;
    Return_EffNormalStress=EffNormalStress;
    Return_Dt=Dt;
    Return_InstabilityThistime= InstabilityThistime;
    
%
% for kk=1:length(V)
%     if V(kk)<1e-20;
%         V(kk)=1e-20;
%         Disp(kk)=DispOld(kk);
%         Theta(kk)=ThetaOld(kk);
%         Friction(kk)=FrictionOld(kk);
%         Dt=Dt/5;
%     end
% end
