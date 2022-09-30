function [Return_V,Return_Friction,Return_Disp,Return_Theta,Return_EffNormalStress,Return_Dt,Return_InstabilityThistime]...
    = Function_OneStepSolve_SmallTimeStep_D(FaultIdx,ConvergenceCrit,DispOld,FrictionOld,ThetaOld,VOld,...
    Friction0,a,b,Dc,V0,K_Self,...
    Dt,Omega,Vl,AccelOld,...
    Mass,ShearModulus,EffNormalStress_Old,Instability,...
    Disp, Friction, V, Theta, DV,  Elastic_Load_Vel, EffNormalStress,Total_Loading_Disp)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


    Iteration=0;
    InstabilityThistime=0;
    maxDiff=10;
    ConvergenceCrit_Iter=ConvergenceCrit;

    
    while maxDiff>ConvergenceCrit_Iter
        Iteration=Iteration+1;
        
        VOldIter=V;
        VTest=V; % Velocity tested in this NR iteration
        %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
        Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
        Friction=Friction0+a*log(V/V0)+b*log(Theta*V0./Dc)+ShearModulus/6000/EffNormalStress*V; % Initial friction
        F=Total_Loading_Disp-Friction*EffNormalStress/K_Self;
        Disp=(DispOld-F)*cos(Omega*Dt)+(VOld/Omega)*sin(Omega*Dt)+F; % Equation (9)
        V=(Disp-DispOld)/Dt*2-VOld;
        
        FOriginal=VOldIter-V; % We are testing this Newton Rhapson Function. Lets send this to zero
        
        % Deviated value for NR
        V=VOldIter+DV;
        VOldIter=V;
        
        %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
        Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
        Friction=Friction0+a*log(V/V0)+b*log(Theta*V0./Dc)+ShearModulus/6000/EffNormalStress*V; % Initial friction
        F=Total_Loading_Disp-Friction*EffNormalStress/K_Self;
        Disp=(DispOld-F)*cos(Omega*Dt)+(VOld/Omega)*sin(Omega*Dt)+F; % Equation (9)
        V=(Disp-DispOld)/Dt*2-VOld;
        
        NRF=VOldIter-V; % We are testing this Newton Rhapson Function. Lets send this to zero
        DF=(NRF-FOriginal)/DV; % tangent of the NR function
        V=VTest-FOriginal/DF; % Update velocity
        
        
        if V<0
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
                    maxDiff=1;
                else
                    %                 V=1e-10;
                    %                 V=VOld;
                    if V>1e-20
                        V=VOld/2;
                    else
                        V=1e-20;
                        Instability=1;
                        InstabilityThistime=1;
                    end
                    Disp=DispOld+(VOld+V)/2*Dt; % Update disp
                    %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
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
                    maxDiff=1;
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
            maxDiff=abs(VTest./V-1);
        end
        
        if rem(Iteration,100)==0;
            ConvergenceCrit_Iter=ConvergenceCrit_Iter*10; % Only used when convergence is failed
            
        end
        if ConvergenceCrit_Iter>1e-7
            InstabilityThistime=2;
        end
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
