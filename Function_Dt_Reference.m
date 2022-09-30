function [ DtRef ] = DT_Reference(V, Accel, LoadingRate,SlowOrFast)
%UNTITLED2 이 함수의 요약 설명 위치
%   자세한 설명
if SlowOrFast==1;
%     DtRef=1e-6;
    
    Dt_Lower=5e-2; % Time step [second]
    Dt_Upper=1e-2; % Time step [second]
    VTreshLower=1e-4;
    VTreshUpper=1e-2;
    
    Vmax=max(V);
    
    if Vmax<VTreshLower
        DtRef=Dt_Lower;
    elseif Vmax > VTreshLower & Vmax < VTreshUpper
        x=(log10(Vmax)-log10(VTreshLower))/(log10(VTreshUpper)-log10(VTreshLower));
        y=3*x^2-2*x^3;
        LogDt=log10(Dt_Lower)+y*(log10(Dt_Upper)-log10(Dt_Lower));
        DtRef=10^LogDt;
        %     DtRef=Dt_Lower-(Dt_Lower-Dt_Upper)*(3*((Vmax-VTreshLower)/(VTreshUpper-VTreshLower))^2-2*((Vmax-VTreshLower)/(VTreshUpper-VTreshLower))^3);
    else
        DtRef=Dt_Upper;
    end
    
    
    
    
else
    
    
    Dt_Lower=10000; % Time step [second]
    Dt_Upper=5e-2; % Time step [second]
    VTreshLower=1e-7;
    VTreshUpper=1e-5;
    
    Vmax=max(V);
    %
    if Vmax<VTreshLower
        DtRef=Dt_Lower;
    elseif Vmax > VTreshLower & Vmax < VTreshUpper
        x=(log10(Vmax)-log10(VTreshLower))/(log10(VTreshUpper)-log10(VTreshLower));
        y=3*x^2-2*x^3;
        LogDt=log10(Dt_Lower)+y*(log10(Dt_Upper)-log10(Dt_Lower));
        DtRef=10^LogDt;
        %     DtRef=Dt_Lower-(Dt_Lower-Dt_Upper)*(3*((Vmax-VTreshLower)/(VTreshUpper-VTreshLower))^2-2*((Vmax-VTreshLower)/(VTreshUpper-VTreshLower))^3);
    else
        DtRef=Dt_Upper;
    end
    
    %
    % if Vmax<VSlowTreshSlow
    %     DtRef=Dt_Slow_Slow;
    % elseif Vmax > VSlowTreshSlow & Vmax < VSlowTreshFast
    %     DtRef=Dt_Slow_Slow-(Dt_Slow_Slow-Dt_Slow_Fast)*(3*((Vmax-VSlowTreshSlow)/(VSlowTreshFast-VSlowTreshSlow))^2-2*((Vmax-VSlowTreshSlow)/(VSlowTreshFast-VSlowTreshSlow))^3);
    % elseif Vmax > VSlowTreshFast & Vmax < VSlowTreshAfterPeak
    %     DtRef=Dt_Slow_Fast-(Dt_Slow_Fast-Dt_Slow_AfterPeak)*(3*((Vmax-VSlowTreshFast)/(VSlowTreshAfterPeak-VSlowTreshFast))^2-2*((Vmax-VSlowTreshFast)/(VSlowTreshAfterPeak-VSlowTreshFast))^3);
    % else
    %     DtRef=Dt_Slow_AfterPeak;
    % end
    
    
    
    
    
end


end
