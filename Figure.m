% clear all
% load 500Fault2


clear EventCount EventTime EventLocation EventSlip EventMagnitude EventMomentMagnitude Tau Distance T R Eta

GRCumMagnitudeX=[0:0.05:6];
MangitudeXCount=length(GRCumMagnitudeX);
EventMomentCriteria=1e9
EventCount=0;
StepUntil=Step;
FaultMomentAccum=zeros(FaultCount,Step);
FaultMomentRate=zeros(FaultCount,Step);
CumEvent=zeros(Step,1);
for i=1:length(FaultNumberForElement)
    for j=1:StepUntil%Step
        FaultMomentAccum(FaultNumberForElement(i),j)=FaultMomentAccum(FaultNumberForElement(i),j)+History_Disp(j,i)*FaultElementLength(i)*500*30e9;
        if j>1
            FaultMomentRate(FaultNumberForElement(i),j)=(FaultMomentAccum(FaultNumberForElement(i),j)-FaultMomentAccum(FaultNumberForElement(i),j-1))/(History_Time(j)-History_Time(j-1));
        end
    end
end

sum(FaultMomentRate(:,2)==0)


EventCount=0;
for i=1:FaultCount
    
    EventCheck=0;
    EventCountSingleFault=0;
    for j=1:Step
        if FaultMomentRate(i,j)>EventMomentCriteria
            if EventCheck==0
                %                 if EventCountSingleFault>1
                %                     if History_Time(EventEnd)-History_Time(j)<1
                %                         EventContinue=1;
                EventBegin=j;
                EventCheck=1;
            end
        else
            if EventCheck==1
                EventCountSingleFault=EventCountSingleFault+1;
                EventCount=EventCount+1
                EventCheck=0;
                EventEnd=j;
                EventTime(EventCount)=History_Time(EventBegin);
                EventLocation(EventCount,:)=FaultCenter_Bulk(i,:);
                EventMagnitude(EventCount)=FaultMomentAccum(i,EventEnd)-FaultMomentAccum(i,EventBegin);
                EventMomentMagnitude(EventCount)=(log10(EventMagnitude(EventCount))-9.05)/1.5;
                %                 if EventMomentMagnitude(EventCount)<2;
                %                     EventCountSingleFault=EventCountSingleFault-1;
                %                     EventCount=EventCount-1;
                %                 end
            end
        end
    end
end


for i=1:Step
    %         CumEvent(i)=0;
    for j=1:EventCount
        if EventTime(j)<=History_Time(i)
            CumEvent(i)=CumEvent(i)+1;
        end
    end
end
InterpCount=10;

InterpInterval=round(History_Time(Step)/InterpCount);
InterpTime=[1:InterpInterval:History_Time(Step)];


CumMagnitudeCount=zeros(1,MangitudeXCount);
for i=1:EventCount
    for j=1:MangitudeXCount
        if EventMomentMagnitude(i) > GRCumMagnitudeX(j)
            CumMagnitudeCount(j)=CumMagnitudeCount(j)+1;
        end
    end
end



figure(2)
hold on
    scatter(EventTime/60/60/24,EventLocation(:,2)/1e3,6.^(EventMomentMagnitude-1),'bo')
set(gcf, 'color', 'w')
ylabel('Distance')
set(gca,'FontSize',13)
ylabel('Y distance (km)')
xlabel('time(days)')
% ylim([DomainMinY,DomainMaxY]/1e3)
box on


figure(1)
hold on
yyaxis left
set(gcf, 'color', 'w')
    scatter(EventTime/60/60/24,EventMomentMagnitude,6.^(EventMomentMagnitude-1),'b')
ylim([2,6])
ylabel('Magnitude')
set(gca,'FontSize',13)
set(gca,'ycolor','k')
yyaxis right
    plot(History_Time(1:Step)/60/60/24,CumEvent,'k','LineWidth',2)
ylabel('Cum. Count')
set(gca,'ycolor','k')
xlabel('time(year)')
box on
%





figure(4)
for i=1:FaultElementCount
    %         if FaultZoneOrNot_FaultElement==FaultTest;
    if History_Disp(Step,i)-History_Disp(1,i)>1e-3
            ColorIndex=[0,0,1];
            LineWidthIndex=2;
    else
        ColorIndex=[0.7,0.7,0.7];
        LineWidthIndex=1;
    end
    line([FaultX1(i),FaultX2(i)]/1000,[FaultY1(i),FaultY2(i)]/1000,'Color',ColorIndex,'LineWidth',LineWidthIndex)

end
set(gcf, 'color', 'w')
set(gca,'fontsize', 13)
xlabel('Distance (km)')
ylabel('Distance (km)')
box on
drawnow







%
% figure(1)
% hold on
% subplot(2,1,1)
% set(gcf, 'color', 'w')
% plot(History_Time/60/60/24/365,History_V)
% ylabel('Velocity (m/s)')
% set(gca, 'YScale', 'log')
% set(gca,'FontSize',13)
% % xlim([0,45])
% ylim([1e-11,1e0])
%
% subplot(2,1,2)
% set(gcf, 'color', 'w')
% scatter(EventTime/60/60/24/365,EventMomentMagnitude,'b')
% ylim([1.5,5])
% ylabel('Magnitude')
% set(gca,'FontSize',13)
% set(gca,'ycolor','b')
% yyaxis right
% plot(History_Time(1:Step)/60/60/24/365,CumEvent,'k','LineWidth',2)
% ylabel('Cum. Count')
% set(gca,'ycolor','k')
% xlabel('time(year)')
%
% figure(2)
% set(gcf, 'color', 'w')
% scatter(EventTime/60/60/24/365,EventMomentMagnitude,'b')
% ylim([1.5,2.5])
%
% ylabel('Magnitude')
% set(gca,'FontSize',13)
% set(gca,'ycolor','b')
% yyaxis right
% plot(History_Time(1:Step)/60/60/24/365,CumEvent,'k','LineWidth',2)
% ylabel('Cum. Count')
% set(gca,'ycolor','k')
% xlabel('time(year)')
% xlim([0,0.15])
%
%
% figure(7)
% hold on
% for i=1:FaultElementCount
%     line([FaultX1(i),FaultX2(i)],[FaultY1(i),FaultY2(i)],'LineWidth',2,'Color','k')
% end
% xlim([DomainMinX,DomainMaxX])
% ylim([DomainMinY,DomainMaxY])
% ylabel('Distance (m)')
% xlabel('Distance (m)')
% set(gcf, 'color', 'w')
% set(gca,'fontsize', 13)
% box on
%
%
%
% figure(4)
% hold on
% scatter(GRCumMagnitudeX,CumMagnitudeCount,'*')
% xlim([1.5,5])
% ylim([0.3,3e3])
% set(gca, 'YScale', 'log')
% set(gcf, 'color', 'w')
% set(gca, 'color', 'none');
% set(gca,'fontsize', 13)
% ylabel('Cum Count')
% xlabel('Magnitude')
% box on
% %
% figure(12)
% scatter(EventTime/60/60/24/365,EventLocation(:,2),'ko')
% set(gcf, 'color', 'w')
% ylabel('Distance')
% set(gca,'FontSize',13)
% ylabel('Distance from injection (m)')
% xlabel('time(year)')
% box on
%
%
% %
% figure(13)
% scatter(EventTime/60/60/24/365,EventLocation(:,2),4.^EventMomentMagnitude,'ko')
% set(gcf, 'color', 'w')
% ylabel('Distance')
% set(gca,'FontSize',13)
% ylabel('Distance from injection (m)')
% xlabel('time(year)')
% ylim([DomainMinY,DomainMaxY])
% box on
%
%
% figure(14)
% scatter(EventTime/60/60/24,EventLocation(:,2),6.^(EventMomentMagnitude-1),'ko')
% set(gcf, 'color', 'w')
% ylabel('Distance')
% set(gca,'FontSize',13)
% ylabel('Distance from injection (m)')
% xlabel('time(year)')
% ylim([DomainMinY,DomainMaxY])
% box on
%
%
% % figure(11)
% % plot(History_Time(1:Step)/60/60/24,CumEventSmooth,'b','LineWidth',2)
% % hold on
% % plot(History_Time(1:Step)/60/60/24,CumEvent,'r','LineWidth',2)
% % %
% % plot(History_Time,History_V)
% % set(gca, 'YScale', 'log')
% % set(gcf, 'color', 'w')
% % set(gca,'FontSize',13)
% % box on
% % xlim([0,800])
% % ylim([1e-11,1e0])
% %
% % figure(4)
% % scatter(EventTime,EventMomentMagnitude,'b')
% % xlim([0,800])
%
%
%
% %
% % RelationCount=0;
% % for i=1:EventCount
% %     for j=1:EventCount
% %         RelationCount=RelationCount+1;
% %         Tau(i,j)=EventTime(j)-EventTime(i);
% %         Distance(i,j)=sqrt((EventLocation(j,1)-EventLocation(i,1))^2+(EventLocation(j,2)-EventLocation(i,2))^2);
% %         if Tau(i,j)<0; Tau(i,j)=0; Distance(i,j)=0; end
% % %         if Tau(i,j)>Tau0; Tau(i,j)=0; Distance(i,j)=0; end
% % %         if Distance(i,j)>Dist0; Tau(i,j)=0; Distance(i,j)=0; end
% %         T(RelationCount)=Tau(i,j)*10^(-EventMomentMagnitude(i)/2);
% %         R(RelationCount)=Distance(i,j)^1.6*10^(-EventMomentMagnitude(i)/2);
% %         Eta(RelationCount)=T(RelationCount)*R(RelationCount);
% %     end
% % end
