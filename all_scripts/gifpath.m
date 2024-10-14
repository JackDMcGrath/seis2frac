function gifpath(gifflg,dir,timestep,pt,exhum,temp,P,gps,step,p,shp)
vel=dir./timestep*1e6;

t2mya=1e6/timestep;

x=[floor(min(pt(:,1))):ceil(max(pt(:,1)))];y=[floor(min(pt(:,2)))-5,ceil(max(pt(:,2)))+5];
[a,b]=meshgrid(x,y);

    gridspace=10;
    Xrange=([(floor(min(pt(:,1)-50)/50)*50+50):gridspace:ceil(max(pt(:,1)))]);
    if ceil(max(pt(:,2)))+5 < 0
    Yrange=([floor(min(pt(:,2)))-5:gridspace:ceil(max(pt(:,2)))+5]);
    else
        Yrange=([floor(min(pt(:,2)))-5:gridspace:0]);
    end
    if size(Yrange,2)==1
        Yrange(2)=Yrange(1)+gridspace;
    end
        if size(Xrange,2)==1
        Xrange(2)=Xrange(1)+gridspace;
    end
    [XX,YY]=meshgrid(Xrange,Yrange);
    ex=scatteredInterpolant(exhum(:,1),exhum(:,2),exhum(:,3));
    ex.Method='natural'; ex.ExtrapolationMethod='linear';
    ex=ex(XX(:),YY(:));
    ex(ex<0.1)=0.1;
    ex=reshape(ex,[length(Yrange),length(Xrange)]);
    load('C:/Jacko/scripts/ScientificColourMaps6/imola/imola.mat')
    
%     figure
% plot3(XX(:),YY(:),ex(:),'.')
% xlabel('X');ylabel('Y');zlabel('Z')
% hold on
% plot3(exhum(:,1),exhum(:,2),exhum(:,3),'r.','MarkerSize',20)

f=figure('Position',[2685 -161 1261 1096],'Visible','On');

if gifflg==1
    s=1;
else
    s=step;
end

for ii=s:step
    
    subplot(5,2,2);
    plot((1:ii)/t2mya,vel(1:ii,1),'.')
    title('Fault Parallel Velocity')
    xlim([0 step/t2mya])
    ylim([floor(min(vel(:,1))) ceil(max(vel(:,1)))])
    ylabel({'Parallel Vel','(mm/yr)'})
    grid
    
    ax2=subplot(5,2,4);
    plot((1:ii)/t2mya,vel(1:ii,2),'.')
    title('Fault Perpendicular Velocity')
    xlim([0 step/10])
    ylim([floor(min(vel(:,2))) ceil(max(vel(:,2)))])
    y2=ylabel({'Perp. Vel','(mm/yr)'});
    y2.Parent.YAxisLocation='right';
    grid
    
    ax3=subplot(5,2,6);
    plot((1:ii)/t2mya,-vel(1:ii,3),'.')
    title('Exhumation Rate')
    xlim([0 step/t2mya])
    ylim([floor(min(-vel(:,3))) ceil(max(-vel(:,3)))])
    xlabel('Time (Mya)')
    ylabel({'Exhumation Rate','(mm/yr)'});
    grid
    
    ax4=subplot(5,2,[1 3 5]);
    plot(P(1:ii)*1e-6,temp(1:ii),'b')
    grid('on')
    hold on
    plot(P(1:ii)*1e-6,temp(1:ii),'k.','MarkerSize',12)
    plot(P(5:5:ii)*1e-6,temp(5:5:ii),'kd','MarkerSize',8,'MarkerFaceColor','k')
    title('P-T-t')
    xlim([0 ceil(max(P)*1e-6)])
    ylim([0 ceil(max(temp)/50)*50])
    xlabel('Pressure (MPa)')
    ylabel({'Temperature (C)'});
    
    
    ax5=subplot(5,2,[7:10]);
    plot3(pt(1:ii,1),pt(1:ii,2),pt(1:ii,3),'b.')
    title('Exhumation Pathway and seismic volume')
    axis equal
    ax5.XLim=([Xrange(1) Xrange(end)]);
    ax5.YLim=([Yrange(1) Yrange(end)]);
    ax5.ZLim=([floor(min(pt(:,3)))-10 ceil(max(pt(:,3)))]);
    xlabel('Along Strike (km)')
    ylabel({'Fault Normal','(km)'},'Rotation',330,'HorizontalAlignment','center');
    zlabel('Depth (km)')
    grid('on')
    hold on
    view([330,20])
    g = hgtransform('Matrix',makehgtform('translate',[0 0 floor(min(pt(:,3)))-10]));
   imagesc(g,ax5.XLim,ax5.YLim,ex,[0 ceil(max(exhum(:,3)))]);
   c=colorbar; c.Label.String='Exhumation Rate (mm/yr)';
   colormap(imola)
    
    ax5.CameraPosition=[-78.0854 -492.6989 191.9934];
    ax5.CameraTarget=[195.0696 -19.5805 -6.8472];
    plot3(pt(1:ii,1),pt(1:ii,2),[(floor(min(pt(:,3)))-9)*ones(ii,1)],'r.')
    quiver3(Xrange(1)+15,Yrange(end),-10,30,0,0,0.2,'k')
    t=text(Xrange(1)+17,Yrange(end),-7,'30 mm/yr','Rotation',12.5,'HorizontalAlignment','center');
    gps_trim=1;
    quiver3(gps(1:gps_trim:end,1),gps(1:gps_trim:end,2),[(floor(min(pt(:,3)))-10)*ones(ceil(length(gps)/gps_trim),1)],...
        gps(1:gps_trim:end,3),gps(1:gps_trim:end,4),zeros(ceil(length(gps)/gps_trim),1),0.2,'k')
    plot3(pt(1:ii,1),pt(1:ii,2),pt(1:ii,3),'o','MarkerEdgeColor','w','MarkerFaceColor','b','MarkerSize',4)
    plot3(pt(1:ii,1),pt(1:ii,2),[(floor(min(pt(:,3)))-9)*ones(ii,1)],'o','MarkerEdgeColor','w','MarkerFaceColor','r','MarkerSize',4)
    plot3(pt(5:5:ii,1),pt(5:5:ii,2),pt(5:5:ii,3),'bd','MarkerFaceColor','b')
    plot3(pt(5:5:ii,1),pt(5:5:ii,2),[(floor(min(pt(:,3)))-9)*ones(length(pt(5:5:ii,1)),1)],'rd','MarkerFaceColor','r')
    legend('Rock Path','Surface Trace','GPS velocity','Location','NorthEast','AutoUpdate','off');
    if gifflg==1
    if ii==1
        gif(['newpath_',num2str(pt(1,1)),'_',num2str(pt(1,2)),'.gif'],'overwrite',true,'DelayTime',0.1,'Frame',f)
    else
        gif
    end
    end
    
end

for ii=s:step
    plot3(p([4*ii-3:4*ii 4*ii-3],1),p([4*ii-3:4*ii 4*ii-3],2),p([4*ii-3:4*ii 4*ii-3],3))
    if gifflg==1
    gif('DelayTime',0.05)
    end
end

%%
plot(shp);
    ax5.XLim=([Xrange(1) Xrange(end)]);
    ax5.YLim=([Yrange(1) Yrange(end)]);
    ax5.ZLim=([floor(min(pt(:,3)))-10 ceil(max(pt(:,3)))]);
if gifflg==1
    gif('DelayTime',5)

savefig(f,['newpath_',num2str(pt(1,1)),'_',num2str(pt(1,2)),'_final_plot.fig'])
end
