% % polyxy=polybuffer(pt(:,1:2),'lines',5);figure();plot(polyxy);hold on; plot(pt(:,1),pt(:,2),'.');axis equal;title('xy')
% % polyxz=polybuffer(pt(:,[1 3]),'lines',5);figure();plot(polyxz);hold on; plot(pt(:,1),pt(:,3),'.');axis equal;title('xz')
% % polyyz=polybuffer(pt(:,2:3),'lines',5);figure();plot(polyyz);hold on; plot(pt(:,2),pt(:,3),'.');axis equal;title('yz')
% % 
% % [m]=AxelRot(45,[0 0 3],pt);
% % pnew=((m(1:3,1:3)*p')+m(1:3,4))';
% % plot3(pnew([1:4 1],1),pnew([1:4 1],2),pnew([1:4 1],3),'b')
% pt=[142.0609,-7.4839,0];
% dir=[-3 -3 -3];
% x0=[];
% deg=0;
% 
% figure
% 
% poly_width=5; %width of poly in km
% [p1,p2,p3,p4]=find_corners(pt(1,1:3),poly_width);
% p=[p1;p2;p3;p4];
% plot3(p([1:4 1],1),p([1:4 1],2),p([1:4 1],3),'r')
%  xlabel('X');ylabel('Y');zlabel('Z')
% hold on
% plot3(pt(1),pt(2),pt(3),'r*')
% [new,R,T]=AxelRot(p',deg,dir,x0);
% R
% T
% % new=new+dir';
% plot3(new(1,[1:4 1]),new(2,[1:4 1]),new(3,[1:4 1]),'b')
% % vector=[pt;pt+dir];
% % plot3(vector(:,1),vector(:,2),vector(:,3),'k')
% try
% vector2=[x0;x0+dir];
% catch
%     x0=[0 0 0];
%     vector2=[x1;x1+dir];
% end
% plot3(vector2(:,1),vector2(:,2),vector2(:,3),'k--')
% plot3(x0(1),x0(2),x0(3),'k*')
% axis equal


pt=[0,0,0];
dir=[1 0 1];
x0=[];
deg=0;
poly_width=5; %width of poly in km
[p1,p2,p3,p4]=find_corners(pt(1,1:3),poly_width);
p=[p1;p2;p3;p4]
p(1:2,2)=p(1:2,2)-1;
p(3:4,2)=p(3:4,2)+1;

figure()
plot3(p([1:4 1],1),p([1:4 1],2),p([1:4 1],3),'r')
xlabel('X');ylabel('Y');zlabel('Z')
axis equal
hold on
x=rotz(atand(dir(2)/dir(1)))*roty(atand(dir(3)/dir(1)))*rotx(atand(dir(3)/dir(2)));
prot=(x\(p-pt)')'+pt;
plot3(prot([1:4 1],1),prot([1:4 1],2),prot([1:4 1],3),'b')
xlabel('X');ylabel('Y');zlabel('Z')