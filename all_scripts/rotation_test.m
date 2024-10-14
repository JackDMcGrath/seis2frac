dir=[1   1   1];

pt=[0 0 0];
poly_width=5; %width of poly in km
p=nan(4,3);
[p(1,:),p(2,:),p(3,:),p(4,:)]=find_corners(pt(1,1:3),poly_width);

 [m]=rotate3d(dir); % Make 3D rot-matrix from vector
    pnew=(m\p')';%+pt(ii,1:3); % Rotate original horizontal square, and translate to new center
    
    figure();  
    
    plot3(p(:,1),p(:,2),p(:,3),'r')
    hold on
    plot3(p(1,1),p(1,2),p(1,3),'k.')

     xlabel('X');ylabel('Y');zlabel('Z')
     plot3(pnew(:,1),pnew(:,2),pnew(:,3),'b')
     plot3(pnew(1,1),pnew(1,2),pnew(1,3),'k.')
     plot3([0 dir(1)],[0 dir(2)],[0 dir(3)])
    