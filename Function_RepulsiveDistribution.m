function [X,L]...
    = Function_RepulsiveDistribution(FaultElementCount)


L = 50; % <-- Choose length of square sides
x0 = L/2; y0 = L/2; % <-- Choose center of square
n = FaultElementCount; % <-- Choose number of points
% Generate 2D epanechnikov-distribution
X = [x0,y0] + (sin(asin(2*rand(n,2)-1)/3))*L;
XYR = [x0,y0]+[[-1;1;1;-1;-1],[-1;-1;1;1;-1]]*L/2;
XB = interp1((0:4)'*L,XYR,linspace(0,4*L,200));
XB(end,:) = [];
nrepulsion = 100;
% Repulsion of seeds to avoid them to be too close to each other
n = size(X,1);
Xmin = [x0-L/2,y0-L/2];
Xmax = [x0+L/2,y0+L/2];
% Point on boundary
XR = x0+[-1,1,1,-1,-1]*L/2;
YR = y0+[-1,-1,1,1,-1]*L/2;
cla;
hold on
% figure(20)
% plot(XR,YR,'r-');
% h = plot(X(:,1),X(:,2),'b.');
axis equal
dmin = 0.5;
d2min = dmin*dmin;
beta = 0.5;
for k = 1:nrepulsion
    XALL = [X; XB];
    DT = delaunayTriangulation(XALL);
    T = DT.ConnectivityList;
    containX = ismember(T,1:n);
    b = any(containX,2);
    TX = T(b,:);
    [r,i0] = find(containX(b,:));
    i = mod(i0+(-1:1),3)+1;
    i = TX(r + (i-1)*size(TX,1));
    T = accumarray([i(:,1);i(:,1)],[i(:,2);i(:,3)],[n 1],@(x) {x});
    maxd2 = 0;
    R = zeros(n,2);
    move = false(n,1);
    for i=1:n
        Ti = T{i};
        P = X(i,:) - XALL(Ti,:);
        nP2 = sum(P.^2,2);
        if any(nP2<2*d2min)
            move(i) = true;
            move(Ti(Ti<=n)) = true;
        end
        maxd2 = maxd2 + mean(nP2);
        b = Ti > n;
        nP2(b) = nP2(b)*5; % reduce repulsion from each point of the border
        R(i,:) = sum(P./max((nP2-d2min),1e-3),1);
    end
    if ~any(move)
        break
    end
    if k==1
        v0 = (L*5e-3)/sqrt(maxd2/n);
    end
    R = R(move,:);
    v = v0/sqrt(max(sum(R.^2,2)));
    X(move,:) = X(move,:) + v*R;
    
    % Project back if points falling outside the rectangle
    X = min(max(X,Xmin),Xmax);
    
%     set(h,'XData',X(:,1),'YData',X(:,2));
    pause(0.01);
end
theta = linspace(0,2*pi,65);
xc = dmin/2*sin(theta);
yc = dmin/2*cos(theta);
% plot circles f diameter dmin around random points
% for i=1:n
%     plot(X(i,1)+xc,X(i,2)+yc,'k');
% end

end