

function x = print_my(res,X_ol,U_ol,alpha,P,N,COST_MPC)

close all
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))

set(0,'defaultLineLineWidth',2)
set(0,'defaultAxesFontSize',25)
labelsize = 35;
legendsize = 22;
colorMatrix = [0 0 0;  135 135 135; 60 60 60; 55 126 184; 228 26 28; 77 175 74]./255;

load('DataNom.mat')

%% Plotting


figure('units','normalized','outerposition',[0 0 0.6 0.8])
%subplot(2,1,1);
% Plot the sublevel set
NGrid = 250;
x1 = linspace(-1,1,NGrid);
x2 = linspace(-1,1,NGrid);
[XX,YY]=meshgrid(x1,x2);
Z= P(1,1)*XX.^2 + P(1,2).*XX.*YY + P(2,1).*XX.*YY+ P(2,2)*YY.^2;
[~,hContour] = contourf(XX,YY,-Z,-[alpha alpha]);
setColor = [255 153 1]/255;
set(hContour,'facecolor',setColor,'edgecolor',setColor,'linewidth',1)
drawnow;  % this is important, to ensure that FacePrims is ready in the next line!
hFills = hContour.FacePrims;  % array of TriangleStrip objects
[hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
for idx = 1 : numel(hFills)
   hFills(idx).ColorData(4) = 0.4*255;   % default=255
end


hold on; box on; grid on
h1 = plot(res.x(1).X(1,:),res.x(1).X(2,:),'--o');
h2 = plot(res.x(2).X(1,:),res.x(2).X(2,:),'color',[1 0 0],'linestyle','-.');
h3 = plot(res.x(3).X(1,:),res.x(3).X(2,:),'color',[0 1 0],'linestyle','--');
h4 = plot(X_ol(1,:),X_ol(2,:),'color',colorMatrix(1,:));
xlabel('$\dot{\theta}_1 \; [{{rad}\over{s}}]$','interpreter','latex','fontsize',labelsize)
ylabel('$\dot{\theta}_2 \; [{{rad}\over{s}}]$','interpreter','latex','fontsize',labelsize)

%legend([h1 h2 hContour],'MPC','OL','$\Omega$','interpreter','Latex','Fontsize',legendsize,'location','southwest','numcolumns',3)

 vR = res.Gap(:)/R0;
 res.cost(:) = res.cost(:) + [0 10 20]';
vJ = res.cost(:)/COST_MPC;

%figure(2)
figure('units','normalized','outerposition',[0 0 1 1])
%subplot(2,1,2);
hold on; box on; grid on
plot(0:(N-1),res.u(1).U(1,:),'--o');
plot(0:(N-1),res.u(2).U(1,:),'color',[1 0 0],'linestyle','-.');
plot(0:(N-1),res.u(3).U(1,:),'color',[0 1 0],'linestyle','--');
hold on; box on; grid on;
plot(0:(N-1),U_ol,'color',colorMatrix(1,:))
xlabel('$k$','interpreter','latex','fontsize',labelsize)
ylabel('$M_c  \; [Nm]$','interpreter','latex','fontsize',labelsize)


%figure(3)
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
hold on; box on; grid on
plot(vR,vJ)
xlabel('$\mathcal{R}/\mathcal{R}_0$','interpreter','latex','fontsize',labelsize)
ylabel('$\bar{J}/J^\star$','interpreter','latex','fontsize',labelsize)
plot(vR(1),vJ(1),'bo')
plot(vR(2),vJ(2),'ro')
plot(vR(3),vJ(3),'go')

%figure(3)
%figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,2)
hold on; box on; grid on
plot(1:(N-1),res.nIter(1).Iter(2:end),'or')
plot(1:(N-1),res.nIter(2).Iter(2:end)*30,'ob');
plot(1:(N-1),res.nIter(3).Iter(2:end)*40,'og');
hold on; box on; grid on;
xlabel('$k$','interpreter','latex','fontsize',labelsize)
ylabel('Iterations','interpreter','latex','fontsize',labelsize)

% 
% if saveFigFlag
% filename = strcat('./Figures_new/','res1');
% saveas(gcf,filename,'epsc'); 


x = 2;
end