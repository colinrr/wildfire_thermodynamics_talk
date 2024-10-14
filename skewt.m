function skewt(pressures,temps,dewpts,varargin)
%
% Required arguments:
%   pressures - Pressure levels in hPa
%   temps - Temperature in C
%   dewpts - Dew point in C
%
% Optional arguments:
%   u - u wind speed (data displayed in given units)
%   v - v wind speed (same unit as u)
%
%   'IsMagDir'  - 1 Specify wind as magnitude and direction (degrees)
%               - 0 (default) Specify wind as U and V components
%   'IsMR'      - 1 Specify moisture as mixing ratio (kg/kg)
%               - 0 (default) Specify moisture as dew points (C)
%   'DrawLines' - 1 (default) Make new figure and draw thermo diagram.
%               - 0 Don't make a new figure and don't draw the thermo diagram.
%                 Draw lines on top of existing figure
%   'PlotBW'    - Plot chart in color (0, default) or in black and white (1)
%   'TempOpts'  - Linestyle options for the temperature profile
%                 e.g. 'ko--' or {'color',[0 .5 0],'linewidth',4}
%   'Dwptopts'  - Linestyle options for the water vapor profile
%
% Examples:
%   Plot temperature and dewpoint data at 2 levels
%       skewt([1000 500],[10 -5],[8 -20])
%
%   Plot temperature, dewpoint, and wind data at 2 levels
%       skewt([1000 500],[10 -5],[8 -20],[5 20],[5 -20])
%
%   Plot temperature, dewpoint, and wind data given in magnitude and 
%   direction at 2 levels
%       skewt([1000 500],[10 -5],[8 -20],[5 20],[45 180],'IsMagDir',1)
%
%   Plot temperature and dewpoint data at 2 levels with custom line specs
%   for temperature
%       skewt([1000 500],[10 -5],[8 -20],'TempOpts',{'color',[0 .5 0],'linewidth',4)
if nargin<3
    error('Too few input arguments')
end
if length(varargin)==1
    error('Only one wind component specified or not enough arguments')
end
defaultDrawlines=1;
defaultGridLinesWidth = 1.0;
defaultDataLineWidth = 1.3;
defaultTempopts={'k-','LineWidth',defaultDataLineWidth};
defaultDwptopts={'-','Color',[0 0.4470 0.7410],'LineWidth',defaultDataLineWidth};
defaultPlotBW=0;
defaultIsMagDir=0;
defaultIsMR=0;
p=inputParser;
addRequired(p,'pressures');
addRequired(p,'temps');
addRequired(p,'dewpts');
if length(varargin)>=2 && isnumeric(varargin{1})
    addRequired(p,'u');
    addRequired(p,'v');
    addOptional(p,'DrawLines',defaultDrawlines,@isscalar);
    addOptional(p,'TempOpts',defaultTempopts);
    addOptional(p,'DwptOpts',defaultDwptopts);
    addOptional(p,'PlotBW',defaultPlotBW,@isscalar); 
    addOptional(p,'IsMagDir',defaultIsMagDir,@isscalar); 
    addOptional(p,'IsMR',defaultIsMR,@isscalar); 
    if length(varargin)==2
        varargin=[];
        parse(p,pressures,temps,dewpts,varargin{1},varargin{2},varargin{:})
    else
        parse(p,pressures,temps,dewpts,varargin{1},varargin{2},varargin{3:end})
    end
    u=varargin{1};
    v=varargin{2};
else
    u=[];
    v=[];
    addOptional(p,'DrawLines',defaultDrawlines,@isscalar);
    addOptional(p,'TempOpts',defaultTempopts);
    addOptional(p,'DwptOpts',defaultDwptopts);
    addOptional(p,'PlotBW',defaultPlotBW,@isscalar);
    addOptional(p,'IsMagDir',defaultIsMagDir,@isscalar);
    addOptional(p,'IsMR',defaultIsMR,@isscalar); 
    parse(p,pressures,temps,dewpts,varargin{:})    
end
drawlines=p.Results.DrawLines;
tempopts=p.Results.TempOpts;
if ~iscell(tempopts)
    tempopts={tempopts};
end
dwptopts=p.Results.DwptOpts;
if ~iscell(tempopts)
    dwptopts={dwptopts};
end
plotbw=p.Results.PlotBW;
isMagDir=p.Results.IsMagDir;
isMR=p.Results.IsMR;
%Some constants
T0  = 273.15;
R   = 287;
Lv	= 2.5e6;
cp  = 1004;
H   = 10;
R_v = 461;
epsilon = R/R_v;
%convert mixing ratio to dewpoint if necessary
if isMR
    dewpts=dewpts./(dewpts+epsilon).*pressures*100;
    dewpts=-5.42e3./log(dewpts/2.53e11)-T0;
end
%chart bottom pressure
cbp=1025;
%chart top pressure
ctp=50;
%Evenly spaced vector of height values, then converted to pressure
z=0:.5:30;  
p=cbp*exp(-z/H);
nz=length(z);

%% CR - redefine chart limits and params to be flexible...
nj = 190; %Range of T
offset = -40; %Minimum chart T

%Isotherm values
Tcontours=-150:10:(nj-offset);
%Dry adiabat values
thcontours=-200:20:300;
% thcontours=2.*min(Tcontours):10:2*max(Tcontours);
%Moist adiabat values
% thecontours=thcontours(1:end-5)-5;
thecontours = [-200:20:80 logspace(2,2.4771,5) logspace(2.4771,3,4)]; % 180:20:300 ]; % 120 150 200 300]; % logspace(2,3,5) logspace(4,6,5)];
%Mixing ratio line values
% rvcontours=[.1,.5,1,2,3,5,7,10,15,22,30,50,100,200,400,800];
rvcontours=[.1,1,3,5,10,22,50,200,800];
%%
%y-axis tick locationspepsi
paxis=[cbp,1000:-50:100];
% nj=85; %Range of T
% offset=-40; %Minimum chart T
[T,es,ws,th,the]=deal(zeros(nj,nz));
for i=1:nj
    for j=1:nz
        %Create temperature matrix
        T(i,j)=i-2*j+offset;
        
        %Saturation vapor pressure
%         es(i,j) = 2.53e11*exp(-5.42e3./(T(i,j)+T0)); % Too inaccurate at high T
        es(i,j) = 100 * 6.112 .* exp(17.67.*(T(i,j))./(T(i,j) + 243.5));
        
        %Saturation mixing ratio
        ws(i,j)=epsilon*es(i,j)./(p(j)*100-es(i,j));
        
        %Dry adiabats
        th(i,j)=(T(i,j)+T0)*(1000/p(j))^(R/cp);
        
        %Moist adiabats
        the(i,j)=th(i,j).*exp(Lv.*ws(i,j)/cp./(T(i,j)+T0));
    end
end
% %Isotherm values
% Tcontours=-100:10:100;
% %Dry adiabat values
% thcontours=-200:10:200;
% %Moist adiabat values
% thecontours=thcontours(1:end-5)-5;
% %Mixing ratio line values
% rvcontours=[.1,.5,1,2,3,5,7,10,15,22,30];
%Begin plotting the chart
if drawlines==1
    figure
    hold on
    if plotbw==0 %plot lines in color
        %isotherms
        contour(T(:,1),p,T',Tcontours,'-','color',[.3 .3 .3])
        %dry adiabats
        contour(T(:,1),p,th'-T0,thcontours,'-','color',[.55 .55 .55])
        %saturation mixing ratio lines
        contour(T(:,1),p,ws'*1000,rvcontours,'--','color',[.48 .75 .48])
        %moist adiabats
        contour(T(:,1),p,the'-T0,thecontours,'-','color',[.7 .85 1.0])
    else
        contour(T(:,1),p,T',Tcontours,'k')
        contour(T(:,1),p,th'-T0,thcontours,'k')
        contour(T(:,1),p,ws'*1000,rvcontours,'--','color',[.3 .3 .3])
        contour(T(:,1),p,the'-T0,thecontours,'-','color',[.7 .7 .7])
    end
    
    %Make mixing ratio labels around the edge of the plot
    rvlabel = interp1(ws(4,:)*1000,p,rvcontours(1));
    text(T(4,1),rvlabel,num2str(rvcontours(1)),...
        'HorizontalAlignment','center',...
        'fontsize',get(gca,'fontsize'),'color',[.3 .3 .3]);
    
    rvlabel = interp1(ws(:,2)*1000,T(:,1),rvcontours);
    for i=2:length(rvcontours)
        text(rvlabel(i),paxis(2),num2str(rvcontours(i)),...
            'HorizontalAlignment','center',...
            'fontsize',get(gca,'fontsize'),'color',[.3 .3 .3]);
    end
    
    %Set y-axis
    set(gca,'YGrid','on')
    set(gca,'GridLineStyle','-')
    set(gca,'ydir','reverse','yscale','log','ylim',[ctp cbp],'YTick',sort(paxis))
    set(gca,'YTickLabel',{'100','','200','','300','','400','','500','','600',...
        '','700','','','850','','','1000',''})
    
    ylabel(' Pressure (mb) ')
    
    %Set x-axis
    xlabel(' Temperature (^oC) ')
    set(gca,'xlim',[offset offset+nj-5])
    
    box on
    
end %End plotting the chart
%--------------------------------------------------------------------------
%Now plot the data
%Do temperature and dewpoint first
inds=find(isnan(pressures)+isnan(temps)+isnan(dewpts)==0);
pressures=pressures(inds);
temps=temps(inds);
dewpts=dewpts(inds);
yax = H*log(cbp./pressures);
xax=temps+4*yax;
xax2=dewpts+4*yax;
plot(xax,pressures,tempopts{:});
plot(xax2,pressures,dwptopts{:});
%--------------------------------------------------------------------------
%Now create wind barbs on the right
%Based on this file:
%  MFILE:   windbarbm.m
%  MATLAB:  7.8.0 (R2009a)
%  VERSION: 1.3 (28 November 2011)
%  AUTHOR:  Nick Siler
%  CONTACT: siler@atmos.washington.edu
if ~isempty(u)
    %This is important only if wind barbs are displayed
    axis square
    
    u=u(pressures>=ctp);
    v=v(pressures>=ctp);
    pressures=pressures(pressures>=ctp);
    
    xloc = offset+nj+3; %change this to change the position of the wind barb axis
    
    xx=ones(size(u))*xloc;
    yy=H*log(cbp./pressures);
    
    if ~isMagDir
        umag = (u.^2+v.^2).^0.5; %wind speed
        %find theta; add pi to atan(v/u) when u<0
        dummy = (u<0)*pi;
        theta = atan(v./u)+dummy;
    else
        umag=u;
        theta=v*pi/180-pi/2;
    end
    
    [a,b] = size(umag);
    
    %create 18 logical matrices for 18 possible barbs. Non-zero when the barb
    %is called for at that gridpoint.
    g{1} = umag > 7.5 & umag <= 47.5;
    g{2} = umag > 17.5 & umag <= 47.5;
    g{3} = umag > 27.5;
    g{4} = (umag > 37.5 & umag <= 47.5) | (umag > 57.5 & umag <= 97.5);
    g{5} = umag > 67.5;
    g{6} = (umag > 77.5 & umag < 97.5) | umag > 107.5;
    g{7} = umag > 87.5 & umag < 97.5 | umag > 117.5;
    g{8} = umag > 127.5;
    g{9} = (umag > 2.5 & umag <= 7.5) | (umag > 12.5 & umag <= 17.5);
    g{10} = umag > 22.5 & umag <= 27.5;
    g{11} = (umag > 32.5 & umag <= 37.5) | (umag > 52.5 & umag <= 57.5);
    g{12} = (umag > 42.5 & umag <= 47.5) | (umag > 62.5 & umag <= 67.5);
    g{13} = (umag > 72.5 & umag <= 77.5) | (umag > 102.5 & umag <= 107.5);
    g{14} = (umag > 82.5 & umag <= 87.5) | (umag > 112.5 & umag <= 117.5);
    g{15} = (umag > 92.5 & umag <= 97.5) | (umag > 122.5 & umag <= 127.5);
    g{16} = umag > 47.5;
    g{17} = umag > 97.5;
    g{18} = true(a,b);
    
    %position of each barb relative to grid point: [x0 y0; x1 y1]
    c{1} = [-1 0;-1.125 .325];
    c{2} = [-.875 0; -1 .325];
    c{3} = [-.75 0; -.875 .325];
    c{4} = [-.625 0; -.75 .325];
    c{5} = [-.5 0; -.625 .325];
    c{6} = [-.375 0; -.5 .325];
    c{7} = [-.25 0; -.375 .325];
    c{8} = [-.125 0; -.25 .325];
    c{9} = [-.875 0; -.9375 .1625];
    c{10} = [-.75 0; -.8125 .1625];
    c{11} = [-.625 0; -.6875 .1625];
    c{12} = [-.5 0; -.5625 .1625];
    c{13} = [-.3750 0; -.4375 .1625];
    c{14} = [-.25 0; -.3125 .1625];
    c{15} = [-.125 0; -.1875 .1625];
    c{16} = [-1 0; -.875 .325];
    c{17} = [-.75 0; -.625 .325];
    c{18} = [0 0; -1 0];
    
    for nn=1:18
        c{nn}=c{nn}*2;
    end
    
    scale2x=3;
    scale2y=scale2x*diff(H*log(cbp./[cbp ctp]))/(nj-5);
    
    for nn = 1:18
        dummy = reshape(g{nn},1,a*b);
        count = sum(dummy); % number of barbs to draw
        if count == 0
            continue
        end
        
        %rotation operations
        x1 = c{nn}(1,1)*cos(theta)-c{nn}(1,2)*sin(theta);
        y1 = c{nn}(1,1)*sin(theta)+c{nn}(1,2)*cos(theta);
        x2 = c{nn}(2,1)*cos(theta)-c{nn}(2,2)*sin(theta);
        y2 = c{nn}(2,1)*sin(theta)+c{nn}(2,2)*cos(theta);
        
        x1 = x1*scale2x+xx;
        x2 = x2*scale2x+xx;
        
        y1 = y1*scale2y+yy;
        y2 = y2*scale2y+yy;
        
        x = [reshape(x1(dummy),1,count);reshape(x2(dummy),1,count)];
        y = [reshape(y1(dummy),1,count);reshape(y2(dummy),1,count)];
        y = cbp*exp(-y/H);
               
        line(x,y,'color','k','linestyle','-','linewidth',1,'clipping','off','marker','none')
        if nn==18
            plot(x(1,:),y(1,:),'o','color','k','MarkerSize',2,'clipping','off')
        end
    end
    plot([xloc xloc],[ctp cbp],'color','k','linestyle','-','linewidth',1,'clipping','off')
end