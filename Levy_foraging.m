clear all;
close all;
clc;
%% Settings
c=1;
alpha=1.1;
N=50000;        % number of steps
M=50;           % number of targets
w=10;
rv=1;           % vision distance
D=200;
m=10000;
maxFlightDis = 10000;

%% setting for tempered Levy initialization
mu=1;
n = 1000;
xmax = 40;      % jump length max
xmin = rv;      % jump length min
xoffset = 0;
lambda = 1;

%% generating tempered Levy distribution
x = linspace(xmin,xmax,n);
p = @(xx) exp(-lambda.*xx).*(xx+xoffset).^(-mu);
px = @(xx) x.*exp(-lambda.*xx).*(xx+xoffset).^(-mu);
pxx = @(xx) x.*exp(-lambda.*xx).*(xx+xoffset).^(-mu);
psum = sum(p(x));
cdf = cumsum(p(x))/psum;
C = integral(p,xmin,xmax);
pdf = p(x)./C;

%% main loop
eta=0;         % eta? for efficiency?
iter=100;       % total number of iterations.
xx = [];
for p=1:iter
    num=0;
    dis=0;
    f=0;
    nn=1;
    
    t1=rand(1,M);
    t2=rand(1,M);
    tarx=D*t1;
    tary=D*t2;
    t1=tarx;
    t2=tary;
    
    % searcher start point
    starx(1)=D*rand(1);
    stary(1)=D*rand(1);
       
    for i=1:m
        % get jump length
        rn(i) = find(cdf >= rand,1,'first');
        
        % check if correct jump made
        while isnan(rn(i))
            rn(i) = find(cdf <= rand,1,'last');
        end
        
    end
    
    % truncate the jump length to max/min
    x=rn*(xmax-xmin)/n+xmin;
    
    % turn angle 
    t=2*pi*rand(1,N);
    
    % search the map
    for i=2:N
        
        % update the searcher location
        starx(i)=starx(i-1)+cos(t(i-1))*x(i-1);
        stary(i)=stary(i-1)+sin(t(i-1))*x(i-1);
        
        % check if targets are found
        for k=1:M
            
            % check if targets are within vision distance
            if point_to_line([tarx(k),tary(k),0],[starx(i-1),stary(i-1),0],[starx(i),stary(i),0])<rv
                starx(i)=tarx(k);
                stary(i)=tary(k);
                xx(nn)=tarx(k);
                yy(nn)=tary(k);
                nn=nn+1;
                tarx(k)=-100000;
                tary(k)=-100000;
                f=1;
                break;
            end
            
        end
        
        % update non periodic boundary conditions 
        if starx(i)>D||starx(i)<0
            starx(i)=starx(i-1);
        end
        
        if stary(i)>D||stary(i)<0
            stary(i)=stary(i-1);
        end
        
        % update flight length
        m=sqrt((starx(i)-starx(i-1))^2+(stary(i)-stary(i-1))^2);
        
        % update total flight distance
        dis=dis+m;
        
        % check if target found
        if f==1
            num=num+1;                      % update total targets found
%             plot(starx(i),stary(i),'*');    % update tar found plot 
%             hold on;
            f = 0;                          % reset f
        end
        
        
        if num==M || dis>maxFlightDis
            break;
        end
        
    end
    
    % add search efficiency for this iteration
    eta=eta+num/dis;
    xx = [xx,x];
end

% get average search efficiency
eta_avg = eta/iter; display(['average search efficiency: ',num2str(eta_avg)])
%%
% plot the last search
figure;
hold on;
plot(t1,t2,'o');
plot(xx(1:nn-1),yy(1:nn-1),'*');
plot(starx(1:i),stary(1:i));
% plot(tarx,tary,'*r');
plot(starx(1),stary(1),'^k','MarkerFaceColor','g');
plot(starx(i),stary(i),'^k','MarkerFaceColor','r');
hold off;

%%
nbins = 100;
nbinsx = 500;
xxx = [xx];
meanx = mean(xx);
stdx = std(xx);
pdfx = (xx-meanx)./stdx;


yyaxis left
hist(x,50);
yyaxis right
[countsx,binsx] = hist(xxx,nbinsx);
countsx = countsx./((binsx(2)-binsx(1))*sum(countsx));
plot(binsx,countsx.*C,'g'); hold on

xn = randn(1,length(xxx));
[counts,bins] = hist(xn,nbins);
counts = counts./((bins(2)-bins(1))*sum(counts));
plot(bins,counts,'g--'); hold off;

legend('hist','tempered levy','Gaussian')
xlim([-10,10])
% ylim([0,1.4])

