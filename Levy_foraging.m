clear all;
close all;
clc;
%% Settings
c=1;
alpha=1.1;
N=50000;
M=50;
w=10;
rv=1;
D=200;
m=10000;

%% setting for tempered Levy initialization
mu=1;
n = 1000;
xmax = 40;
xmin = rv;
xoffset = 0;
lambda = 1;

%% generating tempered Levy distribution
x = linspace(xmin,xmax,n);
p = @(xx) exp(-lambda.*xx).*(xx+xoffset).^(-mu);
psum = sum(p(x));
cdf = cumsum(p(x))/psum;
C = integral(p,xmin,xmax);
pdf = p(x)./C;

%% main loop
aita=0;         % eta? for efficiency?
iter=100;       % total number of iterations.
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
        
        rn(i) = find(cdf >= rand,1,'first');
        
        % check if correct jump made
        while isnan(rn(i))
            rn(i) = find(cdf <= rand,1,'last');
        end
        
    end
    
    x=rn*(xmax-xmin)/n+xmin;
    t=2*pi*rand(1,N);
    
    for i=2:N
        
        starx(i)=starx(i-1)+cos(t(i-1))*x(i-1);
        stary(i)=stary(i-1)+sin(t(i-1))*x(i-1);
        
        for k=1:M
            
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
        
        if starx(i)>D||starx(i)<0
            starx(i)=starx(i-1);
        end
        
        if stary(i)>D||stary(i)<0
            stary(i)=stary(i-1);
        end
        
        m=sqrt((starx(i)-starx(i-1))^2+(stary(i)-stary(i-1))^2);
        dis=dis+m;
        
        if f==1
            num=num+1;
            plot(starx(i),stary(i),'*');
            hold on;
        end
        
        f=0;
        
        if num==M || dis>10000
            break;
        end
        
    end
    
    aita=aita+num/dis;
    
end

%%
aita/iter
figure;
hold on;
plot(t1,t2,'o');
plot(xx(1:nn-1),yy(1:nn-1),'*');
plot(starx(1:i),stary(1:i));
hold off;
% text(tarx,tary,'*','color','r');
% text(starx(1),stary(1),'*','color','r');


