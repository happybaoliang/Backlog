function []=singlebacklog()
% experimentI();
% close all;
% experimentII();
% close all;
% experimentIII();
% close all;
experimentIV();
end

function experimentI()
C=1000;
lambda=800;
maximum=15;
buf_size=10;

tau=0.000;
[X1,Y1]=backlog_bound_calculus_MD1(lambda,C,tau,buf_size,maximum);
plot(X1,Y1,'-.k');
hold on;

tau=0.002;
[X3,Y3]=backlog_bound_calculus_MD1(lambda,C,tau,buf_size,maximum);
plot(X3,Y3,'-.kv');
hold on;

tau=0.004;
[X5,Y5]=backlog_bound_calculus_MD1(lambda,C,tau,buf_size,maximum);
plot(X5,Y5,'-.ko');
hold on;

tau=0.006;
[X7,Y7]=backlog_bound_calculus_MD1(lambda,C,tau,buf_size,maximum);
plot(X7,Y7,'-.kx');

xlabel('x (kb)','fontsize',14)
ylabel('backlog bound','fontsize',14)
axis([0 maximum 0 0.2])
handle=legend('\tau=0.000s','\tau=0.002s','\tau=0.004s','\tau=0.006s');
set(handle,'FontSize',14);
print -depsc ../figures/backlogtau.eps
end

function experimentII()
W=5;
mu=10;
lambda=8;
maximum=20;

[X1,Y1]=sliding_window_bound(lambda,mu,W,maximum);
plot(X1,Y1,'-kx');
hold on;
[X2,Y2]=martingale_MM1_bound(lambda,mu,W,maximum);
plot(X2,Y2,'-ko');
hold on;

W=10;
[X3,Y3]=sliding_window_bound(lambda,mu,W,maximum);
plot(X3,Y3,'-kv');
hold on;
[X4,Y4]=martingale_MM1_bound(lambda,mu,W,maximum);
plot(X4,Y4,'-k+');
hold on;

W=15;
[X5,Y5]=sliding_window_bound(lambda,mu,W,maximum);
plot(X5,Y5,'-k<');
hold on;
[X6,Y6]=martingale_MM1_bound(lambda,mu,W,maximum);
plot(X6,Y6,'-k.');

xlabel('x (kb)','fontsize',14)
ylabel('backlog bound','fontsize',14)
axis([0 maximum 0 0.5])
handle=legend('Eq.(19)-W=5kb','Eq.(20)-W=5kb','Eq.(19)-W=10kb','Eq.(20)-W=10kb','Eq.(19)-W=15kb','Eq.(20)-W=15kb');
set(handle,'FontSize',14);
print -depsc ../figures/backlogcomp.eps
end

function experimentIII()
tau=0.00;

prob=0.75;
maximum=50;
timeslot=0.02;
servicetime=0.025;

buf_size=5;
[X1,Y1]=backlog_bound_martingale_bernoulli(prob,timeslot/servicetime,buf_size,maximum);
plot(X1,Y1,'-.k');
hold on;
filename=sprintf('../dataset/Bernoulli_queuelength(T=%1.2fs,p=%1.2f,L=1kb,mu=%1.2fs,W=%02dkb,tau=%1.3fs).csv',timeslot,prob,timeslot/servicetime,buf_size,tau);
[X2,Y2]=backlog_bound_simulation(filename);
plot(X2,Y2,'-k');
hold on;

buf_size=10;
[X3,Y3]=backlog_bound_martingale_bernoulli(prob,timeslot/servicetime,buf_size,maximum);
plot(X3,Y3,'-.kv');
hold on;
filename=sprintf('../dataset/Bernoulli_queuelength(T=%1.2fs,p=%1.2f,L=1kb,mu=%1.2fs,W=%02dkb,tau=%1.3fs).csv',timeslot,prob,timeslot/servicetime,buf_size,tau);
[X4,Y4]=backlog_bound_simulation(filename);
plot(X4,Y4,'-kv');
hold on;

buf_size=15;
[X5,Y5]=backlog_bound_martingale_bernoulli(prob,timeslot/servicetime,buf_size,maximum);
plot(X5,Y5,'-.ko');
hold on;
filename=sprintf('../dataset/Bernoulli_queuelength(T=%1.2fs,p=%1.2f,L=1kb,mu=%1.2fs,W=%02dkb,tau=%1.3fs).csv',timeslot,prob,timeslot/servicetime,buf_size,tau);
[X6,Y6]=backlog_bound_simulation(filename);
plot(X6,Y6,'-ko');
hold on;

buf_size=20;
[X7,Y7]=backlog_bound_martingale_bernoulli(prob,timeslot/servicetime,buf_size,maximum);
plot(X7,Y7,'-.kx');
filename=sprintf('../dataset/Bernoulli_queuelength(T=%1.2fs,p=%1.2f,L=1kb,mu=%1.2fs,W=%02dkb,tau=%1.3fs).csv',timeslot,prob,timeslot/servicetime,buf_size,tau);
[X8,Y8]=backlog_bound_simulation(filename);
plot(X8,Y8,'-kx');

xlabel('x (kb)','fontsize',14)
ylabel('backlog bound','fontsize',14)
axis([0 maximum 0 1])
handle=legend('Theorem 5-W=5kb','simulation-W=5kb','Theorem 5-W=10kb','simulation-W=10kb','Theorem 5-W=15kb','simulation-W=15kb','Theorem 5-W=20kb','simulation-W=20kb');
set(handle,'FontSize',14);
print -depsc ../figures/backlogbuf.eps
end

function experimentIV()
t=1;
mu=1;
maximum=50;
buf_size=5;

prob=0.8;
[X1,Y1]=backlog_bound_effective_Bernoulli(prob,mu,t,buf_size,maximum);
plot(X1,Y1,'-.k');
hold on;

prob=0.85;
[X3,Y3]=backlog_bound_effective_Bernoulli(prob,mu,t,buf_size,maximum);
plot(X3,Y3,'-.kv');
hold on;

prob=0.9;
[X5,Y5]=backlog_bound_effective_Bernoulli(prob,mu,t,buf_size,maximum);
plot(X5,Y5,'-.ko');
hold on;

prob=0.95;
[X7,Y7]=backlog_bound_effective_Bernoulli(prob,mu,t,buf_size,maximum);
plot(X7,Y7,'-.kx');

xlabel('x (kb)','fontsize',14)
ylabel('backlog bound','fontsize',14)
axis([0 maximum 0 1])
handle=legend('\rho=0.80','\rho=0.85','\rho=0.90','\rho=0.95');
set(handle,'FontSize',14);
print -depsc ../figures/backlogrho.eps
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y]=backlog_bound_calculus_MD1(lambda,C,tau,buf_size,maximum)
if nargin~=5
    error('wrong number of input arguments.');
else
    X=0:1:maximum;
    Y=0:1:maximum;
    
    func=@(x) lambda*(exp(x)-1)/x-C;
    p=fzero(func,1);
    
    if (p<=0.000001)
        error('failed to find a zero point');
    end
    
    for i=1:1:maximum+1
        min=100000000000;
        for theta=0.000001:0.001:p
            temp=exp(-theta*(X(i)+buf_size-tau*lambda*(exp(theta)-1)/theta));
            if temp<min
                min=temp;
            end
        end
        Y(i)=min;
    end
end
end

function [X,Y]=sliding_window_bound(lambda,mu,W,maximum)
temp=0;
rho=lambda/mu;
X=0:1:maximum;
Y=0:1:maximum;

for i=1:1:W
    temp=temp+rho^i;
end

P00=1/(1+temp+(rho^(W+1))/(1-rho));

for i=0:1:maximum
    Y(i+1)=P00*(rho^(W+i))/(1-rho);
end
end

function [X,Y]=martingale_MM1_bound(lambda,mu,W,maximum)
X=0:1:maximum;
Y=0:1:maximum;

theta_max=log(mu/lambda);

for i=1:1:maximum+1
    min=10000000000;
    for theta=0.001:0.001:theta_max+0.0001
        y=exp(lambda*(exp(theta)-1)+mu*(exp(-theta)-1)-theta*(W+X(i)));
        if y<min
            min=y;
        end
    end
    if min>1
        Y(i)=1;
    else
        Y(i)=min;
    end
end
end

function [X,Y]=backlog_bound_simulation(filename)
Delay=csvread(filename,1,1);
len=max(Delay);
X=0:1:len;
Y=0:1:len;
disp('entering simulation');
for i=0:1:len
    Y(i+1)=length(find(Delay>i))/length(Delay);
end
disp('leaving simulation');
end

function [X,Y]=backlog_bound_martingale_bernoulli(p,mu,buf_size,maximum)
X=0:1:maximum;
Y=0:1:maximum;

func=@(x) (1-p+p*exp(x))*exp(mu*(exp(-x)-1))-1;
theta_max=fzero(func,0.2);
if theta_max<=0
    error('theta_max<=0');
end

for i=1:1:maximum+1
    min=100000000000;
    for theta=0.00001:0.00001:theta_max+0.000001
        y=(1-p+p*exp(theta))*exp(mu*(exp(-theta)-1))*exp(-theta*(buf_size+X(i)));
        if min>y
            min=y;
        end
    end
    if min>1
        Y(i)=1;
    else
        Y(i)=min;
    end
end
end

function [X,Y]=backlog_bound_effective_Bernoulli(p,mu,t,buf_size,maximum)
if nargin~=5
    error('wrong number of input arguments.');
else
    X=0:1:maximum;
    Y=0:1:maximum;
    
    func=@(x) log(1-p+p*exp(x))./x+mu*(exp(-x)-1)./x;
    theta_max=fzero(func,5);
    if theta_max<=0
        error('theta_max<=0');
    end
    
    for i=1:1:maximum+1
        min=1000000000000;
        for theta=0.001:0.01:theta_max+0.001
            temp=exp(t*log(1-p+p*exp(theta))-theta*(X(i)+buf_size));
            if temp<min
                min=temp;
            end
        end
        if min>1
            Y(i)=1;
        else
            Y(i)=min;
        end
    end
end
end