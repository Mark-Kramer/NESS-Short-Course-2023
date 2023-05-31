%% Load and visualize data
load MarkedDecodeExample.mat
figure; plot(1:T,x);
figure; plot(x(s),m,'.')

%% Spike sorting
id = kmeans(m',2); id1 = (id==1)'; id2 = (id==2)';
figure; plot(x(s(id1)),m(id1),'.',x(s(id2)),m(id2),'.')

%% Visualize receptive fields for sorted neurons
dx = 1e-2; xs = -8:dx:8; sigk = .25;
lambda1 = sum( normpdf(xs'*ones(1,sum(id1)),ones(size(xs'))*x(s(id1)),sigk),2)./sum( normpdf(xs'*ones(size(x)),ones(size(xs'))*x,sigk),2);
lambda2 = sum( normpdf(xs'*ones(1,sum(id2)),ones(size(xs'))*x(s(id2)),sigk),2)./sum( normpdf(xs'*ones(size(x)),ones(size(xs'))*x,sigk),2);
figure; plot(xs,lambda1,xs,lambda2);

%% Decode with sorted spikes
post = normpdf(xs'*ones(1,T),0,1); 
a = .98; sig = .3; OSM = normpdf(xs'*ones(size(xs)),ones(size(xs'))*a*xs,sig);

for t=2:T, 
  os = OSM*post(:,t-1); os = os/sum(os)/dx; 
  l = lambda1.^any(s(id1)==t).*exp(-lambda1).*lambda2.^any(s(id2)==t).*exp(-lambda2); 
  post(:,t) = os.*l; post(:,t) = post(:,t)/sum(post(:,t))/dx;
end;

figure; imagesc(1:T,xs,post); hold on; plot(1:T,x,'r');

%% Visualize joint mark intensity
dm = .1; ms = 9:dm:14; sigm = .6;
jmi = zeros(length(xs),length(ms));
for i=1:length(s);
    jmi=jmi+(normpdf(xs',x(s(i)),sigk)*normpdf(ms,m(i),sigm));%./(sum(normpdf(xs'*ones(size(x)),ones(size(xs'))*x,sigk),2)*ones(size(ms)));
end;
figure; imagesc(xs,ms,jmi'); set(gca,'YDir','normal'); axis([-4 4 9 14]);

%% Decode with JMI
post = normpdf(xs'*ones(1,T),0,1); 

groundIntensity = sum( normpdf(xs'*ones(size(s)),ones(size(xs'))*x(s),sigk),2)./sum( normpdf(xs'*ones(size(x)),ones(size(xs'))*x,sigk),2);

for t=2:T, 
  os = OSM*post(:,t-1); os = os/sum(os)/dx; 
  l = exp(-groundIntensity); if any(s==t), l=l .* (normpdf(xs'*ones(size(s)),ones(size(xs'))*x(s),sigk) * normpdf(m(find(s==t)),m',sigm))./sum( normpdf(xs'*ones(size(x)),ones(size(xs'))*x,sigk),2); end;
  post(:,t) = os.*l; post(:,t) = post(:,t)/sum(post(:,t))/dx;
end;

figure; imagesc(1:T,xs,post); hold on; plot(1:T,x,'r');
