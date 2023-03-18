function [k1f, k2f, rejRatio] = mcmc(nSimul)
%  N = number of simulations

f1 = fopen('data.txt');
fData = textscan(f1, '%f%f', 'HeaderLines', 1);
x = fData{1};
y = fData{2};


TW = 5; %initial temp of water
TA = 25; %initial temp of air
k1 = 10;
k2 = 5; %initial guesses for k1 k2

% fOfT = @(T,TW,TA) -k1*(T-TW) - k2*(T-TA)
func_ss = @(h) sum((y- - h(1).*(x-5) - h(2).*(x-25)).^2); %y = data, fOfT = function we want to fit

b_0 = [k1; k2]; %initial guesses

data = [x',y'];
n =length(x);

[bmin, ssmin] = fminsearch(func_ss, b_0);

chain = zeros(nSimul,2);
sigma2 = ssmin;
qcov = 0.1e-8*eye(2);
R=chol(qcov);

oldpar1 = bmin(:)';
chain(1,:) = oldpar1;
rej = 0;

ss = funcss(data, oldpar1(1), oldpar1(2));
SS = zeros(1,nSimul); 
SS(1) = ss;

for i=2:nSimul
    newPar1 = oldpar1 + randn(1,2)*R;
    %newPar2 = oldpar2 + randn(1,2)*R;
    ssNew = funcss(data, newPar1(1), newPar1(2));
    ssmin = ss;
    minVal = min(1, exp(-0.5*(ssNew-ssmin)/sigma2));
    u = rand;
    
    if u < minVal
        oldpar1 = newPar1;
        %oldpar2 = newPar2;
        SS(i-1) = SS(i);
        chain(i,:) = newPar1;
        %chain(i,2) = newPar2;
        ss =ssNew;
    else
        chain(i,:) = oldpar1;
        %chain(i,2) = oldpar2;
        rej = rej+1;
        ss= ssmin;
    end
    
end
% 
% li = exp(-1/(2*sigma.^2) * SS(theta))
% 
% mu = rand
% 
% if mu < exp((-1 / sigma.^2) * (SS(thetaNew) - SS(thetaOld))
%     thetaOld = thetaNew
%     SS(thetaOld) = SS(thetaNew)
% 

k1f = mean(chain(:,1));
k2f = mean(chain(:,2));
rejRatio = rej/nSimul;
% 
% plot(x,y);
% hold on
 plot(x, -k1f*(x-25)-k2f*(x-5), 'b')
