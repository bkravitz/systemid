function am=anmean(a);
% compute annual mean of variables
% note all data starts in December... so these aren't aligned with years
nyr=length(a)/12;
am=zeros(nyr,1);
for k=1:nyr,
    am(k)=mean(a((k-1)*12+[1:12]));
end