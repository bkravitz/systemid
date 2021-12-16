% Code to divide the ocean regions into boxes of similar area
% and then create random, uncorrelated sequences for each box
%
% Created by Ben Kravitz (ben.kravitz@pnnl.gov or ben.kravitz.work@gmail.com)
% Last updated 20 January 2015
%
% This is for CESM (2 degree resolution) with 144 longitudes and 96 latitudes.
% Many things are hard-coded.  Port this script at your own risk.

clear all;

%% Input values

nboxes=20; % approximate number of regions (the actual amount will vary slightly)
lfthreshold=0.2; % The threshold of land fraction - above this number, a grid box is
                 % considered land and is omitted
nyrs=5; % how long you want the simulation to be (in years)
nt=7300; % how many elements in the sequence you want (for example, if you want the
         % sequence to shift daily, this number should be nyrs*365)
tmscl=14; % the timescale (in days) of the processes that you care about
ampl=1; % This gives you the option to scale the random sequences so that the largest
           % magnitude in any sequence for any grid box is equal to this

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% MODIFY ANYTHING BELOW AT YOUR OWN RISK %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Creating boxes

ncid=netcdf.open('sftlf_CESM2x2.nc','NC_NOWRITE');
latid=netcdf.inqVarID(ncid,'lat');
lonid=netcdf.inqVarID(ncid,'lon');
fracid=netcdf.inqVarID(ncid,'LANDFRAC');
lat=netcdf.getVar(ncid,latid);
lon=netcdf.getVar(ncid,lonid);
lf=netcdf.getVar(ncid,fracid);
netcdf.close(ncid);

thegrid=zeros(144,96,2);
thegrid(:,:,1)=repmat(lon,[1 96]);
thegrid(:,:,2)=repmat(transpose(lat),[144 1]);

gbspace=diff(lat);gbspace=gbspace(1);
latedges=[-90 ; lat(1:end-1)+gbspace/2 ; 90];
colatedges=pi/2-latedges*pi/180;
dphi=2.5*pi/180;
r=6371000;
areaout=abs(diff(cos(colatedges)))*dphi.*(r.^2);
gba=repmat(transpose(areaout),[144 1]);

oceans=zeros(size(gba));

% We are defining the Southern Ocean as all ocean regions below and including
% latitude 20 (defined as the Southern tip of South America)
% Southern Ocean has index 1
oceans(:,1:20)=1;

% The Arctic Ocean is everything North of the Bering Strait (latitude 83)
% Arctic Ocean has index 2
oceans(:,83:96)=2;

% The Indian Ocean Western Boundary is the middle of the Southern tip of Africa
% (longitude 10)
% The Eastern most boundary is defined by Australia and the Maritime Continent
% straight line from (46,61) to (55,42) [NOTE - INCORRECT]
% South of Australia along longitude 55 from the Southern Ocean to latitude 35
% None of the Indian Ocean is North of latitude 65
% Indian Ocean has index 3
oceans(11:46,21:65)=3;
oceans(45:54,21:42)=3;
for m=42:51;
    for n=41:58;
        if n<(((58-41)/(42-51))*(m-42)+58);
            oceans(m,n)=3;
        end
    end
end

% The Pacific Ocean has an Western boundary that is the same as the Eastern boundary
% of the Indian Ocean
% All Pacific Ocean falls in longitude 46 to 117
% All Pacific Ocean falls in latitude 21 to 82
% Central American boundary line goes from (104,60) to (114,52)
% Pacific Ocean has index 4
oceans(45:104,61:82)=4;
for m=42:55;
    for n=41:61;
        if n>=(((58-41)/(42-51))*(m-42)+58);
            oceans(m,n)=4;
        end
    end
end
for m=104:114;
    for n=52:60;
        if n<=(((60-52)/(104-114))*(m-104)+60);
            oceans(m,n)=4;
        end
    end
end
oceans(55:104,21:61)=4;
oceans(104:117,21:52)=4;

% The Atlantic Ocean is everything else
oceans(oceans==0)=5;

% Hudson Bay is part of the Arctic
oceans(106:116,76:84)=2;

% Mediterranean, Black Sea, etc. get their own label
oceans(1:25,64:75)=6;
oceans(7:10,76:83)=6;

I=find(lf>lfthreshold);
oceans2=oceans;
oceans2(I)=0;

figure;
pcolor(transpose(oceans2));
title('Ocean Basins','FontSize',15);

I=find(isnan(oceans2)==1);
outputmatrix(I)=-1e30;
fid=fopen('oceans2.bin','wb');
fwrite(fid,oceans2,'single','b');
fclose(fid);

J=find(oceans2>0);
totaloceanarea=sum(gba(J));
regareas=zeros(1,6);
for k=1:6;
    J=find(abs(oceans2-k)<0.5);
    regareas(k)=sum(gba(J));
end

regboxes=nboxes*regareas/totaloceanarea;
J=0;
while isempty(J)==0;
    J=find(regboxes<1);
    K=find(regboxes>=1);
    c=sum(1-regboxes(J));
    regboxes(J)=1;
    regboxes(K)=regboxes(K)-c/length(K);
end
regboxes=round(regboxes);

latwts=gba(1,:);
oceans3=oceans2;
for k=1:6;
    % simple case = only one region
    nreg=regboxes(k);
    if nreg==1;
        continue;
    end

    % normal case = multiple regions
    indices=k:(1/nreg):k+1-1/nreg/10;
    
    % first divide each region into equal area latitude bands described by aspect ratio
    widths=zeros(1,96);
    heights=zeros(1,144);
    for j=1:96;
        I=find(abs(oceans2(:,j)-k)<0.5);
        if k~=5;
            leftmost=min(I);
            rightmost=max(I);
            lonwts=gba(:,j);
            widths(j)=sum(lonwts(leftmost:rightmost));
        else
            J=I;
            J(J<72)=J(J<72)+144;
            leftmost=min(J);
            rightmost=max(J);
            lonwts=cat(1,gba(:,j),gba(:,j));
            widths(j)=sum(lonwts(leftmost:rightmost));
        end
    end
    for j=1:144;
        I=find(abs(oceans(j,:)-k)<0.5);
        botmost=min(I);
        topmost=max(I);
        heights(j)=sum(latwts(botmost:topmost));
    end
    medwidth=median(widths(widths>0));
    medheight=median(heights(heights>0));
    aspectratio=medheight/medwidth;
    aspsout=zeros(1,nreg);
    for j=1:nreg;
        aspsout(j)=j^2/nreg;
    end
    [~,nrows]=min(abs(aspsout-aspectratio));
    
    % defining the latitude bands
    I=find(abs(oceans2-k)<0.5);
    baselats=squeeze(thegrid(:,:,2));
    havelats=baselats(I);
    botlat=min(havelats);
    toplat=max(havelats);
    botlatind=find(lat==botlat);
    toplatind=find(lat==toplat);
    gbakeep=zeros(size(gba));
    gbakeep(I)=gba(I)./gba(I);
    baselatareas=sum(gba(:,botlatind:toplatind).*gbakeep(:,botlatind:toplatind),1);
    eachband=regareas(k)/nrows;
    rowinds=cell(1,nrows);
    c=1;
    for j=1:nrows;
        X=cumsum(baselatareas(c:end));
        [~,I]=min(abs(X-eachband));
        rowinds{j}=c:c+I-1;
        c=c+I;
    end
    
    % subdividing the latitude bands into boxes of similar longitude range
    c=0;
    f=1;
    for j=1:nrows;
        boxesremaining=nreg-c;
        boxesthisrow=min([ceil(nreg/nrows) boxesremaining]);
        c=c+boxesthisrow;
        worklats=rowinds{j};
        uselats=botlatind+worklats(1)-1:botlatind+worklats(end)-1;
        baselons=squeeze(thegrid(:,uselats,1));
        I=find(abs(oceans2(:,uselats)-k)<0.5);
        havelons=baselons(I);
        if k~=5;
            leftlon=min(havelons);
            rightlon=max(havelons);
            leftlonind=find(lon==leftlon);
            rightlonind=find(lon==rightlon);
            nlons=length(leftlonind:rightlonind);
            eachbox=nlons/boxesthisrow;
            theboxes=cell(1,boxesthisrow);
            d=1;
            for q=1:boxesthisrow;
                theboxes{q}=d:min([d+ceil(eachbox)-1 144]);
                d=d+ceil(eachbox);
            end
            for q=1:boxesthisrow;
                currindex=indices(f);
                worklons=theboxes{q};
                uselons=leftlonind+worklons(1)-1:leftlonind+worklons(end)-1;
                usepairs=zeros(length(uselons)*length(uselats),2);
                usepairs(:,1)=repmat(transpose(uselons),[length(uselats) 1]);
                blah=repmat(uselats,[length(uselons) 1]);
                usepairs(:,2)=blah(:);
                I=find(abs(oceans2-k)<0.5);
                oceans2pairs=zeros(length(I),2);
                blah=squeeze(thegrid(:,:,1));
                oceans2pairs(:,1)=blah(I);
                blah=squeeze(thegrid(:,:,2));
                oceans2pairs(:,2)=blah(I);
                oceans2pairs2=zeros(size(oceans2pairs));
                for qqq=1:length(I);
                    oceans2pairs2(qqq,1)=find(lon==oceans2pairs(qqq,1));
                    oceans2pairs2(qqq,2)=find(lat==oceans2pairs(qqq,2));
                end
                blahout=ismember(usepairs,oceans2pairs2,'rows');
                usepairs2=usepairs(blahout==1,:);
                for qqq=1:length(usepairs2);
                    oceans3(usepairs2(qqq,1),usepairs2(qqq,2))=currindex;
                end
                f=f+1;
            end 
        else
            havelons(havelons<180)=havelons(havelons<180)+360;
            leftlon=min(havelons);
            rightlon=max(havelons);
            leftlonind=find(lon==leftlon);
            if rightlon>360;
                rightlonind=find(lon==rightlon-360);
                rightlonind=rightlonind+144;
            else
                rightlonind=find(lon==rightlon);
            end
            nlons=length(leftlonind:rightlonind);
            eachbox=nlons/boxesthisrow;
            theboxes=cell(1,boxesthisrow);
            d=1;
            for q=1:boxesthisrow;
                theboxes{q}=d:min([d+ceil(eachbox)-1 144]);
                d=d+ceil(eachbox);
            end
            for q=1:boxesthisrow;
                currindex=indices(f);
                worklons=theboxes{q};
                uselons=leftlonind+worklons(1)-1:leftlonind+worklons(end)-1;
                uselons(uselons>144)=uselons(uselons>144)-144;
                usepairs=zeros(length(uselons)*length(uselats),2);
                usepairs(:,1)=repmat(transpose(uselons),[length(uselats) 1]);
                blah=repmat(uselats,[length(uselons) 1]);
                usepairs(:,2)=blah(:);
                I=find(abs(oceans2-k)<0.5);
                oceans2pairs=zeros(length(I),2);
                blah=squeeze(thegrid(:,:,1));
                oceans2pairs(:,1)=blah(I);
                blah=squeeze(thegrid(:,:,2));
                oceans2pairs(:,2)=blah(I);
                oceans2pairs2=zeros(size(oceans2pairs));
                for qqq=1:length(I);
                    oceans2pairs2(qqq,1)=find(lon==oceans2pairs(qqq,1));
                    oceans2pairs2(qqq,2)=find(lat==oceans2pairs(qqq,2));
                end
                blahout=ismember(usepairs,oceans2pairs2,'rows');
                usepairs2=usepairs(blahout==1,:);
                for qqq=1:length(usepairs2);
                    oceans3(usepairs2(qqq,1),usepairs2(qqq,2))=currindex;
                end
                f=f+1;
            end
        end
    end
end

I=find(lf>lfthreshold);
oceans3(I)=0;

%figure;
%pcolor(transpose(oceans3));

% cleaning up grid boxes that are non-contiguous
oceans4=oceans3;
thegrid2=zeros(size(thegrid));
for k=1:144;
for j=1:96;
    thegrid2(k,j,1)=find(lon==thegrid(k,j,1));
    thegrid2(k,j,2)=find(lat==thegrid(k,j,2));
end
end

alldistinct=unique(oceans4);
I=length(alldistinct);
for k=2:I;
    J=find(oceans4==alldistinct(k));
    baselats=squeeze(thegrid2(:,:,2));
    havelats=baselats(J);
    uniquelats=unique(havelats);
    difflats=diff(uniquelats);
    K=find(difflats>10);
    if isempty(K)==0;
        toobiglats=uniquelats(K+1:end);
        L=find(ismember(baselats,toobiglats)==1);
        M=intersect(J,L);
        Q=find(floor(alldistinct(k))==floor(alldistinct));
        Q1=alldistinct(Q(end));
        oceans4(M)=(Q1+floor(Q1)+1)/2;
    end
end

I=find(lf>lfthreshold);
oceans4(I)=0;
%figure;
%pcolor(transpose(oceans4));

% absorbing grid boxes that are too small
oceans5=oceans4;
I=0;
blah1=squeeze(thegrid2(:,:,1));
blah2=squeeze(thegrid2(:,:,2));
while isempty(I)==0;
    alldistinct2=unique(oceans5);
    sizes=zeros(length(alldistinct2),1);
    for k=1:length(alldistinct2);
        I=find(oceans5==alldistinct2(k));
        sizes(k)=sum(gba(I));
    end
    sizes2=sizes(2:end);
    meanbox=mean(sizes2);
    I=find(sizes<0.5*meanbox);

    J=find(oceans5==alldistinct2(I(1))); % just taking the first one (that is not the
                                         % land mask) and pairing it with a bigger one
    totaldist=zeros(length(alldistinct2),1);
    meanlonJ=squeeze(mean(blah1(J)));
    meanlatJ=squeeze(mean(blah2(J)));
    for j=1:length(alldistinct2);
        K=find(oceans5==alldistinct2(j));
        meanlon=squeeze(mean(blah1(K)));
        meanlat=squeeze(mean(blah2(K)));
        totaldist(j)=sqrt((meanlon-meanlonJ)^2+(meanlat-meanlatJ)^2);
    end
    K=find(alldistinct2==0);
    totaldist(K)=0;
    K=find(floor(alldistinct2)~=floor(alldistinct2(I(1))));
    totaldist(K)=0;
    x=min(totaldist(totaldist>0));
    if isempty(x)==0;
        K=find(totaldist==x); % closest box to the J box that is too small
        wantspace=find(oceans5==alldistinct2(K));
        wantind=oceans5(wantspace(1));
        oceans5(J)=wantind;
    else
        I=I(2:end);
    end
end

%figure;
%pcolor(transpose(oceans5));

% splitting grid boxes that are too big
% absorbing grid boxes that are too small
oceans6=oceans5;
I=0;
blah1=squeeze(thegrid2(:,:,1));
blah2=squeeze(thegrid2(:,:,2));
while isempty(I)==0;
    alldistinct2=unique(oceans6);
    sizes=zeros(length(alldistinct2),1);
    for k=1:length(alldistinct2);
        I=find(oceans6==alldistinct2(k));
        sizes(k)=sum(gba(I));
    end
    sizes2=sizes(2:end);
    meanbox=mean(sizes2);
    I=find(sizes>1.5*meanbox);
    I=I(2:end); % omitting the land mask
    if isempty(I)==1;
        break;
    end

    J=find(oceans6==alldistinct2(I(1))); % just taking the first one (that is not the
                                         % land mask) and pairing it with a bigger one

    % finding the aspect ratio so I know how to split the box
    widths=zeros(1,96);
    heights=zeros(1,144);
    for j=1:96;
        indval=oceans6(J);
        K=find(blah2==j);
        L=intersect(J,K);
        if floor(indval)~=5;
            leftmost=min(blah1(L));
            rightmost=max(blah1(L));
            lonwts=gba(:,j);
            widths(j)=sum(lonwts(leftmost:rightmost));
        else
            K=find(blah2==j);
            L=intersect(J,K);
            blah3=blah1;
            blah3(blah3<72)=blah3(blah3<72)+144;
            leftmost=min(blah3(L));
            rightmost=max(blah3(L));
            lonwts=cat(1,gba(:,j),gba(:,j));
            widths(j)=sum(lonwts(leftmost:rightmost));
        end
    end
    for j=1:144;
        K=find(blah1==j);
        L=intersect(J,K);
        botmost=min(blah2(L));
        topmost=max(blah2(L));
        heights(j)=sum(latwts(botmost:topmost));
    end
    medwidth=median(widths(widths>0));
    medheight=median(heights(heights>0));
    aspectratio=medheight/medwidth;
    if aspectratio<1;
        x=find(widths~=0);
        midpoints=zeros(1,length(x));
        for j=1:length(x);
            indval=oceans6(J(1));
            K=find(blah2==x(j));
            L=intersect(J,K);
            if floor(indval)~=5;
                leftmost=min(blah1(L));
                rightmost=max(blah1(L));
                midpoints(j)=(leftmost+rightmost)/2;
            else
                K=find(blah2==x(j));
                L=intersect(J,K);
                blah3=blah1;
                blah3(blah3<72)=blah3(blah3<72)+144;
                leftmost=min(blah3(L));
                rightmost=max(blah3(L));
                blah=(leftmost+rightmost)/2;
                if blah>144;
                    blah=blah-144;
                end
                midpoints(j)=blah;
            end
        end
        splitlon=round(mean(midpoints));
        Q=unique(oceans6);
        Q1=find(Q>indval);
        indval2=(indval+Q(Q1(1)))/2;
        J1=J(find(blah1(J)<=splitlon));
        J2=J(find(blah1(J)>splitlon));
        oceans6(J1)=indval;
        oceans6(J2)=indval2;
    else
        x=find(heights~=0);
        midpoints=zeros(1,length(x));
        for j=1:length(x);
            indval=oceans6(J(1));
            K=find(blah1==x(j));
            L=intersect(J,K);
            botmost=min(blah2(L));
            topmost=max(blah2(L));
            midpoints(j)=(botmost+topmost)/2;
        end
        splitlat=round(mean(midpoints));
        Q=unique(oceans6);
        Q1=find(Q>indval);
        indval2=(indval+Q(Q1(1)))/2;
        J1=J(find(blah2(J)<=splitlat));
        J2=J(find(blah2(J)>splitlat));
        oceans6(J1)=indval;
        oceans6(J2)=indval2;
    end
end
        
figure;
pcolor(transpose(oceans6));
n=length(unique(oceans6));
title({['Oceans divided into ' num2str(n-1) ' boxes'],...
  ['(Original choice was ' num2str(nboxes) ' boxes)']},'FontSize',15);

I=find(isnan(oceans6)==1);
outputmatrix(I)=-1e30;
fid=fopen('oceans6.bin','wb');
fwrite(fid,oceans6,'single','b');
fclose(fid);


%% Building random sequences

% The frequency w is defined as 2*adjustscale/tmscl, where the 2 is for Nyquist, and
% adjustscale is how often you adjust
% For example, if you are interested in processes that take two weeks to resolve,
% if you adjust every day, w=1/7
% w is nondimensional - w=1 is Nyquist
w=2*(nyrs*365/nt)/tmscl;
m=n-1;
x=randn(m,nt);
[b,a]=butter(1,w);
x=filtfilt(b,a,x')';
[u,s,v]=svd(x*x');
z=diag(diag(s).^(-1/2))*u'*x;
% z'*z is the identity, i.e., no cross correlation and equal variance
I=max(max(z));
z=ampl*z/I;

figure;
plot(transpose(z));

%% Building output files

% Output will be a binary file
% All output will be contained in a single file that is 144*96 (number of grid boxes)
% by nt+5 (number of times + a few extras)
% For each box, there will be the following columns
% lon # (out of 144), lat # (out of 96), actual longitude, actual latitude, box number
% (out of m), sequence associated with that lat/lon box

Q=unique(oceans6);
outputmatrix=zeros(144*96,nt+5);
for k=1:144;
for j=1:96;
    outputmatrix((k-1)*96+j,1)=k;
    outputmatrix((k-1)*96+j,2)=j;
    outputmatrix((k-1)*96+j,3)=lon(k);
    outputmatrix((k-1)*96+j,4)=lat(j);
    Jval=oceans6(k,j);
    Jind=find(Jval==Q)-1;
    outputmatrix((k-1)*96+j,5)=Jind;
    if Jind==0;
        outputmatrix((k-1)*96+j,6:end)=0;
    else
        outputmatrix((k-1)*96+j,6:end)=z(Jind,:);
    end
end
end

I=find(isnan(outputmatrix)==1);
outputmatrix(I)=-1e30;
fid=fopen('sysidout.bin','wb');
fwrite(fid,outputmatrix,'single','b');
fclose(fid);
