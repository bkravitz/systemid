function [h_bar]=Recolor_contourf(hc,cmap,Levls,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs:
%%% hc is handle to contours
%%% Levls is vector of contour levels
%%% The contours must be created on levels [min(data) Levls] otherwise the
%%% program will give wrong results
%%% in addtion to that max(Levls) should be less than max(data)
%%% cmap is the colormap of individual colors 
%%% size(cmap,1)=length(Levls)+1
%%% varargin can have only one argument 'vert/horizontal' for the colorbar 
%%% Outputs:
%%% h_bar is the handle to colorbar
%%% example:
%%% X=1:10;Y=1:10;
%%% C=rand(10,10)*150;
%%% c=colormap(jet(7));
%%% L=[10 20 50 60 70 100];
%%% [cdd hc]=contourf(X,Y,C,[-5 L]);
%%% Recolor_contourf(hc,c,L,'vert');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=[-inf Levls inf];
if(length(hc)==1)
    Cld = get(hc, 'Children');
else
    Cld=hc;
end
for j=1:length(Cld)
    if strcmp(get(Cld(j), 'Type'), 'patch')
        Iso = get(Cld(j), 'CData');
        for k=1:length(L)-1,
            if(Iso>=L(k) && Iso<L(k+1))
                set(Cld(j), 'FaceColor', cmap(k,:));
                if(Iso<0), set(Cld(j), 'LineStyle',':'); end;
            end
        end        
    end
end
h_bar=[];
if(~isempty(varargin))
[h_bar]=colorbar_levels(Levls,varargin{1});
end


