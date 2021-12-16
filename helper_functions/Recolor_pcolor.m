function [h_bar]=Recolor_pcolor(h_pcol,Levls,C,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs:
%%% h_pcol is handle to Pcolor plot
%%% Levls is vector of contour levels
%%% in addtion to that max(Levls) should be less than max(C(:))
%%% C is the matrix used in pcolor function 
%%% varargin- can have only one argument 'vert/horizontal' for the colorbar 
%%% A colormap of individual colors 
%%% where such as(size(colormap(),1)=length(Levls)+1 must be set before
%%% calling this function.
%%% Outputs:
%%% h_bar is the handle to colorbar
%%% example:
%%% X=1:10;Y=1:10;
%%% C=rand(10,10)*150;
%%% c=colormap(jet(7));
%%% L=[10 20 50 60 70 100];
%%% caxis([-2 100]);
%%% h_pcol=pcolor(X,Y,C); shading flat;
%%% Recolor_pcolor(h_pcol,L,C,'vert')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cdt=get(h_pcol,'CData'); cdtsave = cdt;
Lold=Levls;
clm=caxis;
L1=((1:length(Levls))*(diff(caxis)/(length(Levls)+1)))+min(caxis());
Levls=[Levls Inf];
for k=1:length(Levls)-1,
      cdt(C>=Levls(k)&C<Levls(k+1))=L1(k); 
end 
set(h_pcol,'CData',cdt);
caxis(clm);
if(~isempty(varargin))
[h_bar]=colorbar_levels(Lold,varargin{1});
end

