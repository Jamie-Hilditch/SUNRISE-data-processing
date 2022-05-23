function [X, Y, M, SD, N] = binave2(x,y,binsize,xmin,xmax)
% [X, Y, M, SD, SE, SDM, SEM] = binave(x,y,binsize,xmin,xmax)
% This function takes the data in y and sorts it 
% into bins according to x where:
%
% xmin is the smallest anticipated value of x
% xmax is the largest anticipated value of x
% binsize is the size of each bin, i.e.,
% bin 1 = xmin to xmin+binsize
% bin 2 = xmin+binsize to xmin+2*binsize
% bin 3 = xmin+2*binsize to xmin+3*binsize
% etc.
%
% The outputs are:
% X: average value of x in each bin
% Y: average value of y in each bin
% M: median value of y in each bin
% SD: standard deviation about mean of y in each bin
% SE: standard error about mean of y in each bin
% SDM: standard deviation about median of y in each bin
% SE: standard error about median of y in each bin
%


X = [];
Y = [];
M = [];
SD = [];
SE = [];
SDM = [];
SEM = [];
N = [];
bin = 0;
for i = xmin:binsize:xmax-binsize
   j = find(x>i & x<i+binsize & ~isnan(x) );
   if ~isempty(j)
      bin = bin + 1;
      X = [X; mean(x(j))];
      if length(j)>1
          M = [M; nanmedian(y(:,j)')];
          Y = [Y; nanmean(y(:,j)')];
          SD = [SD; nanstd(y(:,j)')];
          N = [N; length(j)];
      else
          M = [M; y(:,j)'];
          Y = [Y; y(:,j)'];
         SD = [SD; 0*y(:,j)'];
         N = [N; length(j)];
      end;
%       SE = [SE; std(y(j),1)/sqrt(length(j))];
%       SDM = [SDM; sqrt(mean((y(j)-median(y(j))).^2))];
%       SEM = [SEM; sqrt(mean((y(j)-median(y(j))).^2))/sqrt(length(j))];
   else
      % X = [X; NaN];
   end
end

