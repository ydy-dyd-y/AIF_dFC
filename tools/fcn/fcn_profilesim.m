function [sq,pd] = fcn_profilesim(mat)
% fcn_profilesim
%   
%   [sq,pd] = fcn_profilesim(mat)
%
%   The community profile similarity measures the similarity of edge
%   communities at each node
% 
%   Inputs:
%       mat, 
%           edge communits in the nxn space
%   Outputs:
%       sp,
%           community similarity, squareform
%       pd,
%           community similarity, vector form
% 

hdl = @(x,y)(mean(x~=y,2));   % 두 행렬이 x,y로 입력될 경우, 각 행렬 행에서 같은 열 끼리 값을 비교 했을 때 다른 것들의 비율
pd = pdist(mat,hdl);   % pd벡터의 첫번째 요소 값은 mat행렬의 1열과 2열을 hdl(x, y)에 각각 넣었을 때 얻어지는 값
sq = 1 - squareform(pd);