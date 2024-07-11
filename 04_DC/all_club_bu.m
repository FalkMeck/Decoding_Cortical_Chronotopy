function [C,Nk,Ek] = all_club_bu(CIJ,SCORE)
%DIVERSE_CLUB_BU        Diverse club coefficients (binary undirected graph)
%
%   D = diverse_club_bu(CIJ,REPMAT)
%   [D,Nk,Ek] = diverse_club_bu(CIJ,REPMAT)
%
%   The rich club coefficient, R, at level k is the fraction of edges that
%   connect nodes of degree k or higher out of the maximum number of edges
%   that such nodes might share.
%
%   Input:      CIJ,       connection matrix, binary and undirected
%            SCORE,        a sequence of weights the nodes are sorted by,
%                          has to have the as many entries as nodes exist
%
%   Output:       D,        vector of cloubness coefficients for levels
%                              1 to klevel.
%                Nk,        number of nodes with participation coefficient >k
%                Ek,        number of edges remaining in subgraph with
%                              coeffcient >k
%
%   Reference: Berolero et al. (2017) Nat. comm. 8:1277.; Alstott et al.
%   (2014) Scientific Reports 4:7258


% if nargin == 2
%     klevel = numel(unique(PC));
% elseif nargin > 2
%     error('number of inputs incorrect. Should be [CIJ,REPMAT], or [CIJ,REPMAT, klevel]')
% end
if numel(SCORE) ~= size(CIJ,1)
    error('Dimensions do not match')
end

if floor(SCORE) == SCORE
    klevel = max(SCORE);
    uniScore = 1:klevel;
else
    uniScore = unique(SCORE);
    klevel = numel(uniScore);
end

C = zeros(1,klevel);
Nk = zeros(1,klevel);
Ek = zeros(1,klevel);


for k = 1:klevel
    SmallNodes=find(SCORE<=uniScore(k));       %get 'small nodes' with degree <=k
    subCIJ=CIJ;                       %extract subnetwork of nodes >k by removing nodes <=k of CIJ
    subCIJ(SmallNodes,:)=[];          %remove rows
    subCIJ(:,SmallNodes)=[];          %remove columns
    Nk(k)=size(subCIJ,2);             %number of nodes with degree >k
    Ek(k)=sum(subCIJ(:));             %total number of connections in subgraph
    C(k)=Ek(k)/(Nk(k)*(Nk(k)-1));     %unweighted rich-club coefficient
end