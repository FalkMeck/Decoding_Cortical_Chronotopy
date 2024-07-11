%% Identify the most consesual partition for partition coefficient and within-module degree z-score

function [repPartTrue, gammaTrue] = consensual_partition(conMat, gammaVec, parti)
% This functional calculation the best gamma (gammaTrue) for louvain algorithm consesus
% matrix and then uses this gamma to give out the most fitting partition to
% calculate the the partition coefficinet or within module degree z score.
%
% INPUT
% conMAT: undirected connectivity matrix
% gammaVec: a vector of gamma values the algorithm is tested on
% parti: number of partitions that are calcualted at the start

nConMat = size(conMat,1);

group_modules = zeros(nConMat,parti,numel(gammaVec));
consensusPart = zeros(nConMat,nConMat, numel(gammaVec));
repPart = zeros(nConMat,numel(gammaVec));
for g = 1:numel(gammaVec)
    for li = 1:parti
        group_modules(:,li,g) = community_louvain(conMat,gammaVec(g));
    end
    consensusPart(:,:,g) = agreement(group_modules(:,:,g),150);
    repPart(:,g) = consensus_und(consensusPart(:,:,g), 0.6,100);
end

repComMat = zeros(numel(gammaVec),numel(gammaVec));
for gi = 1:(numel(gammaVec)-1)
    gamma_y = gi+1:numel(gammaVec);
    for gy = gamma_y
        repComMat(gy,gi) = variation_of_information(repPart(:,gi),repPart(:,gy));
    end
end
repComMatFull = repComMat + repComMat';
repComDiag = ~logical(eye(size(repComMatFull)));

medianVI = zeros(4,numel(gammaVec));
for gi = 1:numel(gammaVec)
    medianVI(1,gi) = gammaVec(gi);
    medianVI(2,gi) = median(repComMatFull(gi,repComDiag(gi,:)));
    [~, medianVI(3,gi)] = modularity_und(consensusPart(:,:,gi));
    %medianVI(3,gi) = mean(repComMatFull(gi,repComDiag(gi,:)));
    medianVI(4,gi) = sum(repComMatFull(gi,repComDiag(gi,:)));
end

gammaTrue = medianVI(1,medianVI(2,:) == min(medianVI(2,:)));
if numel(gammaTrue) ~= 1
    disp(['We have to recalculate. gammaTrue = ', num2str(gammaTrue)]);
    gammaTrue = mean(gammaTrue);
    group_modulesTrue = zeros(nConMat,parti);
    for li = 1:parti
        group_modulesTrue(:,li) = community_louvain(conMat,gammaTrue);
    end
    consensusPartTrue = agreement(group_modulesTrue,150);
    repPartTrue = consensus_und(consensusPartTrue, 0.6,100);
else
    repPartTrue = repPart(:, medianVI(1,:) == gammaTrue);
end
end