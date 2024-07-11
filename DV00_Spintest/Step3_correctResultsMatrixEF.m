cd('...\Results_SpintestMatrix');

load('spinTest_MeanClass_results.mat'); 

resultsEF = results;
k = 3;
n = 219; 

resultsEF(2,[1,3:7,9:18],1) = (results(2,[1,3:7,9:18],1)-k+1)./(n-k);
resultsEF([1,3:7,9:18],2,1) = (results([1,3:7,9:18],2,1)-k+1)./(n-k);

resultsEF(8,[1,3:7,9:18],1) = (results(8,[1,3:7,9:18],1)-k+1)./(n-k);
resultsEF([1,3:7,9:18],8,1) = (results([1,3:7,9:18],8,1)-k+1)./(n-k);

resultsEF(2,8,1) = sqrt(results(2,8,1)/(n*(k-1)));
resultsEF(8,2,1) = sqrt(results(2,8,1)/(n*(k-1)));

save('spinTest_resultsEF.mat', 'resultsEF'); 
writematrix(resultsEF(:,:,1), 'spintestResults_effectSizes.csv')
writematrix(resultsEF(:,:,2),'spintestResults_pValues.csv')


% FDR correction
addpath("...\_Tools"); 
writematrix(resultsEF(:,:,1),'effectSizes_MeanVariables.txt')
writematrix(results(:,:,1),'rawValues_MeanVariables.txt')

imagesc(resultsEF(:,:,2) < .05)
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(resultsEF(:,:,2), 0.05, 'dep', 'yes');

writematrix(adj_p, "spin_pFDRtrue_MeanVariables.txt");
writematrix(adj_p < 0.05, "spin_pFDRtrue_sig_MeanVariables.txt");

pHalf = tril(resultsEF(:,:,2),-1);
pHalfBool = pHalf ~= 0; 
pHalfVec = pHalf(pHalfBool); 
[hh, crit_ph, adj_ci_cvrgh, adj_ph] = fdr_bh(pHalfVec, 0.05, 'dep', 'yes');
pHalfCorr = zeros(size(pHalf)); 
pHalfCorr(pHalfBool) = adj_ph; 
pCompleteCorr = pHalfCorr + pHalfCorr.'; 
pCompleteCorr(1:(size(pCompleteCorr,1)+1):end) = max(pCompleteCorr,[] ,'all'); 
imagesc(pCompleteCorr)

writematrix(pCompleteCorr, "spin_pFDRlower_MeanVariables.txt");
writematrix(pCompleteCorr < 0.05, "spin_pFDRlower_sig_MeanVariables.txt");

pHalf = triu(resultsEF(:,:,2),1);
pHalfBool = pHalf ~= 0; 
pHalfVec = pHalf(pHalfBool); 
[hh, crit_ph, adj_ci_cvrgh, adj_ph] = fdr_bh(pHalfVec, 0.05, 'dep', 'yes');
pHalfCorr = zeros(size(pHalf)); 
pHalfCorr(pHalfBool) = adj_ph; 
pCompleteCorr = pHalfCorr + pHalfCorr.'; 
pCompleteCorr(1:(size(pCompleteCorr,1)+1):end) = max(pCompleteCorr,[] ,'all'); 
imagesc(pCompleteCorr)

writematrix(pCompleteCorr, "spin_pFDRupper_MeanVariables.txt");
writematrix(pCompleteCorr < 0.05, "spin_pFDRupper_sig_MeanVariables.txt");
