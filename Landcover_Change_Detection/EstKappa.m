function [Kappa,OA,ConfMatrix,EstConfMatrix] = EstKappa(GT,map)

% EstKappa: estimates confusion matrix and OA / Kappa statistic for binary
% problems
%
% usage: [Kappa,OA,ConfMatrix,EstCOnfMatrix] = EstKappa(GT,map)
%
% inputs :      - GT = binary ground control map (n x 1)
%               - map = binary detection result (n x 1)
%
% outputs:      - Kappa = Kappa statistic (scalar)
%               - OA = overall accuracy (scalar)
%               - ConfMatrix = observed confusion matrix (2 x 2)
%               - EstConfMatrix = estimated confusion matrix (2 x 2)
%
% devis dot tuia at epfl dot ch 2011

if size (GT,1) ~= size (map,1)
    disp('The two datasets have different size. Aborting');
    return
end

classes = unique(GT);
tot = length(GT);


ConfMatrix(1,1) = length(find(GT + map == 0));
ConfMatrix(2,1) = length(find(GT - map == -1));
ConfMatrix(1,2) = length(find(GT - map == 1));
ConfMatrix(2,2) = length(find(GT + map == 2));

OA = sum(diag(ConfMatrix))/tot*100;

EstConfMatrix(1,1) = sum(ConfMatrix(:,1))*sum(ConfMatrix(1,:))/tot;
EstConfMatrix(2,1) = sum(ConfMatrix(:,2))*sum(ConfMatrix(1,:))/tot;
EstConfMatrix(1,2) = sum(ConfMatrix(:,1))*sum(ConfMatrix(2,:))/tot;
EstConfMatrix(2,2) = sum(ConfMatrix(:,2))*sum(ConfMatrix(2,:))/tot;

Kappa = (OA/100 - sum(diag(EstConfMatrix))/tot)./(1-sum(diag(EstConfMatrix))/tot);
