clc;
close all;
clear all;

load('dataTPcoverChange.mat');

diff(:,:,1) = (aft(:,:,1) - bef(:,:,1)).^2;
diff(:,:,2) = (aft(:,:,2) - bef(:,:,2)).^2;
diff(:,:,3) = (aft(:,:,3) - bef(:,:,3)).^2;

figure;
imagesc(bef(:,:,[2,1,3]));
title('before');
colorbar;
figure;
imagesc(aft(:,:,[2,1,3]));
title('after');
colorbar;
figure;
imagesc(diff(:,:,3));
title('difference');
colorbar;


NDVIbef = (bef(:,:,3)-bef(:,:,2))./(bef(:,:,3)+bef(:,:,2));
NDVIaft = (aft(:,:,3)-aft(:,:,2))./(aft(:,:,3)+aft(:,:,2));
NDVIdiff = (NDVIaft-NDVIbef).^2;
NDVIdiffnorm = (NDVIdiff - min(min(NDVIdiff))) ./ (max(max(NDVIdiff)) - min(min(NDVIdiff)));

hold on
figure;
imagesc(NDVIbef);
title('NDVIbef');
figure;
imagesc(NDVIaft);
title('NDVIaft');
figure;
imagesc(NDVIdiff);
title('NDVIdiff');
colorbar;

threshold = 0.5;
detectionDiff = double(diff(:,:,3) > threshold);
detectionNDVI = double(NDVIdiffnorm > threshold);

figure;
imagesc(detectionDiff);
figure;
imagesc(detectionNDVI);

hold on 
figure;
subplot(1,2,1);
imagesc(detectionDiff(:,:));
title('detectionDiff');
subplot(1,2,2);
imagesc(detectionNDVI(:,:));
title('detectionNDVI');

hold on
figure;
imagesc(GT(:,:));
[Kappa,OA,ConfMatrix,EstCOnfMatrix] = EstKappa(GT(GT<2),detectionNDVI(GT<2))

i = 1 ;
for th = 0.01 : 0.01 : 1
    detectionNDVI = double(NDVIdiffnorm > th);
    [Kappa,OA,ConfMatrix,EstCOnfMatrix] = EstKappa(GT(GT<2),detectionNDVI(GT<2));
    Kappas(i, :) = [th Kappa];
    i = i+1 ;
end
kpx = Kappas(:,1);
kpy = Kappas(:,2);
hold on;
figure;
plot(kpx,kpy);