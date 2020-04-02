function [matTranslationsX,matTranslationsY] = doInterpolateFromList(matDataPoints,matXi,matYi)

	


vecShiftsX = matDataPoints(:,3);
vecShiftsX(getOutliers(vecShiftsX)) = 0;
[dummy,dummy,matVx]=idx2mat(matDataPoints(:,1),matDataPoints(:,2),vecShiftsX);
vecShiftsY = matDataPoints(:,4);
vecShiftsY(getOutliers(vecShiftsY)) = 0;
[dummy,dummy,matVy]=idx2mat(matDataPoints(:,1),matDataPoints(:,2),vecShiftsY);

%extend matrices
vecPointsX = unique(matDataPoints(:,1));
vecPointsY = unique(matDataPoints(:,2));

vecMeshX = [min(matXi(:));vecPointsX(:);max(matXi(:))]';
vecMeshY = [min(matYi(:));vecPointsY(:);max(matYi(:))]';
[matX,matY]= meshgrid(vecMeshX,vecMeshY);

%extend matrixces
matVx = [matVx(1,1) matVx(1,:) matVx(1,end);...
	matVx(:,1) matVx matVx(:,end);...
	matVx(end,1) matVx(end,:) matVx(end,end)];
matVy = [matVy(1,1) matVy(1,:) matVy(1,end);...
	matVy(:,1) matVy matVy(:,end);...
	matVy(end,1) matVy(end,:) matVy(end,end)];

strMethod='cubic';
matTranslationsX = interp2(matX,matY,matVx,matXi,matYi,strMethod);
matTranslationsY = interp2(matX,matY,matVy,matXi,matYi,strMethod);

%{
figure,imagesc(matTranslationsX),colorbar
title('x shift')

figure,imagesc(matTranslationsY),colorbar
title('y shift')

figure,imagesc(sqrt(matTranslationsY.*matTranslationsY + matTranslationsX.*matTranslationsX)),colorbar
title('tot shift')
%}
end