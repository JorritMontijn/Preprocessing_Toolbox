function vecDriftZ = doPlotWithinRecordingDriftCheck(sRegZ,boolDoPlots)
	%doPlotAcrossRecordingDriftCheck Plots z-registration data
	%   Input: sRegZ (doAcrossRecordingDriftCheck output structure)
	
	%retrieve variables
	if nargin < 2 || isempty(boolDoPlots),boolDoPlots=true;end
	matRegistrationZ = sRegZ.matRegistrationZ;
	sData = sRegZ.sData;
	strSession = sRegZ.strSession;
	
	%reshape
	intT = size(matRegistrationZ,1);
	matRegNorm = abs(imnorm(matRegistrationZ(:,:,1),1)-1);
	
	%fit with gaussians
	vecMeans = nan(1,intT);
	matRegFit = nan(size(matRegNorm));
	for intPointT=1:intT
		if mod(intPointT,100)==0 && boolDoPlots,fprintf('Fitting time point %d/%d [%s]\n',intPointT,intT,getTime);end
		mu=size(matRegNorm,2)/2;
		sigma=3;
		peak=1;
		baseline=0;
		vecParamsInit = [mu sigma peak baseline];
		vecY = matRegNorm(intPointT,:);
		[p1, mse] = MLFit('singleGaussian', vecParamsInit, 1:length(vecY), vecY);
		gaussVector = singleGaussian(1:length(vecY),p1);
		matRegFit(intPointT,:) = gaussVector;
		vecMeans(intPointT) = p1(1);
	end
	if boolDoPlots
		% plot similarity heat map
		figure;
		subplot(2,2,1);
		imagesc(rot90(matRegNorm));
		xlabel('Time (frame number)');
		ylabel('Z plane number')
		colormap('hot');freezeColors;
		colorbar;cbfreeze;
		title([strSession '; Color: normalized similarity'])
		
		subplot(2,2,2);
		imagesc(rot90(matRegFit));
		xlabel('Time (frame number)');
		ylabel('Z plane')
		set(gca,'clim',[0 1])
		colormap('hot');freezeColors;
		colorbar;cbfreeze;
		title('Gaussian fits; Color: normalized similarity')
		
		subplot(2,2,3);
		imagesc(rot90(matRegNorm-matRegFit),[-1 1]);%*abs(max(max(matRegNorm-matRegFit))));
		xlabel('Time (frame number)');
		ylabel('Z plane number')
		%set(gca,'clim',[0 1])
		colormap('redblue');freezeColors;
		colorbar;cbfreeze;
		title('Residuals')
		
		%plot mean of gaussian fit
		subplot(2,2,4);
		plot((vecMeans-mean(vecMeans(:))) * sData.dblMicronPerPlane)
		grid on
		ylim([-10 10])
		ylabel('Mean Z drift over recordings (micron)')
		xlabel('Time (frame number)');
		title('Recording stability');
	end
	
	%output
	vecDriftZ = (vecMeans-mean(vecMeans(:))) * sData.dblMicronPerPlane;
end

