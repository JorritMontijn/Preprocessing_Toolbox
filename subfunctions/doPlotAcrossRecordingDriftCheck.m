function doPlotAcrossRecordingDriftCheck(sRegZ)
	%doPlotAcrossRecordingDriftCheck Plots z-registration data
	%   Input: sRegZ (doAcrossRecordingDriftCheck output structure)
	
	%retrieve variables
	matRegistrationZ = sRegZ.matRegistrationZ(:,:,:,1);
	sData = sRegZ.sData;
	strSession = sRegZ.strSession;
	vecRecordings = sRegZ.vecRecordings;

	%reshape
	intRecs = size(matRegistrationZ,1);
	intT = size(matRegistrationZ,2);
	matRegNorm = nan(intRecs*intT,size(matRegistrationZ,3));
	for intRec=1:intRecs
		for intPointT=1:intT
			matRegNorm((intT*(intRec-1)+intPointT),:) = abs(imnorm(squeeze(matRegistrationZ(intRec,intPointT,:)))-1);
		end
	end
	
	%fit with gaussians
	vecMeans = nan(1,size(matRegNorm,1));
	matRegFit = nan(size(matRegNorm));
	for intPointT=1:size(matRegNorm,1)
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
	% plot similarity heat map
	figure;
	subplot(2,2,1);
	imagesc(rot90(matRegNorm));
	xlabel('Recording')
	vecTicksRec = round((intT-1)/2+1):intT:(intRecs*intT);
	set(gca,'XTick',vecTicksRec,'XTickLabel',vecRecordings)
	ylabel('Z plane number')
	colormap('hot');freezeColors;
	colorbar;cbfreeze;
	title([strSession '; Color: normalized similarity'])
	
	subplot(2,2,2);
	imagesc(rot90(matRegFit));
	xlabel('Recording')
	vecTicksRec = round((intT-1)/2+1):intT:(intRecs*intT);
	set(gca,'XTick',vecTicksRec,'XTickLabel',vecRecordings)
	ylabel('Z plane')
	set(gca,'clim',[0 1])
	colormap('hot');freezeColors;
	colorbar;cbfreeze;
	title('Gaussian fits; Color: normalized similarity')
	
	subplot(2,2,3);
	imagesc(rot90(matRegNorm-matRegFit),[-1 1]);%*abs(max(max(matRegNorm-matRegFit))));
	xlabel('Recording')
	vecTicksRec = round((intT-1)/2+1):intT:(intRecs*intT);
	set(gca,'XTick',vecTicksRec,'XTickLabel',vecRecordings)
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
	xlabel('Recording')
	set(gca,'XTick',vecTicksRec,'XTickLabel',vecRecordings)
	title('Recording stability')
end

