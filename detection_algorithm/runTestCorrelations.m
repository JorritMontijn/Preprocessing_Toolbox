%% get figure location
intCaCh = sRec.sProcLib.ch; %channel that contains calcium data
intMaxX = sRec.sProcLib.x;
intMaxY = sRec.sProcLib.y;

strPath = [sRec.sMD.strMasterDir sRec.sMD.strImgTarget sRec.strSession filesep sRec.sProcLib.strRecording];
sDir = dir(strPath);
if isempty(sDir),strPath = strRecPath;end
strImage = [strPath filesep 'average' filesep 'OverlayProcHQ.tif'];
imProc = im2double(imread(strImage));

matMean = squeeze(xmean(xmean(sRec.sPixResp.matPixelCorrMat,1),2));
intSpatWin = size(matMean,1);
intSubSW = floor(intSpatWin/2);
imProc = padarray(imProc,[intSpatWin*2 intSpatWin*2 0],0,'both');

%% build aggregate map to select potential regions of interest
matPixelCorrMat = sRec.sPixResp.matPixelCorrMat;
intSubY = size(sRec.sPixResp.matPixelCorrMat,1);
intSubX = size(sRec.sPixResp.matPixelCorrMat,2);
matROIs = zeros(size(sRec.sPixResp.matPixelCorrMat,1),size(sRec.sPixResp.matPixelCorrMat,2));

for intY=144:256
	intY
	for intX=156:256
		
		subplot(2,2,1)
		mS=reshape(sRec.sPixResp.matPixelCorrMat(intY,intX,:,:),[25 25])-matMean;imagesc(mS,[-0.05 0.05])
		imCorr = zeros(intSpatWin*2,intSpatWin*2);
		imCorr(1:2:end,1:2:end) = mS;
		imCorr(2:2:end,1:2:end) = mS;
		imCorr(1:2:end,2:2:end) = mS;
		imCorr(2:2:end,2:2:end) = mS;
		imshow(imCorr,[0 0.05]);
		title(sprintf('x %d; y %d',intX,intY))
		colorbar
		
		subplot(2,2,2)
		vecY = intSpatWin*2+intY*2+[-intSpatWin:intSpatWin];
		vecX = intSpatWin*2+intX*2+[-intSpatWin:intSpatWin];
		imshow(imProc(vecY,vecX,:));
		hold on
		scatter(intSpatWin+1,intSpatWin+2,'xr');
		hold off
		pause
		
		mS=reshape(sRec.sPixResp.matPixelCorrMat(intY,intX,:,:),[25 25])-matMean;
		vecY = intY+[-intSubSW:intSubSW];
		indUseY = vecY>0 & vecY<=intSubY;
		vecX = intX+[-intSubSW:intSubSW];
		indUseX = vecX>0 & vecX<=intSubX;
		matROIs(vecY(indUseY),vecX(indUseX)) = matROIs(vecY(indUseY),vecX(indUseX)) + double(mS(indUseY,indUseX));
		
	end
end

