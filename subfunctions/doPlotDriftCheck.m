%function doPlotDriftCheck(sRegZ)
	%doPlotDriftCheck Plots z-registration data
	%   Input: sRegZ (doDriftCheck output structure)
	
	%retrieve variables
	matRegistrationZ = sRegZ.matRegistrationZ;
	sData = sRegZ.sData;
	strSession = sRegZ.strSession;
	strRecording = sRegZ.strRecording;
	sRec = sRegZ.sRec;

	% plot similarity heat map
	matRegNorm = abs(imnorm(matRegistrationZ(:,:,1))-1);
	figure,imagesc(matRegNorm);
	ylabel('Frame number')
	xlabel('Z plane number')
	colorbar
	title([strSession strRecording '; Color: normalized similarity'])

	%calculate center of mass
	matDistZ = repmat(1:size(matRegNorm,2),[size(matRegNorm,1) 1]);
	vecMeanZ = mean(matRegNorm .* matDistZ,2);
	vecMeanZ = (vecMeanZ - mean(vecMeanZ)) * sData.dblMicronPerPlane;

	%plot center of mass
	figure
	plot((1:length(vecMeanZ)-1)*sRec.xml.sData.dblFrameDur,vecMeanZ(1:(end-1)))
	ylim([-10 10])
	ylabel('Z drift from mean (micron)')
	xlabel('Time (s)')
%end

