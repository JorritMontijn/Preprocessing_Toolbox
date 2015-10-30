function [dblPhaseCorr,output] = doCheckPhasePrePro(sRec)
	
	%% get data
	fprintf('Retrieving recording data...\n\n')
	strSourcePath = [sRec.sMD.strSourceDir sRec.sMD.strImgSource sRec.strSession filesep sRec.sRawLib.strRecording filesep];
	intMaxT = sRec.sProcLib.t;
	intMaxZ = sRec.sProcLib.z;
	intSizeX = sRec.sProcLib.x;
	intSizeY = sRec.sProcLib.y;
	intMaxCh = sRec.sProcLib.ch;
	intLengthT = length(num2str(intMaxT));
	intImCounter = 0;
	intUseT = min([1000 intLengthT]);
	intImNum = intMaxCh * intUseT * intMaxZ;
	
	%define image locations
	strSourceIm = [sRec.sRawLib.strName '_t%0' num2str(intLengthT) 'd_ch%02d.tif'];
	imAverageRaw = zeros(intSizeY,intSizeX);
	
	%% build raw average
	fprintf('Building average image. Please wait...\n\n')
	h=waitbar(0,'Building average image. Please wait...');
	imAverage = struct;
	for intCh=1:intMaxCh
		for intT=1:intUseT
			intImCounter = intImCounter + 1;
			waitbar(intImCounter/intImNum,h,'Building average image for phase correction. Please wait...')
			
			%load raw image
			strIm = sprintf(strSourceIm,intT-1,intCh-1);
			
			%transform to 2D to save space
			imThis = imread([strSourcePath strIm]);
			if ndims(imThis) == 3
				imThis = mean(imThis,3);
			end
			matIm = im2double(imThis);
			
			%update average raw image
			imAverageRaw = imAverageRaw + matIm;
		end
		%output average images
		imAverage.Ch(intCh).Raw = imAverageRaw / intMaxT;
		
		%make HQ
		IhqRaw = imAverage.Ch(intCh).Raw;
		backGroundRaw = imopen(IhqRaw, strel('disk', 30)) ;
		IhqRaw = imsubtract(IhqRaw, backGroundRaw) ;
		IhqRaw = imadjust(IhqRaw);
		imAverage.Ch(intCh).RawHQ = IhqRaw;
		
		%put in overlay
		imAverage.Overlay.Raw(:,:,intCh) = imAverage.Ch(intCh).Raw;
		imAverage.Overlay.RawHQ(:,:,intCh) = imAverage.Ch(intCh).RawHQ;
	end
	%close wait bar
	delete(h);
	drawnow;
	
	%add dummy channel 3
	imAverage.Overlay.Raw(:,:,3) = zeros(size(imAverage.Ch(intCh).Raw(:,:,1)));
	imAverage.Overlay.RawHQ(:,:,3) = zeros(size(imAverage.Ch(intCh).RawHQ(:,:,1)));
	
	%% calculate phase correction
	[dblPhaseCorr,output] = doPhaseRegistration(imAverage.Overlay.RawHQ);
end
