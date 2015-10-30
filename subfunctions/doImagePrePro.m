function sRec = doImagePrePro(sRec)
	%run pre-processing
	% - smoothing settings
	% - save session data (x, y, z, t, ch dimensions)
	% - load primary registration image as the average of the first 100 imgs
	% - loop through images
	% - - load image
	% - - update raw average image
	% - - apply phase correction (and save some examples)
	% - - apply smoothing (and save some examples)
	% - - apply img registration (and save some examples)
	% - - update processed average image
	% - update flags for processing steps taken
	% save images
	% - raw average
	% - HQ raw average
	% - processed average
	% - HQ processed average
	% - overlay of processed averages
	% - HQ overlay of processed averages
	
	%allImages:
	%->updateRaw
	%phasecorr
	%remove saturated pixels
	%smooth
	%register
	%->updateProc
	%save in temp
	
	%after loop:
	%move ims from temp to proc
	%make overlays
	%make HQs
	%save average ims
	
	%settings
	matSmoothFilter = fspecial(sRec.sPS.strSmoothMethod, sRec.sPS.intSmoothKernelSize);
	
	%get data
	strSourcePath = [sRec.sMD.strSourceDir sRec.sMD.strImgSource sRec.strSession filesep sRec.sRawLib.strRecording filesep];
	strTargetPath = [sRec.sMD.strMasterDir sRec.sMD.strImgTarget sRec.strSession filesep sRec.sProcLib.strRecording filesep];
	intMaxT = sRec.sProcLib.t;
	intMaxZ = sRec.sProcLib.z;
	intSizeX = sRec.sProcLib.x;
	intSizeY = sRec.sProcLib.y;
	intMaxCh = sRec.sProcLib.ch;
	intLengthT = length(num2str(intMaxT));
	intImNum = intMaxCh * intMaxT * intMaxZ;
	intPhaseCorr = round(sRec.sProcLog.dblMustPhaseCorrect);
	
	%define image locations
	strTargetIm = ['t%0' num2str(intLengthT) 'd_ch%02d.tif'];
	strSourceIm = [sRec.sRawLib.strName '_t%0' num2str(intLengthT) 'd_ch%02d.tif'];
	if isfield(sRec.sRawLib,'strRefIm') && ischar(sRec.sRawLib.strRefIm) && ~isempty(sRec.sRawLib.strRefIm)
		strRefIm = sRec.sRawLib.strRefIm;
	else
		strRefIm = strSourceIm;
	end
	if isfield(sRec.sRawLib,'strRefPath') && ischar(sRec.sRawLib.strRefPath) && ~isempty(sRec.sRawLib.strRefPath)
		strRefPath = sRec.sRawLib.strRefPath;
	else
		strRefPath = strSourcePath;
	end
	fprintf('Reference registration image used: %s%s\n',strRefPath,strRefIm);
	
	%make directories if not present
	structDir = dir(strTargetPath);
	if numel(structDir) == 0;
		mkdir([sRec.sMD.strMasterDir sRec.sMD.strImgSource sRec.strSession filesep],sRec.sRawLib.strRecording);
	end
	structDir = dir([strTargetPath 'average' filesep]);
	if numel(structDir) == 0;
		mkdir(strTargetPath,'average');
	end
	structDir = dir([strTargetPath 'images' filesep]);
	if numel(structDir) == 0;
		mkdir(strTargetPath,'images');
	end
	
	%wait msg
	strIntImLength = num2str(length(num2str(intMaxT)));
	fprintf('Please wait for preprocessing of images to finish...\n\n');
	strFormat = ['Preprocessing... Now processed %0' strIntImLength 'd of %0' strIntImLength 'd\n'];
	
	%disp msg
	fprintf(strFormat,0,intMaxT);
	
	%pre-allocate variables
	intSubPixelDepth = sRec.sPS.intSubPixelDepth;
	boolDoRegistration = sRec.sPS.boolDoRegistration;
	boolDoRemSaturated = sRec.sPS.boolDoRemSaturated;
	if sRec.sPS.boolDoRegistration
		performance = nan(1,intMaxT);
		phaseCorr = nan(1,intMaxT);
		xReg = nan(1,intMaxT);
		yReg = nan(1,intMaxT);
	end
	intUseRefIms = 500;
	matRegRefAllCh = zeros(intSizeY,intSizeX);
	intBufferLength = 3;
	matSlidingBuffer = nan(intSizeY,intSizeX,intMaxCh,intBufferLength);
	
	for intCh=1:intMaxCh
		%pre-allocate structs
		imAverage.Ch(intCh).Raw = zeros(intSizeY,intSizeX);
		imAverage.Ch(intCh).RawHQ = zeros(intSizeY,intSizeX);
		imAverage.Ch(intCh).Proc = zeros(intSizeY,intSizeX);
		imAverage.Ch(intCh).ProcHQ = zeros(intSizeY,intSizeX);
		
		imTempAverageRaw = zeros(intSizeY,intSizeX,intMaxCh);
		imTempAverageProc = zeros(intSizeY,intSizeX,intMaxCh);
		
		
		%make registration reference from first 100 images
		matRegRefTemp = zeros(intSizeY,intSizeX);
		
		for intT = (round(intMaxT/2)-intUseRefIms):round(intMaxT/2)
			%load image
			strIm = sprintf(strRefIm,intT-1,intCh-1);
			
			%transform to 2D to save space
			imThis = imread([strRefPath strIm]);
			if ndims(imThis) == 3
				imThis = mean(imThis,3);
			end
			
			matRegRefTemp = matRegRefTemp + im2double(imThis);
		end
		matRegRefTemp = matRegRefTemp / intUseRefIms;
		matRegRefAllCh = matRegRefAllCh + matRegRefTemp;
	end
	matRegRefAllCh = matRegRefAllCh/intMaxCh;
	fftRegRef = fft2(matRegRefAllCh);
	
	%pre-load image buffer
	for intT=1:intBufferLength
		for intCh=1:intMaxCh
			strIm = sprintf(strSourceIm,intT-1,intCh-1);
			
			%transform to 2D to save space
			imThis = imread([strSourcePath strIm]);
			if ndims(imThis) == 3
				imThis = mean(imThis,3);
			end
			matIm2D = im2double(imThis);
			
			%apply phase correction
			if intPhaseCorr > 0 || intPhaseCorr < 0
				matIm2D = doPhaseCorrect(matIm2D,intPhaseCorr);
			end
			
			%remove saturated pixels
			if boolDoRemSaturated
				matIm2D = doRemSaturated(matIm2D);
			end
			
			%put into buffer
			matSlidingBuffer(:,:,intCh,intT) = matIm2D;
		end
	end
	
	%do actual image processing
	x=tic;
	for intT=1:intMaxT
		%set wait bar
		if mod(intT,100) == 0
			fprintf(strFormat,intT,intMaxT);
		end
		
		%get current image
		if intT == 1
			matThisIm = matSlidingBuffer(:,:,:,1);
		elseif intT == 2
			matThisIm = matSlidingBuffer(:,:,:,2);
		else
			matThisIm = matSlidingBuffer(:,:,:,end);
			matSlidingBuffer = circshift(matSlidingBuffer,[0 0 0 -1]); %shift buffer
		end
		
		%load raw images
		intBufferIm = min(intT,intMaxT-1);
		matIm = zeros(intSizeY,intSizeX,intMaxCh);
		for intCh=1:intMaxCh
			
			%get buffer image
			strBufferIm = sprintf(strSourceIm,intBufferIm,intCh-1);
			
			%transform to 2D to save space
			imBuffer = imread([strSourcePath strBufferIm]);
			if ndims(imBuffer) == 3
				imBuffer = mean(imBuffer,3);
			end
			matIm(:,:,intCh) = im2double(imBuffer);
			
			%update average raw image
			imTempAverageRaw(:,:,intCh) = imTempAverageRaw(:,:,intCh) + matIm(:,:,intCh);
			
			%apply phase correction
			if intPhaseCorr > 0 || intPhaseCorr < 0
				matIm(:,:,intCh) = doPhaseCorrect(matIm(:,:,intCh),intPhaseCorr);
			end
			
			%remove saturated pixels
			if boolDoRemSaturated
				matIm(:,:,intCh) = doRemSaturated(matIm(:,:,intCh));
			end
			
			%apply smoothing
			%if sRec.sPS.boolDoSmooth
			%	matIm{intCh} = doSmoothingPrePro(matIm{intCh},sRec);
			%end
		end
		
		%put into buffer
		if intT ~= 1, matSlidingBuffer(:,:,:,end) = matIm;end
		
		%make mean image from buffer for registration
		matRegIm = mean(mean(matSlidingBuffer,4),3);
		
		%apply image registration to multi-chan image
		if boolDoRegistration
			%calculate required registration
			output = dftregistration( fftRegRef, fft2(matRegIm), intSubPixelDepth);
			
			%put output into log
			performance(intT) = output(1);
			phaseCorr(intT) = output(2);
			xReg(intT) = output(3);
			yReg(intT) = output(4);
		end
		
		%apply correction to all channels
		matImProc = zeros(intSizeY,intSizeX,intMaxCh);
		for intCh=1:intMaxCh
			matImProc(:,:,intCh) = circshift(matThisIm(:,:,intCh),round([output(3) output(4)]));
			
			%save image
			strFilename = sprintf(strTargetIm,intT-1,intCh-1);
			imwrite(matImProc(:,:,intCh),[strTargetPath 'images' filesep strFilename],'tiff', 'Compression', 'lzw');
		end
		
		%update average processed image
		imTempAverageProc = imTempAverageProc + matImProc;
	end
	toc(x)
	for intCh=1:intMaxCh
		%put output into log
		sRec.sProcLog.registration.Ch(intCh).performance = performance;
		sRec.sProcLog.registration.Ch(intCh).phaseCorr = phaseCorr;
		sRec.sProcLog.registration.Ch(intCh).xReg = xReg;
		sRec.sProcLog.registration.Ch(intCh).yReg = yReg;
		
		%output average images
		imAverage.Ch(intCh).Raw = imTempAverageRaw(:,:,intCh) / intMaxT;
		imAverage.Ch(intCh).Proc = imTempAverageProc(:,:,intCh) / intMaxT;
		
		%make HQ
		IhqRaw = imAverage.Ch(intCh).Raw;
		backGroundRaw = imopen(IhqRaw, strel('disk', 30)) ;
		IhqRaw = imsubtract(IhqRaw, backGroundRaw) ;
		IhqRaw = imadjust(IhqRaw);
		imAverage.Ch(intCh).RawHQ = IhqRaw;
		
		IhqProc = imAverage.Ch(intCh).Proc;
		backGroundProc = imopen(IhqProc, strel('disk', 30)) ;
		IhqProc = imsubtract(IhqProc, backGroundProc) ;
		IhqProc = imadjust(IhqProc);
		imAverage.Ch(intCh).ProcHQ = IhqProc;
		
		%save images
		imwrite(imAverage.Ch(intCh).Raw,[strTargetPath sprintf('average%sRaw_ch%02d.tif',filesep,intCh-1)],'tiff', 'Compression', 'lzw');
		imwrite(imAverage.Ch(intCh).Proc,[strTargetPath sprintf('average%sProc_ch%02d.tif',filesep,intCh-1)],'tiff', 'Compression', 'lzw');
		imwrite(imAverage.Ch(intCh).RawHQ,[strTargetPath sprintf('average%sRawHQ_ch%02d.tif',filesep,intCh-1)],'tiff', 'Compression', 'lzw');
		imwrite(imAverage.Ch(intCh).ProcHQ,[strTargetPath sprintf('average%sProcHQ_ch%02d.tif',filesep,intCh-1)],'tiff', 'Compression', 'lzw');
		
		%put into overlay
		imAverage.Overlay.Raw(:,:,intCh) = imAverage.Ch(intCh).Raw;
		imAverage.Overlay.RawHQ(:,:,intCh) = imAverage.Ch(intCh).RawHQ;
		imAverage.Overlay.Proc(:,:,intCh) = imAverage.Ch(intCh).Proc;
		imAverage.Overlay.ProcHQ(:,:,intCh) = imAverage.Ch(intCh).ProcHQ;
	end
	%add dummy channel 3
	imAverage.Overlay.Raw(:,:,3) = zeros(size(imAverage.Ch(intCh).Raw(:,:,1)));
	imAverage.Overlay.RawHQ(:,:,3) = zeros(size(imAverage.Ch(intCh).RawHQ(:,:,1)));
	imAverage.Overlay.Proc(:,:,3) = zeros(size(imAverage.Ch(intCh).Proc(:,:,1)));
	imAverage.Overlay.ProcHQ(:,:,3) = zeros(size(imAverage.Ch(intCh).ProcHQ(:,:,1)));
	
	%save overlay images
	imwrite(imAverage.Overlay.Raw,[strTargetPath sprintf('average%sOverlayRaw.tif',filesep)],'tiff', 'Compression', 'lzw');
	imwrite(imAverage.Overlay.Proc,[strTargetPath sprintf('average%sOverlayProc.tif',filesep)],'tiff', 'Compression', 'lzw');
	imwrite(imAverage.Overlay.RawHQ,[strTargetPath sprintf('average%sOverlayRawHQ.tif',filesep)],'tiff', 'Compression', 'lzw');
	imwrite(imAverage.Overlay.ProcHQ,[strTargetPath sprintf('average%sOverlayProcHQ.tif',filesep)],'tiff', 'Compression', 'lzw');
	sRec.imAverage = imAverage;
	
	%save file
	strNewDir = [sRec.sMD.strMasterDir sRec.sMD.strImgTarget sRec.strSession filesep sRec.sProcLib.strRecording filesep];
	strFile = [strNewDir sRec.strSession sRec.sProcLib.strRecording '_prepro.mat'];
	save(strFile,'sRec');
	fprintf('Saved recording structure to %s\n',strFile)
end

