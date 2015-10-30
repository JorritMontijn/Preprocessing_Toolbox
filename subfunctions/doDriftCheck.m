function sRegZ = doDriftCheck(sSourceData)
	%doDriftCheck Performs z-stack registration of a time series recording
	%   Outputs a similarity measure per z plane per time point
	%
	%Syntax:
	%sRegZ = doDriftCheck(sSourceData)
	%
	%Input must be a structure with the following fields:
	%- strMasterPath; path of master directory (session)
	%- strDirRawStackZ; path containing z-stack images
	%- strSession; name of session (date)
	%- strRecording; name of recording (xyt01 ... xyt99)
	
	%get source data
	strMasterPath = sSourceData.strMasterPath;
	strDirRawStackZ = sSourceData.strDirRawStackZ;
	strSession = sSourceData.strSession;
	strRecording = sSourceData.strRecording;
	
	%calc paths
	strOldPath = cd();
	strRecPath = [strMasterPath strSession filesep strRecording];
	vecFileSeps = strfind(strDirRawStackZ(1:(end-1)),filesep);
	strDirZ = strtok(strDirRawStackZ((vecFileSeps(end)+1):end),filesep);
	
	%parameters
	intUseChannelSwitch = 2; %1 only red; 2 both
	intBufferLength = 3;
	intSubPixelDepth = 1;
	
	%msg
	fprintf('Loading data and preparing Z-stack registration...\n')
	
	%% load prepro file
	load([strRecPath filesep strSession strRecording '_prepro.mat']);
	
	%% prepare z-stack
	cd(strDirRawStackZ);
	sDirZ = dir('*z*_ch*.tif');
	if intUseChannelSwitch == 2
		matStackZ = nan(sRec.sProcLib.y,sRec.sProcLib.x,numel(sDirZ)/2,2);
	else
		matStackZ = nan(sRec.sProcLib.y,sRec.sProcLib.x,numel(sDirZ));
	end
	
	%% load z-stack
	%get raw data
	for intFile=1:numel(sDirZ)
		%get info
		strFile = sDirZ(intFile).name;
		intZ = str2double(getFlankedBy(strFile,'_z','_ch'))+1;
		if intUseChannelSwitch == 2
			intCh = str2double(getFlankedBy(strFile,'_ch','.tif'))+1;
		else
			intCh = 1;
		end
		
		%load image & put in stack
		imTemp = imread(strFile);
		matStackZ(:,:,intZ,intCh) = im2double(imTemp);
	end
	matStackZ = mean(matStackZ,4);
	intSizeStackZ = size(matStackZ,3);
	
	%make fourier-transformed stack
	matStackZ_fft = nan(size(matStackZ));
	for intZ=1:intSizeStackZ
		matStackZ_fft(:,:,intZ) = fft2(matStackZ(:,:,intZ));
	end
	
	%get z-stack slice interval from xml file
	strStackZXML = [strDirRawStackZ filesep 'MetaData' filesep strDirZ '_Properties.xml'];
	sData = loadXMLPrePro(strStackZXML);
	sData.dblMicronPerPlane = abs(sData.dblActualImageSizeZ)/intSizeStackZ;
	
	%% define z-registration variables
	%get data
	intMaxT = sRec.sProcLib.t;
	intMaxZ = sRec.sProcLib.z;
	intSizeX = sRec.sProcLib.x;
	intSizeY = sRec.sProcLib.y;
	intMaxCh = sRec.sProcLib.ch;
	intMinCh = 1;
	intLengthT = length(num2str(intMaxT));
	if intUseChannelSwitch == 1
		intMaxCh = 1;
	end
	
	%define image locations
	strTargetIm = ['t%0' num2str(intLengthT) 'd_ch%02d.tif'];
	
	%pre-allocate
	matSlidingBuffer = nan(intSizeY,intSizeX,intMaxCh,intBufferLength);
	matRegistrationZ = nan(intMaxT,intSizeStackZ,4);
	
	%% loop through images & compare to z-stack
	%pre-load image buffer
	for intT=1:intBufferLength
		for intCh=intMinCh:intMaxCh
			strIm = sprintf(strTargetIm,intT-1,intCh-1);
			
			%transform to 2D to save space
			imThis = imread([strRecPath filesep 'images' filesep strIm]);
			if ndims(imThis) == 3
				imThis = mean(imThis,3);
			end
			matIm2D = im2double(imThis);
			
			%put into buffer
			matSlidingBuffer(:,:,intCh,intT) = matIm2D;
		end
	end
	
	%wait msg
	strIntImLength = num2str(length(num2str(intMaxT)));
	fprintf('Finished z-stack initialization; Please wait for registration of images to finish...\n\n');
	strFormat = ['Now processed %0' strIntImLength 'd of %0' strIntImLength 'd; time is [%s]\n'];
	
	%disp msg
	fprintf(strFormat,0,intMaxT,getTime);
	
	%do actual image processing
	ptrTimer=tic;
	for intT=1:intMaxT
		%set wait bar
		if mod(intT,100) == 0
			fprintf(strFormat,intT,intMaxT,getTime);
			
			dblTperSec = intT/toc(ptrTimer);%time points per second
			dblTotDurSecs = intMaxT / dblTperSec;
			dblRemaining = dblTotDurSecs - toc(ptrTimer);
			
			fprintf('\b; time remaining is [%.0fs] (%.2f time points / second)\n',dblRemaining,dblTperSec);
		end
		
		%shift buffer
		if intT > 2
			matSlidingBuffer = circshift(matSlidingBuffer,[0 0 0 -1]);
		end
		
		%load raw images
		intBufferIm = min(intT,intMaxT-1);
		matIm = zeros(intSizeY,intSizeX,intMaxCh);
		for intCh=intMinCh:intMaxCh
			
			%get buffer image
			strBufferIm = sprintf(strTargetIm,intBufferIm,intCh-1);
			
			%transform to 2D to save space
			imBuffer = imread([strRecPath filesep 'images' filesep strBufferIm]);
			if ndims(imBuffer) == 3
				imBuffer = mean(imBuffer,3);
			end
			matIm(:,:,intCh) = im2double(imBuffer);
		end
		
		%put into buffer
		if intT ~= 1, matSlidingBuffer(:,:,:,end) = matIm;end
		
		%make mean image from buffer for registration
		matRegIm = mean(mean(matSlidingBuffer,4),3);
		matReg_fft = fft2(matRegIm);
		
		%loop through z-stack to register images
		for intZ=1:intSizeStackZ
			%calculate required registration
			matRegistrationZ(intT,intZ,:) = dftregistration(matStackZ_fft(:,:,intZ), matReg_fft, intSubPixelDepth);
		end
	end
	toc(ptrTimer)
	
	%make structure
	sRegZ = struct;
	sRegZ.matRegistrationZ = matRegistrationZ;
	sRegZ.sData = sData;
	sRegZ.strSession = strSession;
	sRegZ.strRecording = strRecording;
	sRegZ.sRec = sRec;
	
	%save data
	strFile = [strRecPath filesep strSession strRecording '_zreg.mat'];
	save(strFile,'matRegistrationZ','sData','strSession','strRecording','sRec');
	fprintf('Saved recording structure to %s; time is [%s]\n',strFile,getTime)
	cd(strOldPath);
end

