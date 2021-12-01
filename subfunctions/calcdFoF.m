function [vecdFoF,vecF,vecF0] = calcdFoF(vecF, dblSamplingFreq,intType,dblFraction,dblSecWindowSize)
	%calcdFoF Calculates dF/F
	%   Syntax: [vecdFoF] = calcdFoF(vecF, dblSamplingFreq, intType)
	%
	%	Version history:
	%	1.0 - August 4 2013 [by Jorrit Montijn]
	%	Created separate calcdFoF subfunction; also vectorized all
	%	calculations to increase speed, and transformed smoothing kernel
	%	computation to operate using built-in conv() function instead of
	%	the custom-made smooth1D() function. Also changed the smoothing
	%	window so it now depends on the frame acquisition speed 
	
	
	%check input, otherwise default to type 2
	if ~exist('intType','var') || isempty(intType),intType=2;end
	if ~exist('dblFraction','var') || isempty(dblFraction),dblFraction=0.5;end
	if ~exist('dblSecWindowSize','var') || isempty(dblSecWindowSize),dblSecWindowSize=30;end
	
	%perform calculation
	if intType==1 
		%% loop
		% smooth and baseline the data
		dblWindowSecs = min( [ (0.4 * (length(vecF)/dblSamplingFreq)) dblSecWindowSize] ); %number of seconds for F0 baselining
		intWindowFrames = round(dblSamplingFreq*dblWindowSecs) ; %number of frames for F0 baselining
		
		vecF = smooth1D(vecF', 3, 'gaussian') ;
		vecF0 = zeros(size(vecF)) ;
		% calculate F0 (baseline) per frame
		sortF = sort( vecF( 1:intWindowFrames ) ) ;
		vecF0(1:intWindowFrames) = mean( sortF(1:round(intWindowFrames/2)) ) ;
		for i=(intWindowFrames+1):length(vecF)
			sortF = sort( vecF( i-intWindowFrames:i ) ) ;
			vecF0(i) = mean( sortF(1:round(intWindowFrames/2)) ) ;
		end
		
		%get dF/F
		vecdFoF = (vecF-vecF0)./vecF0; %calculate dF/F by subtracting F0 trace from F trace and dividing by F0 trace
	elseif intType == 2
		%% new version by jorrit [vectorized]
		%make smoothing kernel
		intKernelSteps = min([5 (round(dblSamplingFreq/4)*2-1)]); %kernel size has to be odd; set to be approximately half a second [or maximum of 5 to avoid over-smoothing]
		intKernelStepSize = 2/(intKernelSteps-1); %required size per step to get correct number of steps
		vecKernel = normpdf(-1:intKernelStepSize:1,0,1); %get gaussian kernel
		vecKernel = vecKernel / sum(vecKernel); %set integral to 1
		
		% smooth F
		vecF = imfilt(vecF(:)',vecKernel);
		
		% calculate F0 window size
		dblWindowSecs = min( [ (0.4 * (length(vecF)/dblSamplingFreq)) dblSecWindowSize] ); %number of seconds for F0 baselining
		intWindowFrames = round(dblSamplingFreq*dblWindowSecs) ; %number of frames for F0 baselining
		
		%calculate F0
		vecSelectBase=1:intWindowFrames; %base vector from which to build full selection matrix
		matSelect=repmat(vecSelectBase,[floor(intWindowFrames/2) 1]); %create static (non-shifting window) selection matrix for first part of trace
		intSizeSecond = (length(vecF)-intWindowFrames); %calculate size of second part of selection matrix
		matSelect=[matSelect;repmat(vecSelectBase,[intSizeSecond 1])+repmat((1:intSizeSecond)',[1 intWindowFrames])]; %add incrementally increasing selection trace matrix (with shifting window) to static first part
		matSelect=[matSelect;repmat(vecSelectBase,[ceil(intWindowFrames/2) 1])]; %create static (non-shifting window) selection matrix for last part of trace
		matSortedWindowTraces=sort(vecF(matSelect),2,'ascend');%select F values from trace with matSelect; then sort F values (ascending) per trace
		vecF0=mean(matSortedWindowTraces(:,1:round(intWindowFrames*dblFraction)),2)'; %take first half sorted values, so that lowest 50% are selected; then calculate mean per trace for each time point
		
		%get dF/F
		vecdFoF = (vecF-vecF0)./vecF0; %calculate dF/F by subtracting F0 trace from F trace and dividing by F0 trace
	elseif intType == 3
		%% newer version by jorrit [exponential smoothing; based on Jia et al. (2011), nature protocols]
		%% do not use; oversmooths traces with exponential filter
		
		%pre-allocate
		vecFsm = zeros(size(vecF));
		vecF0 = zeros(size(vecF));
		vecdFoF = zeros(size(vecF));
		
		%constants
		intEndFrame = length(vecF);
		dblTau0 = 0.2;
		dblTau1 = 0.75;
		dblTau2 = 3; %F0 window size in seconds
		
		%smooth F
		for intFrame=1:intEndFrame
			t=intFrame/dblSamplingFreq;
			intStartFrame = max(round((t-(dblTau1/2))*dblSamplingFreq),1);
			intStopFrame = min(round((t+(dblTau1/2))*dblSamplingFreq),intEndFrame);
			dblNormFac = 1/(intStopFrame - intStartFrame + 1);
			vecFsm(intFrame) = dblNormFac/sum(vecF(intStartFrame:intStopFrame));
		end
		
		%calculate F0
		for intFrame=1:intEndFrame
			t=intFrame/dblSamplingFreq;
			intStartWindow = max(round((t-dblTau2)*dblSamplingFreq),1);
			intStopWindow = min(round(t*dblSamplingFreq),intEndFrame);
			vecF0(intFrame) = min(vecFsm(intStartWindow:intStopWindow));
		end
		
		%calculate dF/F0
		vecR = (vecF-vecF0)./vecF0;
		
		%apply exponential filtering
		intExpCutOff = round((dblTau0*5)*dblSamplingFreq);
		vecExponential = (1:intExpCutOff)/dblSamplingFreq;
		dblIntegral2 = sum(vecExponential);
		for intFrame=1:intEndFrame
			intStopWindow = min(intFrame+intExpCutOff-1,intEndFrame);
			intStartWindow = min(intFrame,intStopWindow-intExpCutOff+1);
			
			dblIntegral1 = sum(vecR(intStartWindow:intStopWindow).*vecExponential);
			
			vecdFoF(intFrame) = dblIntegral1/dblIntegral2;
		end
	end
end

