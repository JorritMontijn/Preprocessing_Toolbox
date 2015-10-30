function spikeData = doDetectCalciumTransients(sRec,boolDetectSpikes,dblTau)
	% script that detects spikes in dFoF data
	
	%check input
	if nargin < 2,boolDetectSpikes = true;end %set switch to detect spikes
	if nargin < 3,dblTau = 0.500;end % set shape of exponential
	
	%generate name
	strRec = sprintf('%sxyt%02d',sRec.strSession,sRec.intRecording);
	
	%calculate actual duration of experiment; saved in format: 10m12.01s
	strT = sRec.xml.sData.strActualImageSizeT;
	intLocM = strfind(strT,'m');
	intMins = str2double(strT(1:(intLocM-1)));
	dblTotDurSecs = str2double(strT((intLocM+1):(end-1))) + intMins*60;
	dblFrameTime = dblTotDurSecs/sRec.sProcLib.t;
	dblSamplingFreq = 1/dblFrameTime;
	
	% calculate activity traces for all neurons
	fprintf('Calculating dFoF and transients for [%s], please wait...\n',strRec);
	for intObject = 1:length(sRec.timeseries.roi)
		strType = sRec.sDC.metaData.cellType{sRec.sDC.ROI(intObject).intType};
		vecC=clock;
		fprintf('\nNow at object %d of %d [type: %s; time: %02.0f:%02.0f:%02.0f]\n',intObject,length(sRec.timeseries.roi),strType,vecC(4),vecC(5),vecC(6));
		if ismember(sRec.sDC.ROI(intObject).intType,sRec.sDC.metaData.vecNeurons) %if is neuron
			%msg
			ptrTime=tic;
			fprintf('	Calculating dF/F0...\n')
			
			%get fluorescence without neuropil subtraction
			F = sRec.timeseries.roi(intObject).F;
			
			% calculate dF/F0
			dFoF = calcdFoF( F, dblSamplingFreq ) ;
			fprintf('\b	Done! Took %.1f seconds\n',toc(ptrTime))
			
			if boolDetectSpikes
				%msg
				ptrTime=tic;
				fprintf('	Performing spike detection...\n')
				
				%detect spikes
				[apFrames, apSpikes, vecSpikes, expFit] = doDetectSpikes(dFoF,dblSamplingFreq,dblTau);
				fprintf('\b	Done! Took %.1f seconds; %d transients dectected; mean spiking rate is %.2f Hz\n',toc(ptrTime),sum(apSpikes),sum(apSpikes)/dblTotDurSecs)
			else
				expFit = [];
				apFrames = [];
				apSpikes = [];
				vecSpikes = [];
			end
			spikeData(intObject).dFoF = dFoF;
			spikeData(intObject).expFit = expFit;
			spikeData(intObject).apFrames = apFrames;
			spikeData(intObject).apSpikes = apSpikes;
			spikeData(intObject).vecSpikes = vecSpikes;
		end
	end
	fprintf('\ndFoF and transient detection completed\n')
end



function display_F_and_spikes( dFoF, samplingFreq, apFrames, apSpikes, expFit, numRows )
	% Output the dFoF and spikes that survived the selection
	figure;
	xvalues = (1:length(dFoF))/samplingFreq;
	for r = 1:numRows
		subplot(numRows,1,r)
		datachop = length(dFoF)/numRows;
		datarange = round((r-1)*datachop)+1:round(r*datachop);
		
		plot(xvalues(datarange), dFoF(datarange));
		hold on;
		plot(xvalues(datarange), expFit(datarange), 'color', [0 1 0]);
		
		for f= 1:length(apFrames)
			if apFrames(f) >= datarange(1) && ...
					apFrames(f) <= datarange(end)
				line( [xvalues(apFrames(f)) xvalues(apFrames(f))], ...
					[-(max(dFoF)-min(dFoF))/10 0], 'color', [1 0 0], ...
					'linewidth', apSpikes(f) );
			end
		end
		set(gca, 'XLim', [xvalues(datarange(1)) xvalues(datarange(end))]);
		set(gca, 'YLim', [min([min(dFoF) min(expFit)]) ...
			max([max(dFoF) max(expFit)])]);
	end
end