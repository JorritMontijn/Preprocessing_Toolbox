function [ses,indKeepList] = doRecalcdFoF(ses,intSwitch,indKeepList,strType,dblFraction,dblSecWindowSize,dblNeuropilSubtractionFactor)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	if nargin < 4 || isempty(strType)
		strType = 'neuron';
	end
	if nargin < 3 || isempty(indKeepList)
		indKeepList = true(1,numel(ses.(strType)));
	end
	if nargin < 2 || isempty(intSwitch)
		intSwitch = 1;
	end
	if ~exist('dblFraction','var') || isempty(dblFraction),dblFraction=0.5;end
	if ~exist('dblSecWindowSize','var') || isempty(dblSecWindowSize),dblSecWindowSize=30;end
	if ~exist('dblNeuropilSubtractionFactor','var'),dblNeuropilSubtractionFactor=[];end
	
	if intSwitch == 1 %rectify non-positive values to 0
		for intNeuron=1:numel(ses.(strType))
			ses.(strType)(intNeuron).dFoF(~(ses.(strType)(intNeuron).dFoF > 0)) = 0;
		end
	elseif intSwitch == 2 %set outliers to nan
		for intNeuron=1:numel(ses.(strType));
			indexList = getOutliers(ses.(strType)(intNeuron).dFoF,5);
			ses.(strType)(intNeuron).dFoF(indexList) = nan;
		end
	elseif intSwitch == 3 %dF/F without neuropil subtraction
		dblSamplingFreq = ses.samplingFreq;
		for intNeuron=1:numel(ses.(strType));
			fprintf('Recalculating dF/F0 for neuron %d/%d [%s]\n',intNeuron,numel(ses.(strType)),getTime);
			%get F
			F = ses.(strType)(intNeuron).F;
			
			% calculate dFoF
			ses.(strType)(intNeuron).dFoF = calcdFoF(F,dblSamplingFreq,[],dblFraction,dblSecWindowSize);
		end
	elseif intSwitch == 4 %dF/F from neuropil annulus
		dblSamplingFreq = ses.samplingFreq;
		for intNeuron=1:numel(ses.(strType));
			
			%get F
			F = ses.(strType)(intNeuron).npF;
			
			% calculate dFoF
			ses.(strType)(intNeuron).dFoF = calcdFoF(F,dblSamplingFreq,[],dblFraction,dblSecWindowSize);
		end
	elseif intSwitch == 5 %dF/F with neuropil subtraction
		dblSamplingFreq = ses.samplingFreq;
		for intNeuron=1:numel(ses.(strType));
			
			%get F
			if isempty(dblNeuropilSubtractionFactor)
				dblFactor = corr(ses.(strType)(intNeuron).F',ses.(strType)(intNeuron).npF');
			else
				dblFactor = dblNeuropilSubtractionFactor;
			end
			F = ses.(strType)(intNeuron).F - dblFactor*ses.(strType)(intNeuron).npF;
			
			% calculate dFoF
			ses.(strType)(intNeuron).dFoF = calcdFoF(F,dblSamplingFreq,[],dblFraction,dblSecWindowSize);
		end
	elseif intSwitch == 6 %remove neurons
		if nargin > 2 %supplied list
			intAddCount = 0;
			for intNeuron=1:numel(ses.(strType));
				if indKeepList(intNeuron)
					intAddCount = intAddCount + 1;
					neuron(intAddCount) = ses.(strType)(intNeuron);
				end
			end
		else %with extreme dF/F values
			intAddCount = 0;
			for intNeuron=1:numel(ses.(strType));
				if any(ses.(strType)(intNeuron).dFoF > 10)
					indKeepList(intNeuron) = false;
				else
					intAddCount = intAddCount + 1;
					neuron(intAddCount) = ses.(strType)(intNeuron);
				end
			end
			fprintf('Removed %d neurons due to extreme dF/F0 values\n',sum(~indKeepList));
		end
		ses = rmfield(ses,'neuron');
		ses.(strType) = neuron;
	elseif intSwitch == 7 %use raw F
		for intNeuron=1:numel(ses.(strType));
			
			%get F
			ses.(strType)(intNeuron).dFoF = ses.(strType)(intNeuron).F;
		end
	elseif intSwitch == 8 %post-dF/F neuropil subtraction
		dblSamplingFreq = ses.samplingFreq;
		for intNeuron=1:numel(ses.(strType));
			fprintf('Recalculating dF/F0 for neuron %d/%d [%s]\n',intNeuron,numel(ses.(strType)),getTime);
			
			%get F
			vecSoma = calcdFoF(ses.(strType)(intNeuron).F,dblSamplingFreq,[],dblFraction,dblSecWindowSize);
			vecNeuropil = calcdFoF(ses.(strType)(intNeuron).npF,dblSamplingFreq,[],dblFraction,dblSecWindowSize);
			
			%get F
			if isempty(dblNeuropilSubtractionFactor)
				dblFactor = corr(vecSoma',vecNeuropil');
			else
				dblFactor = dblNeuropilSubtractionFactor;
			end
			
			% calculate dFoF
			ses.(strType)(intNeuron).dFoF = vecSoma - dblFactor*vecNeuropil;
		end
	end
end

