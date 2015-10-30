function doShowF0(ses,intNeuron)
	[vecdFoF,vecF,vecF0] = calcdFoF(ses.neuron(intNeuron).F, ses.samplingFreq,2);
	intOffset = 2;
	intStart = 8000;
	intStop = 15000;
	dblSmoothSigma = 4;
	
	vecX = ((intStart+intOffset):(intStop-intOffset)) / ses.samplingFreq;
	vecF = conv(vecF,normpdf(-(dblSmoothSigma*3):(dblSmoothSigma*3),0,dblSmoothSigma),'same');
	vecPlotF = vecF((intStart+intOffset):(intStop-intOffset));
	subplot(2,2,1)
	plot(vecX,vecPlotF)
	ylabel('Fluorescence')
	xlabel('Time (s)')
	
	subplot(2,2,2)
	vecF0 = conv(vecF0,normpdf(-(dblSmoothSigma*3):(dblSmoothSigma*3),0,dblSmoothSigma),'same');
	vecPlotF0 = vecF0((intStart+intOffset):(intStop-intOffset));
	plot(vecX,vecPlotF0);
	ylabel('Sliding Baseline Fluorescence')
	xlabel('Time (s)')
	
	subplot(2,2,3)
	vecdFoF = conv(vecdFoF,normpdf(-(dblSmoothSigma*3):(dblSmoothSigma*3),0,dblSmoothSigma),'same');
	vecPlotdFoF = vecdFoF((intStart+intOffset):(intStop-intOffset));
	plot(vecX,vecPlotdFoF)
	ylabel('dF/F0')
	xlabel('Time (s)')
	
	subplot(2,2,4)
	plot(vecX,vecPlotF-vecPlotdFoF)
	ylabel('F - dF/F0 (noise removed)')
	xlabel('Time (s)')
end