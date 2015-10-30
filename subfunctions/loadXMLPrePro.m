function [sData,intOutFlag] = loadXMLPrePro(strXMLfile)
	
	%% pre-allocate variables
	%actual size
	sData.dblActualImageSizeX = [];
	sData.dblActualImageSizeY = [];
	sData.dblActualImageSizeZ = [];
	sData.strActualImageSizeT = '';
	sData.dblActualVoxelSizeX = [];
	sData.dblActualVoxelSizeY = [];
	
	%settings
	sData.dblSettingZoom = [];
	sData.dblSettingPhase = [];
	sData.boolSettingPMT1_Active = false;
	sData.dblSettingPMT1_Offset = [];
	sData.dblSettingPMT1_Gain = [];
	sData.strSettingPMT1_Unit = '';
	sData.boolSettingPMT2_Active = false;
	sData.dblSettingPMT2_Offset = [];
	sData.dblSettingPMT2_Gain = [];
	sData.strSettingPMT2_Unit = '';
	sData.dblSettingWavelength = [];
	sData.strSettingObjective = '';
	
	%image variables
	sData.intImageChannels = 0;
	sData.intImageBits = 0;
	sData.intImageSizeX = 0;
	sData.intImageSizeY = 0;
	sData.intImageSizeT = 0;
	sData.intImageSizeN = 0;
	
	%general
	sData.dblFrameDur = [];
	sData.strStartTime = '';
	sData.strStopTime = '';
	
	%meta data
	cellFieldnames = fieldnames(sData);
	vecVariablesRetrieved = false(1,length(cellFieldnames));
	
	%% get data
	%open file
	ptrFile = fopen(strXMLfile);
	while ~feof(ptrFile) && ~all(vecVariablesRetrieved)
		strLine = fgets(ptrFile);
		if strfind(strLine,'<ScannerSettingRecord')
			if strfind(strLine,'Identifier="dblSizeX"')
				%<ScannerSettingRecord Identifier="dblSizeX" Unit="µm" Description="Size-Width" Data="0" Variant="364.7" VariantType="5" />
				strOut = getFlankedBy(strLine,'Variant="','"');
				strField = 'dblActualImageSizeX';
				
				sData.(strField) = str2double(strOut);
				vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
			elseif strfind(strLine,'Identifier="dblSizeY"')
				%dblActualImageSizeY		<ScannerSettingRecord Identifier="dblSizeY" Unit="µm" Description="Size-Height" Data="0" Variant="364.7" VariantType="5" />
				strOut = getFlankedBy(strLine,'Variant="','"');
				strField = 'dblActualImageSizeY';
				
				sData.(strField) = str2double(strOut);
				vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
			elseif strfind(strLine,'Identifier="dblSizeZ"')
				%dblActualImageSizeZ		<ScannerSettingRecord Identifier="dblSizeZ" Unit="s" Description="Size-Depth" Data="0" Variant="493.6" VariantType="5" />
				strOut = getFlankedBy(strLine,'Variant="','"');
				strField = 'dblActualImageSizeZ';
				
				sData.(strField) = str2double(strOut);
				vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
			elseif strfind(strLine,'Identifier="dblVoxelX"')
				%dblActualVoxelSizeX		<ScannerSettingRecord Identifier="dblVoxelX" Unit="nm" Description="Voxel-Width" Data="0" Variant="713.7" VariantType="5" />
				strOut = getFlankedBy(strLine,'Variant="','"');
				strField = 'dblActualVoxelSizeX';
				
				sData.(strField) = str2double(strOut);
				vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
			elseif strfind(strLine,'Identifier="dblVoxelY"')
				%dblActualVoxelSizeY		<ScannerSettingRecord Identifier="dblVoxelY" Unit="nm" Description="Voxel-Height" Data="0" Variant="713.7" VariantType="5" />
				strOut = getFlankedBy(strLine,'Variant="','"');
				strField = 'dblActualVoxelSizeY';
				
				sData.(strField) = str2double(strOut);
				vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
			elseif strfind(strLine,'Identifier="dblZoom"')
				%dblSettingZoom			<ScannerSettingRecord Identifier="dblZoom" Unit="" Description="Zoom" Data="0" Variant="1.7" VariantType="5" />
				strOut = getFlankedBy(strLine,'Variant="','"');
				strField = 'dblSettingZoom';
				
				sData.(strField) = str2double(strOut);
				vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
			elseif strfind(strLine,'Identifier="nBit"')
				%intImageBits			<ScannerSettingRecord Identifier="nBit" Unit="bits" Description="Resolution" Data="0" Variant="12" VariantType="3" />
				strOut = getFlankedBy(strLine,'Variant="','"');
				strField = 'intImageBits';
				
				sData.(strField) = str2double(strOut);
				vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
			elseif strfind(strLine,'Identifier="nChannels"')
				%intImageChannels		<ScannerSettingRecord Identifier="nChannels" Unit="" Description="Channels" Data="0" Variant="2" VariantType="3" />
				strOut = getFlankedBy(strLine,'Variant="','"');
				strField = 'intImageChannels';
				
				sData.(strField) = str2double(strOut);
				vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
			end
		elseif strfind(strLine,'<FilterSettingRecord')
			strPMT = getFlankedBy(strLine,'ObjectName="PMT NDD','"');
			intPMT = str2double(strPMT);
			if intPMT > -1
				if strfind(strLine,'Attribute="State"')
					%boolSettingPMTx_Active	<FilterSettingRecord ObjectName="PMT NDD1" ClassName="CDetectionUnit" Attribute="State" Description="PMT NDD1" Data="1000" Variant="Active" VariantType="8" />
					strField = ['boolSettingPMT' strPMT '_Active'];
					if strfind(strLine,'Variant="Active"')
						sData.(strField) = true;
					else
						sData.(strField) = false;
					end
					vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
				elseif strfind(strLine,'Attribute="VideoOffset"')
					%dblSettingPMTx_Offset	<FilterSettingRecord ObjectName="PMT NDD1" ClassName="CDetectionUnit" Attribute="VideoOffset" Description="PMT NDD1 (Offs.)" Data="1000" Variant="-3.1" VariantType="5" Unit="%" />
					strOut = getFlankedBy(strLine,'Variant="','"');
					strField = ['dblSettingPMT' strPMT '_Offset'];
					
					sData.(strField) = str2double(strOut);
					vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
				elseif strfind(strLine,'Attribute="HighVoltage"')
					%dblSettingPMTx_Gain		<FilterSettingRecord ObjectName="PMT NDD1" ClassName="CDetectionUnit" Attribute="HighVoltage" Description="PMT NDD1 (HV)" Data="1000" Variant="851.9" VariantType="5" />
					strOut = getFlankedBy(strLine,'Variant="','"');
					strField = ['dblSettingPMT' strPMT '_Gain'];
					
					sData.(strField) = str2double(strOut);
					vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
				elseif strfind(strLine,'Attribute="HighVoltageUnit"')
					%strSettingPMTx_Unit		<FilterSettingRecord ObjectName="PMT NDD1" ClassName="CDetectionUnit" Attribute="HighVoltageUnit" Description="PMT NDD1 (HV_Unit)" Data="1000" Variant="V" VariantType="8" />
					strOut = getFlankedBy(strLine,'Variant="','"');
					strField = ['strSettingPMT' strPMT '_Unit'];
					
					sData.(strField) = strOut;
					vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
				end
			elseif strfind(strLine,'Attribute="Wavelength"')
				%dblSettingWavelength	<FilterSettingRecord ObjectName="Laser (MP, MP)" ClassName="CLaser" Attribute="Wavelength" Description="Laser wavelength" Data="0" Variant="810" VariantType="3" />
				strOut = getFlankedBy(strLine,'Variant="','"');
				strField = 'dblSettingWavelength';
				
				sData.(strField) = str2double(strOut);
				vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
			elseif strfind(strLine,'Attribute="Phase"')
				%dblSettingPhase			<FilterSettingRecord ObjectName="Scan Head" ClassName="CScanCtrlUnit" Attribute="Phase" Description="Phase" Data="0" Variant="2.62302586404212" VariantType="5" />
				strOut = getFlankedBy(strLine,'Variant="','"');
				strField = 'dblSettingPhase';
				
				sData.(strField) = str2double(strOut);
				vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
			elseif strfind(strLine,'Attribute="Objective"')
				%strSettingObjective		<FilterSettingRecord ObjectName="DM6000 Turret" ClassName="CTurret" Attribute="Objective" Description="Objective" Data="0" Variant="HC PL FLUOTAR  25.0x0.95 WATER " VariantType="8" />
				strOut = getFlankedBy(strLine,'Variant="','"');
				strField = 'strSettingObjective';
				
				sData.(strField) = strOut;
				vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
			end
		elseif strfind(strLine,'<DimensionDescription')
			if strfind(strLine,'DimID="X"')
				%intImageSizeX			<DimensionDescription DimID="X" NumberOfElements="512" Origin="0.00" Length="364.71" Unit="µm" BitInc="0" BytesInc="2" LogicalUnit="pixels" Voxel="0.714">
				strOut = getFlankedBy(strLine,'NumberOfElements="','"');
				strField = 'intImageSizeX';
				
				sData.(strField) = str2double(strOut);
				vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
			elseif strfind(strLine,'DimID="Y"')
				%intImageSizeY			<DimensionDescription DimID="Y" NumberOfElements="512" Origin="0.00" Length="364.71" Unit="µm" BitInc="0" BytesInc="1024" LogicalUnit="pixels" Voxel="0.714">
				strOut = getFlankedBy(strLine,'NumberOfElements="','"');
				strField = 'intImageSizeY';
				
				sData.(strField) = str2double(strOut);
				vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
			elseif strfind(strLine,'DimID="T"')
				%intImageSizeT			<DimensionDescription DimID="T" NumberOfElements="3381" Origin="0 s" Length="8m19.116s" Unit="" BitInc="0" BytesInc="1048576" LogicalUnit="" Voxel="0.146s">
				%strActualImageSizeT
				%dblFrameDur
				strField = 'intImageSizeT';
				sData.(strField) = str2double(getFlankedBy(strLine,'NumberOfElements="','"'));
				vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
				
				
				strField = 'strActualImageSizeT';
				sData.(strField) = getFlankedBy(strLine,'Length="','"');
				vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
				
				
				strField = 'dblFrameDur';
				sData.(strField) = str2double(getFlankedBy(strLine,'Voxel="','s">'));
				vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
			end
		elseif strfind(strLine,'<StartTime>')
			%strStartTime			<StartTime>6/13/2012 5:53:12 PM.445</StartTime>
			strOut = getFlankedBy(strLine,'<StartTime>','</StartTime>');
			strField = 'strStartTime';
			
			sData.(strField) = strOut;
			vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
		elseif strfind(strLine,'<EndTime>')
			%strStopTime				<EndTime>6/13/2012 6:01:31 PM.561</EndTime>
			strOut = getFlankedBy(strLine,'<EndTime>','</EndTime>');
			strField = 'strStopTime';
			
			sData.(strField) = strOut;
			vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
		elseif strfind(strLine,'<FrameCount>')
			%intImageSizeN      <FrameCount>6762 (2 channels, 3381 frames)</FrameCount>
			strOut = getFlankedBy(strLine,'<FrameCount>',' (');
			strField = 'intImageSizeN';
			
			sData.(strField) = str2double(strOut);
			vecVariablesRetrieved(strcmp(cellFieldnames,strField)) = true;
		end
	end
	if all(vecVariablesRetrieved);
		intOutFlag = 1;
	elseif any(vecVariablesRetrieved)
		intOutFlag = 0;
	else
		intOutFlag = -1;
	end
end

%{










 
%}