#pragma TextEncoding = "UTF-8"
//Parameters used in the experiment
//Copy paste these lines to the command window
//K19=1 //How much voltage drop is considered as a peak
//K18=10 //Datapoint delay between data and start of the curve
//K17=30 //Datapoint delay between data and end of the curve
//K11=1 //Start of the range of the Time dependence scan
//K10=70 //End of the range of the Time dependence scan


Function/WAVE MaximumLocations(input) 
//Returns the peak locations of the input
//Returns array of indexes of last location in curve where y>GlobalMaximum-Variability
	Wave input
	make/O/N=0 Peaks=0 //Initialize array of Peak index
	Variable GlobalMaximum=WaveMax(input)
	Variable inPeak=0 //Boolean value, if it is in a peak
	Variable Maxes=0 //Counter on how many peaks have been detected
	Variable i
	for(i=0; i<=numpnts(input); i+=1) //Loop through every datapoint in input
		if (!inPeak & (input[i]>=floor(GlobalMaximum-K19))) 
		//If it hasn't passed a peak yet, but passes upper threshold, then start peak
			inPeak=1
		elseif (inPeak & (input[i]<floor(GlobalMaximum-K19)))
		//If it has been through a peak, but passes lower threshold, then stop peak
			inPeak=0
			InsertPoints Maxes,1,Peaks //Expand the Peak array
			Peaks[Maxes]=i-1
			Maxes+=1
		endif
	endfor
	return Peaks
end

Function TausVsTemperature(voltage, times, temperature,lux)
	Wave voltage,times,temperature,lux
//Array of Datapoints of average temperature in curve
	make/O/N=0 TempLocations=0 
//Array of Time Constant values resulted from Curve Fitting
	make/O/N=0 Taus=0 
	variable n=numpnts(voltage)
	Wave Peak=MaximumLocations(voltage)
	make/O W_coef=0
	make/O W_sigma=0
	make/O W_fitConstants=0
	variable i
	variable EndCurve //Location of the end of the Curve
	for (i=0; i+1<numpnts(Peak); i+=1)
//Checks if the current Peak isn't the last peak
		if ((i+1)<numpnts(Peak))
			EndCurve=(Peak[i+1]-K17)
		else
			EndCurve=n
		endif
//Exponential curve-fitting between the peaks
		CurveFit exp_XOffset lux[Peak[i]+K18,EndCurve] /X=times /D
		insertPoints i,1,TempLocations,Taus
//Set the x Coordinate to average of temperature during curve
		TempLocations[i]=faverage(temperature, Peak[i]+K18,EndCurve) 
		Taus[i]=W_coef[2]
	endfor
	Display Taus vs TempLocations
END

Function IntensityTemp(voltage, times, temperature, lux)
//Intensity vs Temperature
//Modify K18 to change delay
	Wave voltage,times,temperature,lux
//Array of Datapoints of average temperature in curve
	make/O/N=0 TempLocations=0 
//Array of Time Constant values resulted from Curve Fitting
	make/O/N=0 Intensity=0 
	variable n=numpnts(voltage)
	Wave Peak=MaximumLocations(voltage)
	variable i
	variable EndCurve //Location of the end of the Curve
	for (i=0; i+1<numpnts(Peak); i+=1)
//Checks if the current Peak isn't the last peak
		if ((i+1)<numpnts(Peak))
			EndCurve=(Peak[i+1]-K17)
		else
			EndCurve=n
		endif
		insertPoints i,1,TempLocations,Intensity
//Set the x Coordinate to average of temperature during curve
		TempLocations[i]=faverage(temperature, Peak[i]+K18,EndCurve) 
		Intensity[i]=lux[Peak[i]+K18]
	endfor
	Display Intensity vs TempLocations
END

//Time distance between time datapoints
Function TimeInterval(times)
	Wave times
	make/O/N=0 Point
	make/O/N=0 Interval
	variable n=numpnts(times)
	variable i
	for (i=1; i+1<n; i+=1)
		insertPoints i,1,Point,Interval
		Interval[i]=times[i]-times[i-1]
		Point[i]=i
	endfor
	Display Interval vs Point
END

//Optimum Temperature for Intensity vs Time
Function IntensityVsTime(voltage,times,temperature)
	Wave voltage,times,temperature
	make/O/N=0 TimeSteps=0
	make/O/N=0 MaxTemp=0

	variable i
	for (i=K11; i<K10; i++)
		make/O/N=0 TempLocations=0 
		make/O/N=0 Intensity=0
		variable n=numpnts(voltage)
		Wave Peak=MaximumLocations(voltage)
		variable j
		variable EndCurve
		for (j=0; j+1<numpnts(Peak); j+=1)
			if ((j+1)<numpnts(Peak))
				EndCurve=(Peak[j+1]-K17)
			else
				EndCurve=n
			endif
			insertPoints j,1,TempLocations,Intensity
			TempLocations[j]=faverage(temperature, Peak[j]+i,EndCurve) 
			Intensity[j]=voltage[Peak[j]+i]
		endfor
		make/O W_coef=0
		insertPoints i-K11,1,TimeSteps,MaxTemp
		CurveFit poly 3, Intensity /X=TempLocations /D
		TimeSteps[i-K11]=i
		MaxTemp[i-K11]=(-1*W_coef[1])/(2*W_coef[2])
	endfor
	Display MaxTemp vs TimeSteps
END

//AuC (of Lux vs Time) vs Temperature
//Change K10, to extend right end of the x axis
Function AuCvsTemp(voltage,times,temperature,lux)
	Wave voltage,times,temperature,lux
	Integrate /METH=1 lux /X=times /D=AuC
	make/O/N=0 TempLocations=0 
	make/O/N=0 RelativeAuC=0 

	variable n=numpnts(voltage)
	Wave Peak=MaximumLocations(voltage)
	variable i
	variable EndCurve 
	for (i=0; i+1<numpnts(Peak); i+=1)
		if ((i+1)<numpnts(Peak))
			EndCurve=(Peak[i+1]-K17)
		else
			EndCurve=n
		endif
		insertPoints i,1,TempLocations,RelativeAuC
		TempLocations[i]=faverage(temperature, Peak[i]+K18,EndCurve) 
		RelativeAuC[i]=AuC[EndCurve]-AuC[Peak[i]+K18]
	endfor
	Display RelativeAuC vs TempLocations
END

//Optimum Temperature for AuC (of Lux vs Time) vs Time
//Change K10, to extend right end of the x axis
Function AuCVsTime(voltage,times,temperature,lux)
	Wave voltage,times,temperature,lux
	Integrate /METH=1 lux /X=times /D=AuC
	make/O/N=0 TimeSteps=0
	make/O/N=0 MaxTemp=0

	variable i
	for (i=K11; i<K10; i++)
		make/O/N=0 TempLocations=0 
		make/O/N=0 RelativeAuC=0 
		variable n=numpnts(voltage)
		Wave Peak=MaximumLocations(voltage)
		variable j
		variable EndCurve
		for (j=0; j+1<numpnts(Peak); j+=1)
			if ((j+1)<numpnts(Peak))
				EndCurve=(Peak[j+1]-K17)
			else
				EndCurve=n
			endif
			insertPoints j,1,TempLocations,RelativeAuC
			TempLocations[j]=faverage(temperature, Peak[j]+i,EndCurve) 
			RelativeAuC[j]=AuC[EndCurve]-AuC[Peak[j]+i]
		endfor
		make/O W_coef=0
		insertPoints i-K11,1,TimeSteps,MaxTemp
		CurveFit poly 3, RelativeAuC /X=TempLocations /D
		TimeSteps[i-K11]=i
		MaxTemp[i-K11]=(-1*W_coef[1])/(2*W_coef[2])
	endfor
	Display MaxTemp vs TimeSteps
END
