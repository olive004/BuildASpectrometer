% STUDENT: Olivia Gallupova

% Build a spectrometer: Day 3 Matlab "self-directed learning"
%% Pt. 1: Plotting data from white LED light 
% During measurements, we saw that the peak is around 280nm 
clear

WhiteSpectrum = csvread('WhiteLED.csv');
% plot(WhiteSpectrum);


%% Pt 2: Creating calibration curve (for white light data)
Purple = 405;           % excitation wavelengths  
%Green = 532;               we did not actually use a green laser
Red = 650;

% Extracting Laser spectrum data
RedSpectrum = csvread('RedLaser.csv'); 
PurpleSpectrum = csvread('UVLaser.csv'); 


% Storing indices (pixel values) of the laser spectra
% find gets the index of the maximum vale in the spectrum 
pixelValueRed = find(RedSpectrum == max(RedSpectrum));
pixelValuePurple = find(PurpleSpectrum == max(PurpleSpectrum));

% defining the simultaneous equations ( Y = m*X +c )
syms c m;
purpleEqn = m*pixelValuePurple + c == Purple;
redEqn = m*pixelValueRed + c == Red;

% Solving simultaneous equations for m and  
% From Matlab doc: sol = solve([eqn1, eqn2, eqn3], [x, y, z]);
unknowns = solve([purpleEqn, redEqn], [m, c]);
mFound = unknowns.m;
cFound = unknowns.c;

% Initially, the following gave a wavelength calibration that reversed the
% intensities spectrum, but reversing the order of extraction of the data
% in the Whitespectrum fixes this problem
% --> Wavelengths = (1:length(WhiteSpectrum)) * mFound + cFound;
% The following is the corrected version:
WLWhite = (length(WhiteSpectrum):-1:1) * mFound + cFound;

figure(2);
plot(WLWhite, WhiteSpectrum);
title('White LED Intensity Spectrum');
xlabel('Wavelength (nm)');
ylabel('Intensity (arbitrary units)');

%% Pt.3: Calculating Absorbance

% Extracting the dye spectra and their blanks from csv files
RedDyeSpectrum = csvread('RedDye.csv');
RedBlank = csvread('BlankForRed.csv');
BlueDyeSpectrum = csvread('BlueDye.csv');
BlueBlank = csvread('BlankForBlue.csv');

% Calculating the absorbances for both dyes
AbsorbancesRed = 1 - RedDyeSpectrum ./ RedBlank;
AbsorbancesBlue = 1 - BlueDyeSpectrum ./ BlueBlank;

% Getting the right calibration of wavelengths:
% Although we have a vector of wavelengths already from the white LED 
% spectrum, it is better to recreate these for the spectrum data of the 
% dyes so that it is certain the matrix dimensions will be the same 
% (in case the measurement taken was slightly different for whatever reason).
% In this case the resolution is the same for all our data, so one could
% leave out the reextraction of calibrated wavelengths for the dyes.
% WLRed -- Wavelength spectrum for red dye
WLRed = (1:length(AbsorbancesRed)) * mFound + cFound; 
WLBlue = (1:length(AbsorbancesBlue)) * mFound + cFound;

% Plotting Absorbance vs Wavelength: 
figure(3);
subplot(2, 2, 1);
plot(WLRed, AbsorbancesRed);
title('Red Dye Absorbances');
xlabel('Wavelength (nm)');
ylabel('Absorbance');
redx=xlim;
redy=ylim;

subplot(2, 2, 2);
plot(WLBlue, AbsorbancesBlue);
title('Blue Dye Absorbances');
xlabel('Wavelength (nm)');
ylabel('Absorbance');
bluex=xlim;
bluey=ylim;

% Removing noisy data btwn 420 - 680nm from spectrum:

highidxRed = min(find(WLRed < 420));      % find(A <420) returns an array of 
% all the indices of elements in A less than 420, we only want the one
% right next to 420 though so the minimum of those indexes is extracted.
% The same principle applies to the index of the wavelength just above
% 680nm, which will be the biggest index of all the indices for
% wavelengths greater than 680nm, hence the use of max() instead of min().
lowidxRed = max(find(WLRed >680));

highidxBlue = min(find(WLBlue < 420));
lowidxBlue = max(find(WLBlue >680));

indexesRangeWh = find( WLWhite > 420 & WLWhite <680);

% The absorbance and wavelength values in the 420-680nm range
filteredAbsRed = AbsorbancesRed( lowidxRed : highidxRed );
filteredWLRed = WLRed( lowidxRed : highidxRed );

filteredAbsBlue = AbsorbancesBlue( lowidxBlue : highidxBlue );
filteredWLBlue = WLBlue( lowidxBlue : highidxBlue );

% Getting the white LED intensity and wavelength vectors btwn 420 and 680nm
filteredIntsWh = WhiteSpectrum(indexesRangeWh);
filteredWLWh = WLWhite( WLWhite > 420 & WLWhite <680);


% Plotting the filtered dye absorbance and wavelength values
subplot(2,2,3)
plot(filteredWLRed, filteredAbsRed, 'r')
title('Filtered Red Dye Abs');
xlabel('Wavelength (nm)');
ylabel('Absorbance');
xlim(redx);     % Making the axes the same as the unfiltered absorbances
ylim(redy);     % so that it's easier to see the filtering
subplot(2,2,4)
plot(filteredWLBlue, filteredAbsBlue, 'r')
title('Filtered Blue Dye Abs');
xlabel('Wavelength (nm)');
ylabel('Absorbance');
xlim(bluex);
ylim(bluey);

%% Pt. 4: Smoothing the data

% We have a moving average with a window going from -k to +k. 

k = 100;                             % k*2+1 point moving ave: this is the moving window
nR = length(filteredAbsRed);         % The amount of points we are sampling
nB = length(filteredAbsBlue);        % from the red (R) and blue (B) dye spectra, and 
nW = length(filteredIntsWh);         % for the white (W) LED intensity spectrum.

% To smooth the data, I am looping through each point in the data and
% summing the value k to the left and k to the right of the index i, then
% dividing by the size of the summation, ie 2*k+1. 
% To make sure the moving average doesn't index out of bounds, the lower
% and upper indices of the current window are determined by max and min
% functions that check to make sure the indices to the left are greater than
% 1 and those to the right of window's center are less than the size of the vector holding
% the smoothened data.

% initializing the smoothed absorbances with a zero vector the size of the
% original vector with the spectrum data
smoothedAbsRed = zeros(size(filteredAbsRed));

for i = 1:nR      % truncated index
% The next two variables truncate the window on each iteration
    lowerHalfMinimumIndex = max(i-k, 1);
    upperHalfMaximumIndex = min(i+k, nR);
    
    % Finding the next set of values to sum
    nextSum = filteredAbsRed(lowerHalfMinimumIndex:upperHalfMaximumIndex);       
    currentMean = sum(nextSum)/length(nextSum);     % calculating the mean
    
    smoothedAbsRed(i) = currentMean; % Making the value at the current index equal to the mean
end


% REPEATING THE SMOOTHING for the blue dye
% See the code for the red dye for comments
smoothedAbsBlue = zeros(size(filteredAbsBlue));

for i = 1:nB      
    lowerHalfMinimumIndex = max(i-k, 1);
    upperHalfMaximumIndex = min(i+k, nR);
    nextSum = filteredAbsBlue(lowerHalfMinimumIndex:upperHalfMaximumIndex);
    currentMean = sum(nextSum)/length(nextSum);
    smoothedAbsBlue(i) = currentMean;
end

% REPEATING SMOOTHING FOR THE WHITE INTENSITY DATA
smoothedAbsWh = zeros(size(filteredIntsWh));
for i = 1:nW      % truncated index
    lowerHalfMinimumIndex = max(i-k, 1);
    upperHalfMaximumIndex = min(i+k, nW);
    nextValue = filteredIntsWh(lowerHalfMinimumIndex:upperHalfMaximumIndex);
    currentMean = sum(nextValue)/length(nextValue);
    smoothedAbsWh(i) = currentMean;
end


% Plot of the smoothed and truncated data
figure(4)
subplot(2, 2, 1);
plot(filteredWLRed, smoothedAbsRed)
title('Smoothed Red Dye Absorbances');
xlabel('Wavelength (nm)');
ylabel('Absorbance');

subplot(2, 2, 2);
plot(filteredWLBlue, smoothedAbsBlue)
title('Smoothed Blue Dye Absorbances');
xlabel('Wavelength (nm)');
ylabel('Absorbance');


% To compare, here is a graph using movmean() instead of the manual
% solution; the is just to check for accuracy of the data

smoothedAbsRedMOVMEAN = movmean(filteredAbsRed, [k k]);
smoothedAbsBlueMOVMEAN = movmean(filteredAbsBlue, [k k]);

subplot(2,2,3);
plot(filteredWLRed, smoothedAbsRedMOVMEAN, 'g');
title('Smoothed Red Dye Abs: Movmean()');
xlabel('Wavelength (nm)');
ylabel('Absorbance');

subplot(2, 2, 4);
plot(filteredWLBlue, smoothedAbsBlueMOVMEAN, 'g');
title('Smoothed Blue Dye Abs: Movmean()');
xlabel('Wavelength (nm)');
ylabel('Absorbance');


% Plotting smoothed white LED intensities
figure(6)
plot(filteredWLWh, smoothedAbsWh, 'b')
title('Smoothed White Intensities');
xlabel('Wavelength (nm)');
ylabel('Intensity (arbitrary units)');
