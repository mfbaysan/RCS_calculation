%%DEFINE RADAR SCENARIO
% Radar parameters      -> define radar parameters in our case
lambda      = 10e-2;        % Wavelength (m)
numPulses   = 4000;         % Number of pulses
prf         = 1e3;          % PRF (Hz)
antBw       = 8;            % Antenna beamwidth (deg)
pulseBw     = 10e6;         % Bandwidth (Hz)
fs          = pulseBw;      % Sampling frequency (Hz)
tpd         = 1e-6;         % Pulse width (sec)
maxRange    = 20e3;         % Maximum range (m)

% Wind turbine parameters  -> change to object later on
distToRadar = 10e3;         % Range (m)
towerHeight = 48;           % Turbine tower height (m)

% Calculate range resolution (m) 
rangeRes = bw2rangeres(pulseBw)

% Create a radar scenario
pri = 1/prf; 
dur = pri*(numPulses - 1); 
scene = radarScenario('UpdateRate',prf,'StopTime',dur);

%%CONFIGURE RADAR SENSOR (??) 
% Create a radar transceiver
[freq,c] = wavelen2freq(lambda); 
ant = phased.SincAntennaElement('Beamwidth',antBw);
rdr = radarTransceiver('MountingAngles',[90 0 0],'RangeLimits',[0 maxRange]); 
rdr.TransmitAntenna.Sensor = ant;
rdr.TransmitAntenna.OperatingFrequency = freq;
rdr.ReceiveAntenna.Sensor = ant;
rdr.ReceiveAntenna.OperatingFrequency = freq;
rdr.Receiver.SampleRate = fs; 

% Estimate the antenna gain
antennaGain = beamwidth2gain(antBw,'CosineRectangular'); 
rdr.Transmitter.Gain = antennaGain;
rdr.Receiver.Gain = antennaGain;

% Configure the LFM signal of the radar
rdr.Waveform = phased.LinearFMWaveform('SampleRate',fs,'PulseWidth',tpd, ...
    'PRF',prf,'SweepBandwidth',pulseBw); 

% Configure the LFM signal of the radar
rdr.Waveform = phased.LinearFMWaveform('SampleRate',fs,'PulseWidth',tpd, ...
    'PRF',prf,'SweepBandwidth',pulseBw); 

% Mount radar sensor on a stationary platform
time = 0:pri:(numPulses - 1)*pri; 
rdrplatPos = [0 -distToRadar towerHeight];
rdrplat = platform(scene,'Position',rdrplatPos);
rdrplat.Sensors = rdr;

%% OBJECT DEFINITION

% Wind turbine blade parameters
numBlades   = 3;            % Number of blades
bladeLength = 30;           % Blade length (m)
bladeWidth  = 2;            % Blade width (m)
rotRateDeg  = 36;           % Rotation rate (deg/sec)
rotAng0     = [55 175 295]; % Initial blade rotation angle (deg)

% Wind turbine tower parameters
towerPos = [-bladeWidth;0;towerHeight]; 

% Approximate far field distance in km assuming the entire wind turbine is
% modeled
RkmFarField = 2*towerHeight^2/lambda*1e-3 % Far-field range (km)

Rhorizon = horizonrange(towerHeight)*1e-3 % Radar horizon (km)

% Approximate far field distance in km using smaller components
res         = 2;                  % Rectangular plate length (m)
RkmFarField = 2*res^2/lambda*1e-3 % Far-field range (km)

% Calculate the RCS of a small rectangular plate 
a = bladeWidth;
b = res; 
rcsSigBlade = helperRcsBladePart(a,b,c,freq);  

% Determine the number of segments for the wind turbine blade
bladeSegmentHeights = 0:res:bladeLength;
numBladeSegments = numel(bladeSegmentHeights) - 1

% Create wind turbine parts and their trajectories
bladePos = helperCreateWindTurbine(scene,rcsSigBlade,towerHeight,bladeLength,res,rotAng0,rotRateDeg);

% Plot wind turbine movement
helperPlotWindTurbine(rdrplatPos,antBw,distToRadar,bladePos)

% Initialize output datacube
ib = 0; 

timeVec = 0:1/fs:(pri - 1/fs); 
rngVec = [0 time2range(timeVec(2:end))];
idx = find(rngVec > maxRange,1,'first'); 
numSamplesPulse = tpd*fs + 1;
numSamples = idx + numSamplesPulse - 1; 
iqsig = zeros(numSamples,numPulses);

% Create range time scope
rti = phased.RTIScope('IQDataInput',true,'RangeLabel','Range (km)', ...
    'RangeResolution',time2range(1/fs), ...
    'TimeResolution',pri,'TimeSpan',1,'SampleRate',fs);

% Collect data 
disp('Simulating IQ...')

while advance(scene)
    % Collect IQ
    ib = ib + 1; 
    tmp = receive(scene); 
    iqsig(:,ib) = tmp{1};

    % Update scope
    matchingCoeff = getMatchedFilter(rdr.Waveform); 
    rti(iqsig(:,ib),matchingCoeff); 

    % Update progress 
    if mod(ib,500) == 0
        fprintf('\t%.1f %% completed.\n',ib/numPulses*100); 
    end
end

function rcsSigBlade = helperRcsBladePart(a,b,c,fc)
% Helper to return rcsSignature object of a small rectangular PEC plate

% Calculate the RCS of a small rectangular plate 
az = -180:1:180; 
el = -90:1:90; 
rcsBlade = helperRcsPlate(a,b,c,fc,az,el); 

% Rotate the rectangular plate such that it is vertical in the y-z plane
el = el + 90; 
el = mod(el + 90,180) - 90; 
[el,idxSort] = unique(el); 
rcsBlade = rcsBlade(idxSort,:); 
rcsSigBlade = rcsSignature('Pattern',pow2db(rcsBlade),'Azimuth',az,'Elevation',el); 
end

function rcsPat = helperRcsPlate(a,b,c,fc,az,el)
%rcsplate   Radar cross section for flat rectangular plate
%   RCSPAT = helperRcsPlate(A,B,C,FC,AZ,EL) returns the radar cross section
%   (RCS) (in square meters) for a flat rectangular plate at given
%   frequencies. C is a scalar representing the propagation speed (in m/s),
%   and FC is an L-element vector containing operating frequencies (in Hz)
%   at which the RCS patterns are computed.
%
%   AZ are the azimuth angles (in degrees) specified as a Q-element array,
%   and EL are the elevation angles (in degrees) specified as a P-element
%   array. The pairs between azimuth and elevation angles define the
%   directions at which the RCS patterns are computed. This helper assumes
%   a monostatic configuration. HH = VV polarization for this case. 
%
%   The flat rectangular plate is assumed to be lying in the x-y plane
%   centered at the origin. Input A is the length of the plate along the
%   x-axis (long side) and B is the width of the plate along the y-axis
%   (short side).
%
%   The RCS calculation assumes a perfectly conducting plate and is an
%   approximation that is only applicable in the optical region, where the
%   target extent is much larger than the wavelength.
%
%   RCSPAT is a P-by-Q-by-L array containing the RCS pattern (in square
%   meters) of the cylinder at given azimuth angles in AZ, elevation angles
%   in EL, and frequencies in FC.

% Reference
% [1] C. A. Balanis, "Advaced Enginnering Electromagnetics", John
%     Wiley&Sons, USA, 1989.

[elg,azg,fcg] = ndgrid(el(:),az(:),fc(:));
lambda = c./fcg;
theta = pi/2 - deg2rad(elg);
phi = deg2rad(azg);

k = 2*pi./lambda; % Wave number

% Balanis
X = k*a/2.*sin(theta).*cos(phi);
idxPos = phi>=0;
Y = zeros(size(X));
Y(idxPos) = k(idxPos)*b/2.*(sin(theta(idxPos)).*sin(phi(idxPos)) + sin(theta(idxPos)));
Y(~idxPos) = k(~idxPos)*b/2.*(sin(theta(~idxPos)).*sin(phi(~idxPos)) - sin(theta(~idxPos)));
rcsPat = 4*pi*(a*b./lambda).^2.*(cos(theta).^2.*sin(phi).^2 + cos(phi).^2) ...
    .*sincFcn(X).^2.*sincFcn(Y).^2 + eps; % TE
end

function y = sincFcn(x)
%sincFcn   Computes the sinc function
%   Y = sincFcn(X) computes the sinc function 
%      Y = sin(X)/X.
%   The output Y is the same length as the input X. X is assumed to be in
%   radians.

y = sin(x)./x;
y(isnan(y)) = 1; 
end

function bladePos = helperCreateWindTurbine(scene,rcsSigBlade,towerHeight,bladeLength,res,rotAng0,rotRateDeg)
% Helper to create a wind turbine for a radar scenario

% Determine the segment heights and number of segments
bladeSegmentHeights = 0:res:bladeLength;
numBladeSegments = numel(bladeSegmentHeights) - 1;
numBlades = numel(rotAng0); 

% Determine the blade segment center points 
bladeSegmentCtrs = zeros(3,numBladeSegments);  
bladeSegmentCtrs(3,:) = (bladeSegmentHeights(1:end-1) + bladeSegmentHeights(2:end))./2; 
rotRateRad = deg2rad(rotRateDeg); 

% Simulation time
prf = scene.UpdateRate; 
time = scene.SimulationTime:1/prf:scene.StopTime;
numPulses = numel(time); 

% Initialize values
platCount = 0; 
bladePos = cell(1,numBlades*numBladeSegments); 

% Create the moving blades
for ib = 1:numBlades % Loop over blades
    for is = 1:numBladeSegments % Blade segments
        % Blade trajectory
        bladeOrient = zeros(3,3,numPulses);
        bladeWaypoints = zeros(numPulses,3);
        bladeVel = zeros(numPulses,3);
        platCount = platCount + 1;
        for tt = 1:numPulses % Loop over time 
            rotAng = rotAng0(ib) + rotRateDeg*time(tt);
            bladeOrient(:,:,tt) = rotx(rotAng);
            bladeWaypoints(tt,:) = bladeOrient(:,:,tt)*bladeSegmentCtrs(:,is);
            bladeVel(tt,:) = cross([rotRateRad;0;0],bladeOrient(:,:,tt)*bladeSegmentCtrs(:,is)); % v = w x r
        end
        % Shift blade waypoints by tower height
        bladeWaypoints(:,3) = bladeWaypoints(:,3) + towerHeight;

        % Create trajectory from waypoints
        bladePos{platCount} = bladeWaypoints;
        traj = waypointTrajectory('SampleRate',prf,'TimeOfArrival',time, ...
            'Waypoints',bladeWaypoints,'Velocities',bladeVel, ...
            'Orientation',bladeOrient);

        % Assign trajectory and rectangular plate RCS to platform
        platform(scene,'Trajectory',traj,'Signatures',rcsSigBlade); 
    end
end
end

function helperPlotWindTurbine(rdrplatPos,antBw,distToRadar,bladePos)
% Helper to plot position of wind turbine

% Plot radar
hFig = figure;
currentPos = hFig.Position;
hFig.Position = [currentPos(1) currentPos(2) 1000 currentPos(4)];
tiledlayout(1,2)
hAx1 = nexttile;
co = colororder;
plot3(rdrplatPos(1),rdrplatPos(2),rdrplatPos(3),'o', ...
    'MarkerFaceColor',co(1,:),'MarkerEdgeColor',co(1,:))
hold on
grid on
view([40 40])
axis('equal')
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Scene Geometry')

% Plot radar beam
distBeam = distToRadar*1.1;
coneRad = distBeam.*sind(antBw/2);
[x,y,z] = cylinder([0 coneRad]);
z = z*distBeam;
M = makehgtform('translate',rdrplatPos,'zrotate',pi/2,'yrotate',pi/2);
surf(x,y,z,'Parent',hgtransform('Matrix',M),'LineStyle','none','FaceAlpha',0.2);
hAx2 = nexttile; 

% Plot wind turbine
numPulses = size(bladePos{1},1);
for ii = 1:100:numPulses
    x = cellfun(@(x) x(ii,1),bladePos);
    y = cellfun(@(x) x(ii,2),bladePos);
    z = cellfun(@(x) x(ii,3),bladePos);
    if ii == 1
        hTurbine1 = plot3(hAx1,x,y,z,'o', ...
            'MarkerFaceColor',co(2,:),'MarkerEdgeColor',co(2,:));
        hTurbine2 = plot3(hAx2,x,y,z,'o', ...
            'MarkerFaceColor',co(2,:),'MarkerEdgeColor',co(2,:));
        view(hAx2,[40 40])
        axis(hAx2,'equal')
        xlim(hAx2,[-100 100])
        ylim(hAx2,[-100 100])
        zlim(hAx2,[0 100])
        grid(hAx2,'on')
        xlabel(hAx2,'X')
        ylabel(hAx2,'Y')
        zlabel(hAx2,'Z')
        title(hAx2,'Wind Turbine Movement')
    else
        hTurbine1.XData = x;
        hTurbine1.YData = y;
        hTurbine1.ZData = z;
        hTurbine2.XData = x;
        hTurbine2.YData = y;
        hTurbine2.ZData = z;
    end
    pause(0.1)
end

end