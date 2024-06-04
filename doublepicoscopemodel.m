%% This is a modified form of the following code for the purpose of using two osiliscopes with all eight chanelles at the highest sampling rate
% key values
TB=3;          % min TB value of 3 
samplingDelayInMicroSeconds= 100;          % sampling inteval = (TB-2)/125000000
 
sampleDurationInMilliSeconds=0.2;       trigval=1000;  % trigval= the trigger threshold
% time base will need to be increased around 33 milliseconds  

%% calculating number of samples for delay and sampling and posttrigsamples
delay=round(samplingDelayInMicroSeconds/(8e-3)); 
posttrigsamps=round(sampleDurationInMilliSeconds/(8e-6))+1;
pretrigsamps=0;

%% PicoScope 5000 Series (A API) Instrument Driver Oscilloscope Block Data Capture With Two Oscilloscopes Example

% These steps are:
%    
% # Create a device object   
% # Connect to the instrument 
% # Configure properties 
% # Invoke functions 
% # Disconnect from the instrument 

%% Clear command window and close any figures

clc;
close all;

%% Load configuration information

PS5000aConfig;

%% Device connection

% Check if an Instrument session using the device object |ps5000aDeviceObj|
% is still open, and if so, disconnect if the User chooses 'Yes' when prompted.
if (exist('ps5000aDeviceObj1', 'var') && ps5000aDeviceObj1.isvalid && strcmp(ps5000aDeviceObj1.status, 'open'))
    
    openDevice = questionDialog(['Device object ps5000aDeviceObj1 has an open connection. ' ...
        'Do you wish to close the connection and continue?'], ...
        'Device Object Connection Open');
    
    if (openDevice == PicoConstants.TRUE)
        
        % Close connection to device.
        disconnect(ps5000aDeviceObj1);
        delete(ps5000aDeviceObj1);
        
    else

        % Exit script if User selects 'No'.
        return;
        
    end
    
end

if (exist('ps5000aDeviceObj2', 'var') && ps5000aDeviceObj2.isvalid && strcmp(ps5000aDeviceObj2.status, 'open'))
    
    openDevice = questionDialog(['Device object ps5000aDeviceObj2 has an open connection. ' ...
        'Do you wish to close the connection and continue?'], ...
        'Device Object Connection Open');
    
    if (openDevice == PicoConstants.TRUE)
        
        % Close connection to device.
        disconnect(ps5000aDeviceObj2);
        delete(ps5000aDeviceObj2);
        
    else

        % Exit script if User selects 'No'.
        return;
        
    end
    
end

% Create a device object - provide the serial number as a second argument if required. 
ps5000aDeviceObj1 = icdevice('picotech_ps5000a_generic.mdd', '');
ps5000aDeviceObj2 = icdevice('picotech_ps5000a_generic.mdd', '');

% Connect device object to hardware.
connect(ps5000aDeviceObj1);
connect(ps5000aDeviceObj2);

%% Set channels
% Default driver settings used - use the |ps5000aSetChannel() function| to
% turn channels on or off and set voltage ranges, coupling, as well as
% analogue offset.
    % [status1.setChC] = invoke(ps5000aDeviceObj1, 'ps5000aSetChannel', 2, 1, 1, 8, 0.0);
    % [status1.setChD] = invoke(ps5000aDeviceObj1, 'ps5000aSetChannel', 3, 1, 1, 8, 0.0);
    % 
    % [status2.setChC] = invoke(ps5000aDeviceObj2, 'ps5000aSetChannel', 2, 1, 1, 8, 0.0);
    % [status2.setChD] = invoke(ps5000aDeviceObj2, 'ps5000aSetChannel', 3, 1, 1, 8, 0.0);

%% Set device resolution

% Max. resolution with 2 channels enabled is 15 bits.
[status1.setResolution, scope1.resolution] = invoke(ps5000aDeviceObj1, 'ps5000aSetDeviceResolution', 15);
[status2.setResolution2, scope2.resolution] = invoke(ps5000aDeviceObj2, 'ps5000aSetDeviceResolution', 15);

%% GET TIMEBASE and max samples
scope1.timebaseIndex = TB;
 [status1.getTimebase, scope1.timeIntervalNanoSeconds, scope1.maxSamples] = invoke(ps5000aDeviceObj1, 'ps5000aGetTimebase', scope1.timebaseIndex, 0);
scope2.timebaseIndex = TB;
[status2.getTimebase, scope2.timeIntervalNanoSeconds, scope2.maxSamples] = invoke(ps5000aDeviceObj2, 'ps5000aGetTimebase', scope2.timebaseIndex, 0);
maxsamples = min([scope2.maxSamples scope1.maxSamples]);
% check if the samples is less than the maximum if not provide a new
% timebase or reduce the number of samples
if (posttrigsamps+pretrigsamps)>= maxsamples+1 %max number of samples =   33554304 and the second scope makes it 4194240 
    newTB=3;
    newposttrigsamps=posttrigsamps;
    while maxsamples<=(newposttrigsamps+pretrigsamps)
    newTB=newTB+1;
    sampint=(newTB-2)/125000;
    newposttrigsamps=sampleDurationInMilliSeconds/sampint;
    end
    msg1=('your sampling duration requires a higher sampling interval of'+string(sampint)+' seconds\n would you like to procedde with this sampling rate?');
    opts=["yes","no"];
choice=menu(msg1,opts);
    if choice==1
        TB=newTB;
        posttrigsamps=newposttrigsamps;
    else
        posttrigsamps=maxsamples-pretrigsamps;
    end
end   


% Use the |ps5000aGetTimebase2()| function to query the driver as to the
% suitability of using a particular timebase index and the maximum number
% of samples available in the segment selected, then set the |timebase|
% property if required.
%
% To use the fastest sampling interval possible, enable one analog
% channel and turn off all other channels.
%
% Use a while loop to query the function until the status indicates that a
% valid timebase index has been selected. In this example, the timebase
% index of 65 is valid.

% Initial call to ps5000aGetTimebase2() with parameters:
%
% timebase      : 65
% segment index : 0

scope1.timebaseIndex = TB;

status1.getTimebase = PicoStatus.PICO_INVALID_TIMEBASE;

while (status1.getTimebase == PicoStatus.PICO_INVALID_TIMEBASE)
    
    [status1.getTimebase, scope1.timeIntervalNanoSeconds, scope1.maxSamples] = invoke(ps5000aDeviceObj1, 'ps5000aGetTimebase', scope1.timebaseIndex, 0);
    
    if (status1.getTimebase == PicoStatus.PICO_OK)
       
        break;
        
    else
        
        scope1.timebaseIndex = scope1.timebaseIndex + 1;
        
    end    
    
end

set(ps5000aDeviceObj1, 'timebase', scope1.timebaseIndex);

% Repeat for second device

scope2.timebaseIndex = TB;

status2.getTimebase = PicoStatus.PICO_INVALID_TIMEBASE;

while (status2.getTimebase == PicoStatus.PICO_INVALID_TIMEBASE)
    
    [status2.getTimebase, scope2.timeIntervalNanoSeconds, scope2.maxSamples] = invoke(ps5000aDeviceObj2, 'ps5000aGetTimebase', scope2.timebaseIndex, 0);
    
    if (status2.getTimebase == PicoStatus.PICO_OK)
       
        break;
        
    else
        
        scope2.timebaseIndex = scope2.timebaseIndex + 1;
        
    end    
    
end

set(ps5000aDeviceObj2, 'timebase', scope2.timebaseIndex);

%% Set simple trigger
% Set a trigger on channel A on both oscilloscopes, with an auto timeout -
% the default value for delay is used.

% Trigger properties and functions are located in the Instrument
% Driver's Trigger group.

triggerGroupObj1 = get(ps5000aDeviceObj1, 'Trigger');
triggerGroupObj1 = triggerGroupObj1(1);

triggerGroupObj2 = get(ps5000aDeviceObj2, 'Trigger');
triggerGroupObj2 = triggerGroupObj2(1);

% Set the |autoTriggerMs| property in order to automatically trigger the
% oscilloscope after 2 seconds if a trigger event has not occurred. Set to 0
% to wait indefinitely for a trigger event.

set(triggerGroupObj1, 'autoTriggerMs', 0);
set(triggerGroupObj2, 'autoTriggerMs', 0);
% set trigger delay where the value = the number of samples selected at the
% top
set(triggerGroupObj1, 'delay',delay)
set(triggerGroupObj2, 'delay',delay)
%set simple trigger goes channel,threshold,direction
% Channel     : channel 0 is A 1 is B.... and channel 4 is external trigger
% Threshold   : selected at the top of code
% Direction   : 2 (ps5000aEnuminfo.enPS5000AThresholdDirection.PS5000A_RISING)

[status1.setSimpleTrigger] = invoke(triggerGroupObj1, 'setSimpleTrigger', 4, trigval, 2);
[status2.setSimpleTrigger] = invoke(triggerGroupObj2, 'setSimpleTrigger', 4, trigval, 2);

%% Set block parameters and capture data
% Capture a block of data and retrieve data values for all enabled analog
% channels.

% Block data acquisition properties and functions are located in the 
% Instrument Driver's Block group.

blockGroupObj1 = get(ps5000aDeviceObj1, 'Block');
blockGroupObj1 = blockGroupObj1(1);

blockGroupObj2 = get(ps5000aDeviceObj2, 'Block');
blockGroupObj2 = blockGroupObj2(1);

% Set pre-trigger and post-trigger samples at the top of code
set(ps5000aDeviceObj1, 'numPreTriggerSamples', pretrigsamps);
set(ps5000aDeviceObj1, 'numPostTriggerSamples', posttrigsamps);

set(ps5000aDeviceObj2, 'numPreTriggerSamples', pretrigsamps);
set(ps5000aDeviceObj2, 'numPostTriggerSamples', posttrigsamps);
%% selecting values for storage
f=256;  % total values
vals=1:8:f;  % 8 is the number of values test
j=1;
k=1;
for h=1:length(vals) 
    j=1;
while j==1
currentvals=string(vals(k))+"-"+string(vals(k)+7);
name="vals"+currentvals;
msg=("do you want to record "+name);
opts=["yes","no"];
choice=menu(msg,opts);
if choice==1
    j=2;
elseif choice==2
newk=input("what is the first value you need to record\n>>>");
k=find(vals>=newk-7,1);
end
end
% Start the devices collecting data

% Capture a block of data:
%
% segment index: 0

status1.ps5000aRunBlock = invoke(blockGroupObj1, 'ps5000aRunBlock', 0);

status2.ps5000aRunBlock = invoke(blockGroupObj2, 'ps5000aRunBlock', 0);

% Poll the device driver to see if data is available

scope1.ready = PicoConstants.FALSE;
scope2.ready = PicoConstants.FALSE;

while (scope1.ready == PicoConstants.FALSE || scope2.ready == PicoConstants.FALSE)
    
    [status1.ready, scope1.ready] = invoke(blockGroupObj1, 'ps5000aIsReady');
    [status2.ready, scope2.ready] = invoke(blockGroupObj2, 'ps5000aIsReady');
    
    pause(0.01);
    
end


% Retrieve data values:

startIndex              = 0;
segmentIndex            = 0;
downsamplingRatio       = 1;
downsamplingRatioMode   = ps5000aEnuminfo.enPS5000ARatioMode.PS5000A_RATIO_MODE_NONE;

[scope1.numSamples, scope1.overflow, scope1.chA, scope1.chB, scope1.chC, scope1.chD] = invoke(blockGroupObj1, 'getBlockData', startIndex, segmentIndex, ... 
                                                                                            downsamplingRatio, downsamplingRatioMode);

[scope2.numSamples, scope2.overflow, scope2.chA, scope2.chB, scope2.chC, scope2.chD] = invoke(blockGroupObj2, 'getBlockData', startIndex, segmentIndex, ...
                                                                                            downsamplingRatio, downsamplingRatioMode);


%% new proccess data 
% format time to nano seconds
scope1.timeNs = double(scope1.timeIntervalNanoSeconds) * downsamplingRatio * double(0:scope1.numSamples - 1);
scope1.timeMs = scope1.timeNs / 1e6;
scope2.timeNs = double(scope2.timeIntervalNanoSeconds) * downsamplingRatio * double(0:scope2.numSamples - 1);
scope2.timeMs = scope2.timeNs / 1e6;

%% generate all plots from data
%picoscope 1 plots will have a black line and the picoscope 2 plots will
%have a red line
% Channel A 
close all
set(figure,'name',name)
subplot(4,2,1)
plot(scope1.timeMs, scope1.chA, 'b');
title('Scope B Channel 1');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
grid('on');

subplot(4,2,2)
plot(scope1.timeMs, scope1.chB, 'b');
title('Scope B Channel 1');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
grid('on');

subplot(4,2,3)
plot(scope1.timeMs, scope1.chC, 'b');
title('Channel C1');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
grid('on');

subplot(4,2,4)
plot(scope1.timeMs, scope1.chD, 'b');
title('Channel D1');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
grid('on');

subplot(4,2,5)
plot(scope2.timeMs, scope2.chA, 'r');
title('Channel A2');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
grid('on');

subplot(4,2,6)
plot(scope2.timeMs, scope2.chB, 'r');
title('Channel B2');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
grid('on');

subplot(4,2,7)
plot(scope2.timeMs, scope2.chC, 'r');
title('Channel C2');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
grid('on');

subplot(4,2,8)
plot(scope2.timeMs, scope2.chD, 'r');
title('Channel D2');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
grid('on');


%% stop device

[status1.stop] = invoke(ps5000aDeviceObj1, 'ps5000aStop');
[status2.stop] = invoke(ps5000aDeviceObj2, 'ps5000aStop');


% save data
datafile1=[scope1.timeMs',scope1.chA,scope1.chB,scope1.chC,scope1.chD,scope2.timeMs',scope2.chA,scope2.chB,scope2.chC,scope2.chD];
 
namef=append(name,'.csv');
writematrix(datafile1,namef);
k=k+1;
end
%% Disconnect device
% Disconnect device object from hardware

disconnect(ps5000aDeviceObj1);
delete(ps5000aDeviceObj1);

disconnect(ps5000aDeviceObj2);
delete(ps5000aDeviceObj2);