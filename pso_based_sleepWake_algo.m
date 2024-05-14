numSensors = 200;                   % Number of sensors
areaSize = 100;                     % Size of the area
numIterations = 5;                  % Number of iterations
numParticles = 20;                  % Number of particles
weights = [0.5, 0.5, 0.1, 0.2];     % Weights for fitness function
energyConsumptionRate = 240;        % Energy consumption rate
sensing_range = 6;                  % Sensing range of sensors
energySumOverTime = [];             % Array to store energy over time
coverageOverTime = [];              % Array to store coverage over time
numActiveSensorsTime = [];
numComponentsTime= [];
timePoints = [];                    % Array to store time points
transmissionRange = 7;              % Transmission range of sensors
thresholdEnergy = 500;
standbyEnergyConsumptionRate = 2;
% 24 Vdc
% mAh X voltage x 3.6 = Joules of energy
% 30 mAh 
% 2592 Joules 
% power consumption =  2 W  ==> ( 240 J/ 2min)

% Randomly generate sensor locations
%sensorLocations = rand(numSensors, 2) * areaSize;

% Randomly generate sensor energy
initialEnergyLevels = ones(numSensors, 1) * 2592;



sensorLocations = [1.6658 5.4770; 1.1176 11.8469; 1.7714 15.7595; 2.2875 26.7268; 3.3909 31.4062; 3.9778 39.2143; 0.9256 44.6298; 2.6698 54.7955; 3.9772 62.2063; 5.2175 65.4518; 1.8337 75.1116; 5.9640 79.7707; 5.7480 90.3292; 4.6323 95.8552; 8.4507 1.0801; 8.5410 9.5469; 8.9490 18.1031; 11.5520 27.6087; 9.0901 29.9643; 13.2426 42.0906; 13.5210 47.1681; 12.5417 51.2799; 8.1360 60.2408; 12.9196 67.5292; 9.1161 74.8550; 10.3059 82.5946; 9.3490 89.8907; 12.1318 97.9224; 17.4576 5.1982; 17.2373 11.9059; 15.7042 20.6407; 18.4026 26.9091; 19.3402 31.0193; 16.6739 40.2404; 16.9100 48.8924; 16.2533 54.9758; 15.1725 62.1504; 16.9602 70.0000; 16.5524 77.7790; 16.7029 79.8314; 18.0645 88.6034; 19.8305 97.5184; 26.2592 5.1518; 23.7468 11.7440; 24.5016 18.5470; 24.7426 24.7188; 22.4300 32.6852; 27.5787 36.7121; 25.9532 44.3479; 25.8404 56.2369; 27.8029 61.4418; 24.3509 70.1452; 26.1825 73.4195; 26.8053 80.2688; 23.3847 87.8046; 24.0117 94.6557; 34.9589 5.6725; 30.3666 11.5671; 31.5489 16.5782; 29.5571 27.1892; 32.8956 31.7696; 33.0064 37.8636; 33.6855 48.5059; 30.3955 55.9801; 33.6502 61.1935; 31.8078 67.5373; 31.1078 72.3087; 30.8624 81.3435; 34.3931 86.7970; 34.0881 94.5201; 40.2303 6.3869; 38.4029 10.7754; 36.8347 16.0509; 39.2500 27.7634; 39.2332 31.6099; 38.6616 37.8419; 39.7125 45.1592; 40.9704 53.2792; 38.4496 62.6469; 40.5280 69.2378; 37.8468 75.5712; 40.5202 84.6295; 40.9544 90.0226; 39.0981 99.1193; 47.2799 3.2174; 45.3129 8.5437; 46.5530 18.8766; 44.5049 24.1150; 48.7385 31.6900; 47.8778 41.3787; 48.7971 45.3208; 44.5093 51.2705; 44.7000 59.6412; 43.6242 70.5173; 47.0121 75.8969; 48.4237 81.7324; 49.0288 90.5281; 46.2329 97.3211; 54.2954 5.6982; 53.5932 10.2086; 52.1319 20.4199; 53.5550 23.4059; 55.5122 30.3365; 51.5064 37.5056; 53.6232 44.2729; 54.0419 53.5307; 53.3056 60.0140; 55.4033 67.3452; 50.8712 77.6481; 50.8473 83.9251; 54.9871 90.1175; 54.7996 98.3597; 59.6455 2.9722; 60.0482 10.2360; 61.0871 18.5349; 59.5110 26.8106; 61.7171 30.5066; 60.2300 39.0418; 61.1153 46.4862; 60.6320 51.8960; 61.5934 61.9154; 60.7024 68.7448; 60.5375 74.1033; 61.8321 83.6544; 58.9937 91.0541; 59.4479 98.7620; 67.1493 4.4686; 66.4311 13.0621; 67.9823 17.8682; 67.7159 27.5260; 69.6295 29.9742; 66.1196 38.0827; 66.5684 47.3796; 68.0040 51.6396; 65.5265 63.0138; 66.7099 67.8694; 70.2084 76.3180; 65.8635 83.8547; 69.6659 90.7940; 70.3035 97.7957; 76.4221 3.3879; 77.7414 11.9223; 75.1761 17.9224; 76.9867 23.7757; 77.2501 34.0654; 77.8025 38.6208; 76.2937 46.1588; 77.1004 55.8019; 75.4758 58.9190; 73.6909 70.4088; 76.5905 73.3044; 73.4387 81.7291; 74.8182 87.1209; 72.8099 97.8621; 83.4120 3.2774; 84.6712 9.6669; 82.6190 17.8824; 82.4878 26.6596; 81.8789 30.1430; 82.8889 38.0193; 84.6949 44.2821; 83.9276 55.5853; 80.1616 61.2984; 83.5956 66.9589; 81.7378 76.9871; 81.4994 83.7391; 82.9886 86.5475; 83.2990 94.4323; 88.9757 4.3947; 91.4604 8.1077; 87.7578 18.7628; 90.1182 25.6281; 89.9832 29.5658; 86.7720 41.0677; 89.2144 46.5450; 87.4193 56.2124; 88.8687 58.1101; 87.8521 69.9216; 88.9609 72.2781; 91.9848 81.0787; 87.3901 88.0255; 89.1241 94.2225; 96.9373 0.8132; 99.0221 11.8056; 94.4176 17.4256; 95.5849 25.8913; 94.3973 30.9401; 94.4635 37.4246; 94.3959 47.8975; 94.9507 56.0973; 98.9144 60.8109; 95.6718 67.8040; 98.5931 73.3443; 98.3823 84.4277; 96.7430 88.4249; 95.5038 99.2114; 0 0; 3.3909 31.4062; 3.9778 39.2143; 0.9256 44.6298];




% Randomly generate sensor energy
%initialEnergyLevels = rand(numSensors, 1) * 90 + thresholdEnergy; 
%initialEnergyLevels = [71.5344;21.8875;75.0452;19.9318;20.5744;67.6646;39.5933;68.8431;77.4218;62.4867;76.6029;31.1344;76.1462;97.3539;88.0237;17.7611;42.9793;43.2279;71.6526;63.8147;81.0428;43.0888;28.5425;17.8000;79.4741;28.5107;44.9444;59.6601;30.6058;67.7747;53.6032;23.6661;80.3739;19.0546;36.4660;31.3636;57.7785;18.2349;46.4784;19.4362;20.1056;80.5985;36.2413;64.3180;96.7980;48.9236;72.5277;78.2289;48.9378;68.9948;19.8780;94.0384;26.8715;33.9561;81.8047;53.8843;79.2062;45.6406;34.5645;13.3511;70.5965;48.6608;50.6565;64.8871;15.3463;38.4230;79.5450;72.6790;21.2799;21.7136;18.3117;10.7038;48.0798;69.0016;75.0630;57.8088;19.7936;66.8590;21.3850;22.0873;18.8735;22.7825;25.1426;27.6624;38.5732;38.4786;29.5807;32.5938;90.3630;73.2901;60.0164;26.5990;29.0828;16.9612;92.2420;73.6044;60.2010;38.2086;24.9583;66.0248;98.9141;25.3389;33.2013;45.7119;16.6595;71.5686;46.2149;98.4552;46.1966;65.8605;23.8933;44.3211;24.5021;78.2301;88.4000;41.5699;71.6982;36.4734;57.7566;84.9181;63.7741;40.1780;36.9303;50.7333;48.0381;42.3646;60.2487;76.8291;48.1901;48.6420;21.2385;12.1991;36.1167;38.5769;68.8321;96.1242;94.2158;51.2098;31.6431;78.7508;78.3395;76.6583;76.9320;19.5328;71.3404;51.6935;29.0947;18.8667;84.1217;25.7509;24.7213;69.9388;90.4950;56.4902;73.2432;23.8231;95.8111;58.6796;71.1761;13.2907;82.8283;77.3757;20.8168;57.2541;39.3250;59.1804;45.8993;47.3584;26.2664;32.9848;11.8482;93.1308;68.8330;93.9352;24.7161;92.8988;81.5192;61.9655;49.6032;33.1852;77.6752;30.5803;15.7768;79.0597;70.4082;74.3691;67.7855;47.7143;45.1686;83.4526;38.5685;83.3086;81.0166;86.7038;55.5073;67.2095;95.5805;49.9568;15.4017;88.0075];




% Initialize particles for optimization
particles = rand(numParticles, numSensors) > 0.5;

% cognition influence factor
c1 = 2; 
% social factor
c2 = 2; 
%inertia weight
w = 0.2; 

% velocity intialization
velocity = zeros(numParticles, numSensors);



time = 0;
allSensorsDepleted = false;


% loop runs until network is alive
while ~allSensorsDepleted
    
    
    
    % after every time initialize the new personal and global best
    
               particles = rand(numParticles, numSensors) > 0.5;
               personalBest = particles;
               for i = 1:numParticles
                    personalBestFitness(i) = fitnessFunction(particles(i, :), initialEnergyLevels, weights,sensorLocations,areaSize,sensing_range,numSensors,energyConsumptionRate,standbyEnergyConsumptionRate,transmissionRange,thresholdEnergy);
               end 
               [globalBestFitness, globalBestIndex] = max(personalBestFitness);
               globalBest = personalBest(globalBestIndex, :);
    
   
    %Main loop of PSO
    for iter = 1:numIterations
        
        % in every iteration update particles and velocities as per finess
        for i = 1:numParticles
            r1 = rand(1, numSensors);
            r2 = rand(1, numSensors);

            velocity(i, :) = w * velocity(i, :) + c1 * r1 .* (personalBest(i, :) - particles(i, :)) + ...
                c2 * r2 .* (globalBest - particles(i, :));

            
            % disp(velocity(i, :))
            max_velocity = 0.5; 
            % maximum velocity constraint
            velocity(i, :) = max(min(velocity(i, :), max_velocity), -max_velocity);
            
            % update particle position
            particles(i, :) = particles(i, :) + velocity(i, :);

            
            particles(i, particles(i, :) >= 0.5) = 1;
            particles(i, particles(i, :) < 0.5) = 0;

            % Evaluate fitness of particles
            fitness = fitnessFunction(particles(i, :), initialEnergyLevels, weights,sensorLocations,areaSize,sensing_range,numSensors,energyConsumptionRate,standbyEnergyConsumptionRate,transmissionRange,thresholdEnergy);

            % update personal best
            if fitness > personalBestFitness(i)
                personalBest(i, :) = particles(i, :);
                personalBestFitness(i) = fitness;
            end
        end

        % update global best
        [iterBestFitness, iterBestIndex] = max(personalBestFitness);
        if iterBestFitness > globalBestFitness
            globalBestFitness = iterBestFitness;
            globalBest = personalBest(iterBestIndex, :);
        end
    end
    
     % Update energy levels based on active sensors
    activeSensors = globalBest; % Identify active sensors
    %disp(activeSensors);
    
    
    %initialEnergyLevels(activeSensors,:) = max(thresholdEnergy, initialEnergyLevels(activeSensors,:) - energyConsumptionRate);
    % initialEnergyLevels(~activeSensors) = max(0, initialEnergyLevels(~activeSensors) - standbyEnergyConsumptionRate); 
    
    for i = 1: numSensors
        if( globalBest(i) == 1 && (initialEnergyLevels(i) - energyConsumptionRate) >= thresholdEnergy)
            initialEnergyLevels(i) = initialEnergyLevels(i) - energyConsumptionRate;
        end
    end
   
    if all(initialEnergyLevels - energyConsumptionRate < thresholdEnergy)
        allSensorsDepleted = true;
    end
    
    disp(['Time: ', num2str(time), ' sec, Fitness: ', num2str(globalBestFitness)]);
    
    % calculate total energy and coverage
    if mod(time, 2) == 0
        energy =0;
        for i = 1:numSensors
            if ((initialEnergyLevels(i) - energyConsumptionRate) >= thresholdEnergy)
                energy= energy+ initialEnergyLevels(i);
            end
        end
        coverage = coverageCalculatorFunction(globalBest, initialEnergyLevels, weights, sensorLocations, areaSize, sensing_range, thresholdEnergy,energyConsumptionRate);
        
       numActiveSensors = 0;
        for i = 1:numel(globalBest)
            if globalBest(i) == 1 && ((initialEnergyLevels(i) - energyConsumptionRate) >= thresholdEnergy)
                 numActiveSensors = numActiveSensors + 1;
            end
        end

        
       [numComponents, sensorSets] = calculateConnectivity(sensorLocations, numSensors, transmissionRange,globalBest,initialEnergyLevels,thresholdEnergy,energyConsumptionRate);
        
        % Store values
        disp(['Number of Active Sensors : ' num2str(numActiveSensors)]);
        disp(['Coverage : ' num2str(coverage)]);
        disp(['Number of Connected components: ' num2str(numComponents)]);
        disp(['Total Residual Energy : ' num2str(energy)]);
        disp("***********************************************");

        numComponentsTime = [numComponentsTime,numComponents];
        numActiveSensorsTime = [numActiveSensorsTime, numActiveSensors];
        energySumOverTime = [energySumOverTime, energy];
        coverageOverTime = [coverageOverTime, coverage];
        timePoints = [timePoints, time];
    end
    
    % identify active and stand by sensor for plot
    activeSensorsIndices = [];
    standbySensorsIndices = [];
    for i = 1:numel(globalBest)
        if globalBest(i) == 1 && ((initialEnergyLevels(i) - energyConsumptionRate) >= thresholdEnergy)
            activeSensorsIndices = [activeSensorsIndices , i];
        end
    end

    for i = 1:numel(globalBest)
        if globalBest(i) == 0 && ((initialEnergyLevels(i) - energyConsumptionRate) >= thresholdEnergy)
            standbySensorsIndices = [standbySensorsIndices , i];
        end
    end

    %disp(numel(activeSensorsIndices));
    
    %figure(1);
    %hold off;
    %axis equal;
    %plot
    %scatter(sensorLocations(activeSensorsIndices, 1), sensorLocations(activeSensorsIndices, 2), 25, 'red', 'filled');
    %hold on;
    %scatter(sensorLocations(standbySensorsIndices, 1), sensorLocations(standbySensorsIndices, 2), 25, 'blue', 'filled');
   
     %for i = 1:size(sensorLocations, 1)
     %    viscircles(sensorLocations(i, :), sensing_range, 'EdgeColor', 'r');
     %end
     %viscircles(sensorLocations(activeSensorsIndices, :), sensing_range, 'EdgeColor', 'r','LineWidth',1);
   
    
    %xlabel('X');
    %ylabel('Y');
    %legend('Active Sensors', 'Standby Sensors');
    %title(['Sensor Network (Time: ', num2str(time), ' min)']);
    %axis([0 areaSize 0 areaSize]);
    %xticks(0:20:100);
    %drawnow;
    

    PlotSensorNetwork(sensorLocations, activeSensorsIndices, standbySensorsIndices);
    time = time + 2; 
    pause(2);
end

finalTime = time;


% Plot energy sum and coverage over time
figure;
%subplot(4, 1, 1);
plot(timePoints, energySumOverTime,'LineWidth',3);
xlabel('Time (min)');
ylabel('Total Residual Energy (J)');
title('Residual Energy vs Time');

figure;
%subplot(4, 1, 2);
plot(timePoints, coverageOverTime,'LineWidth',3);
xlabel('Time (min)');
ylabel('Coverage (m*m)');
title('Coverage vs Time');

figure;
%subplot(4, 1, 3);
plot(timePoints, numActiveSensorsTime,'LineWidth',3);
xlabel('Time (min)');
ylabel('Number of Active Sensors');
title('Active Sensors vs Time');


%subplot(4, 1, 4);
figure;
plot(timePoints, numComponentsTime,'LineWidth',3);
xlabel('Time (min)');
ylabel('Number of Clusters');
title('Clusters vs Time');

TotalCoverageArea = sum(coverageOverTime);
disp(['Total Coverage : ', num2str(TotalCoverageArea)]);
disp(['Residual Energy : ', num2str(energySumOverTime)]);
disp(['Coverage : ', num2str(coverageOverTime)]);
disp(['Active Sensors : ', num2str(numActiveSensorsTime)]);


% fitness function
function fitness = fitnessFunction(particle, energyLevels, weights,sensorLocations,areaSize,sensing_range,numSensors,energyConsumptionRate,standbyEnergyConsumptionRate,transmissionRange,thresholdEnergy)
    
    coverageArea = coverageCalculatorFunction(particle, energyLevels, weights,sensorLocations,areaSize,sensing_range,thresholdEnergy,energyConsumptionRate);
    %disp(coverageArea);
    coverage = (coverageArea * 100) / (areaSize^2);
    %disp(coverage);

    num_active_sensor = 0;
    total_residual_energy = 0;
    for i = 1:numel(particle(1,:))
        if particle(1,i) == 1 && ((energyLevels(i) - energyConsumptionRate) >= thresholdEnergy)
            num_active_sensor = num_active_sensor + 1;
            total_residual_energy =total_residual_energy + energyLevels(i);
        end
    end

    numTotalSensors = 0;
     for i = 1:numel(energyLevels)
        if ((energyLevels(i) - energyConsumptionRate) >= thresholdEnergy)
            numTotalSensors = numTotalSensors + 1;
        end
    end


    activePercentage = (num_active_sensor * 100)/ numSensors; 
    %lifetime_fitness = total_residual_energy / (energyConsumptionRate*num_active_sensor);
    avgBatteryLifeFitness = calculateAvgBatteryLifeFitness(energyLevels,particle,energyConsumptionRate,thresholdEnergy);
      %disp(avgBatteryLifeFitness);


    [numSets, sensorSets] = calculateConnectivity(sensorLocations, numSensors, transmissionRange,particle,energyLevels,thresholdEnergy,energyConsumptionRate);
    numSetsFitness = (numSets*100)/numSensors;
    %coveragePercentage = coverageCalculater(num_active_sensor, areaSize,sensing_range);
    %minBatteryLevel = minBatteryLevelCalculate(energyLevels,particle,energyConsumptionRate,thresholdEnergy);
  
     % Define scaling factors for each term
    coverageScale = 0.4;
    activePercentageScale = 0.3;
    avgBatteryLifeScale = 0.1;
    numSetsScale = 0.1;
    %disp([ num2str(weights(1) * coverage * coverageScale) , " ",num2str(weights(3)*avgBatteryLifeFitness*avgBatteryLifeScale)," ",num2str(weights(2) * activePercentage*activePercentageScale)," ",num2str(weights(4) * numSetsFitness*numSetsScale)]);

    % fitness calculation based on above parameters
    
    answer = 0;
    if activePercentage >= 30
        answer = (weights(1) * coverage * coverageScale) - (weights(2) * activePercentage*activePercentageScale) + (weights(3)*avgBatteryLifeFitness*avgBatteryLifeScale) - (weights(4) * numSetsFitness*numSetsScale)  ;
        %answer = (weights(1) * coverage * coverageScale) + (weights(2) * lifetime_fitness*1000) - (weights(4) * numSetsFitness*numSetsScale);
    else
        answer = coverage ;
    end
    fitness = answer;
    %disp([fitness, " ", answer, " ", activePercentage]);
end


% coverage calculation
function coverageArea = coverageCalculatorFunction(particle, energyLevels, weights,sensorLocations,areaSize,sensing_range,thresholdEnergy,energyConsumptionRate)
    sensors = particle;
    range = sensing_range;
    gridLength = areaSize;

    % Initialize a grid to represent the coverage area
    coverageGrid = zeros(gridLength);

    % Iterate over each active sensor
    activeSensorIndices = find(sensors == 1);
    %if(length(activeSensorIndices) > 0)
        for i = 1:length(activeSensorIndices)
            sensorIndex = activeSensorIndices(i);
            sensorLocation = sensorLocations(sensorIndex, :);
            % Check if the sensor has sufficient battery level to be active
            if ((energyLevels(sensorIndex) - energyConsumptionRate) >= thresholdEnergy)
                % Mark the cells in the coverage grid that the sensor can sense
                minX = max(1, floor(sensorLocation(1) - range));
                maxX = min(gridLength, ceil(sensorLocation(1) + range));
                minY = max(1, floor(sensorLocation(2) - range));
                maxY = min(gridLength, ceil(sensorLocation(2) + range));
        
                for x = minX:maxX
                    for y = minY:maxY
                        % Check if the cell is within the sensor's range
                        if norm([x, y] - sensorLocation) <= range
                            coverageGrid(x, y) = 1;
                        end
                    end
                end
            end
        end
    %end


    % Calculate the coverage area by counting the distinct 1*1 square
    coverageArea = sum(coverageGrid(:)); 
end




% average battery life of active sensors
function avgBatteryLifeFitness = calculateAvgBatteryLifeFitness(energyLevels,particle,energyConsumptionRate,thresholdEnergy)

    energyLevelOfActiveSensors = [];
    for i = 1:numel(particle(1,:))
        if particle(1,i) == 1 && ((energyLevels(i) - energyConsumptionRate) >= thresholdEnergy)
           energyLevelOfActiveSensors = [energyLevelOfActiveSensors,energyLevels(i)];
        end
    end
    avgBatteryLife = mean(energyLevelOfActiveSensors);
    avgBatteryLifeFitness = avgBatteryLife;

end

function minBatteryLevel = minBatteryLevelCalculate(energyLevels,particle,energyConsumptionRate,thresholdEnergy)

    energyLevelOfActiveSensors = [];
    for i = 1:numel(particle(1,:))
        if particle(1,i) == 1 && ((energyLevels(i) - energyConsumptionRate) >= thresholdEnergy)
           energyLevelOfActiveSensors = [energyLevelOfActiveSensors,energyLevels(i)];
        end
    end

    minBatteryLevel = min(energyLevelOfActiveSensors);
end 

function coveragePercentage = coverageCalculater(activeSensors, area,sensingRadius)
    %area = 1 - e^(-lamda * pie * r^2)
    %lamda = no of active sensor / total geographic area

    lambda = activeSensors / (area*area);
    pi_value = pi;
    r = sensingRadius;

    coveragePercentage = (1 - exp(-lambda * pi_value * r^2))*100;

end

% connected sets calculation
function [numSets, sensorSets] = calculateConnectivity(sensorLocations, numSensors, transmissionRange, particle, energyLevels,thresholdEnergy,energyConsumptionRate)
    sensorSets = {};
    numSets = 0;
    
    
    visited = zeros(numSensors, 1);
    % BFS based approach 
    for i = 1:numSensors
        if visited(i) == 0 && particle(i) == 1 && ((energyLevels(i) - energyConsumptionRate) >= thresholdEnergy)
            
            numSets = numSets + 1;
            currentSet = [i];
            visited(i) = 1;
            
            queue = [i];  % Initialize the queue with the current sensor
            
            while ~isempty(queue)
                currentSensor = queue(1);
                queue(1) = [];  % Dequeue
                
                % Check neighbors within transmission range
                for j = 1:numSensors
                    if i ~= j && particle(j) == 1 && ((energyLevels(i) - energyConsumptionRate) >= thresholdEnergy) && norm(sensorLocations(currentSensor, :) - sensorLocations(j, :)) <= transmissionRange && visited(j) == 0
                        currentSet = [currentSet, j];
                        visited(j) = 1;
                        queue = [queue, j];  % Enqueue the neighbor
                    end
                end
            end
            
            sensorSets{numSets} = currentSet;
        end
    end
end



function PlotSensorNetwork(sensorLocations, activeSensorsIndices, standbySensorsIndices)
    % Update the existing figure
    figure(1);
    cla;

    hold on;
    axis equal;
    
    scatter(sensorLocations(activeSensorsIndices, 1), sensorLocations(activeSensorsIndices, 2), 25, 'red', 'filled');
    scatter(sensorLocations(standbySensorsIndices, 1), sensorLocations(standbySensorsIndices, 2), 25, 'blue', 'filled');
    viscircles(sensorLocations(activeSensorsIndices, :), 6, 'EdgeColor', 'r','LineWidth',1);
    axis([0 100 0 100]);
    title('Sensor Network');
    xlabel('X meters');
    ylabel('Y meters');
    legend('Active Sensors', 'Stand by Sensors');
    hold off;
    pause(2);
    
end
