

numSensors = 200;
areaSize = 100;
numIterations = 5;
numParticles = 20;
weights = [0.5, 0.15, 0.15, 0.2]; 
energyConsumptionRate = 20; 
standbyEnergyConsumptionRate = 2; 
sensing_range = 10;
energySumOverTime = [];
coverageOverTime = [];
timePoints = [];
transmissionRange = 7;


sensorLocations = rand(numSensors, 2) * areaSize;


initialEnergyLevels = rand(numSensors, 1) * 100; 


particles = rand(numParticles, numSensors) > 0.5;


c1 = 2; 
c2 = 2; 
w = 0.5; 

velocity = zeros(numParticles, numSensors);
personalBest = particles;
personalBestFitness = zeros(numParticles, 1);


for i = 1:numParticles
    minMaxValue = minMaxValueCalculate(particles,numParticles,initialEnergyLevels,sensorLocations,areaSize,sensing_range,numSensors,transmissionRange);
    personalBestFitness(i) = fitnessFunction(particles(i, :), initialEnergyLevels, weights,sensorLocations,areaSize,sensing_range,numSensors,energyConsumptionRate,standbyEnergyConsumptionRate,transmissionRange,minMaxValue);
end


[globalBestFitness, globalBestIndex] = max(personalBestFitness);
globalBest = personalBest(globalBestIndex, :);


time = 0;
allSensorsDepleted = false;

while ~allSensorsDepleted
    
    time = time + 2; 
    
    if(time > 2)
               particles = rand(numParticles, numSensors) > 0.5;
               personalBest = particles;
               for i = 1:numParticles
                    minMaxValue = minMaxValueCalculate(particles,numParticles,initialEnergyLevels,sensorLocations,areaSize,sensing_range,numSensors,transmissionRange);
                    personalBestFitness(i) = fitnessFunction(particles(i, :), initialEnergyLevels, weights,sensorLocations,areaSize,sensing_range,numSensors,energyConsumptionRate,standbyEnergyConsumptionRate,transmissionRange,minMaxValue);
               end 
               [globalBestFitness, globalBestIndex] = max(personalBestFitness);
               globalBest = personalBest(globalBestIndex, :);
    end
   
    for iter = 1:numIterations
        
        for i = 1:numParticles
            r1 = rand(1, numSensors);
            r2 = rand(1, numSensors);

            velocity(i, :) = w * velocity(i, :) + c1 * r1 .* (personalBest(i, :) - particles(i, :)) + ...
                c2 * r2 .* (globalBest - particles(i, :));

            
            % disp(velocity(i, :))
            max_velocity = 0.25; 
            velocity(i, :) = max(min(velocity(i, :), max_velocity), -max_velocity);
            
            
            particles(i, :) = particles(i, :) + velocity(i, :);

           
            particles(i, particles(i, :) >= 0.5) = 1;
            particles(i, particles(i, :) < 0.5) = 0;

           minMaxValue = minMaxValueCalculate(particles,numParticles,initialEnergyLevels,sensorLocations,areaSize,sensing_range,numSensors,transmissionRange);
            fitness = fitnessFunction(particles(i, :), initialEnergyLevels, weights,sensorLocations,areaSize,sensing_range,numSensors,energyConsumptionRate,standbyEnergyConsumptionRate,transmissionRange,minMaxValue);

            
            if fitness > personalBestFitness(i)
                personalBest(i, :) = particles(i, :);
                personalBestFitness(i) = fitness;
            end
        end

        
        [iterBestFitness, iterBestIndex] = max(personalBestFitness);
        if iterBestFitness > globalBestFitness
            globalBestFitness = iterBestFitness;
            globalBest = personalBest(iterBestIndex, :);
        end
    end
    
     
    activeSensors = globalBest; % Identify active sensors
    

    for i = 1:numel(globalBest)
        if globalBest(i) == 1
            initialEnergyLevels(i) = max(0, initialEnergyLevels(i) - energyConsumptionRate);
        end
    end
    % initialEnergyLevels(activeSensors,:) = max(0, initialEnergyLevels(activeSensors,:) - energyConsumptionRate); 
    % initialEnergyLevels(~activeSensors) = max(0, initialEnergyLevels(~activeSensors) - standbyEnergyConsumptionRate); 
    
   
    if all(initialEnergyLevels <= 0)
        allSensorsDepleted = true;
    end
    
    disp(['Time: ', num2str(time), ' sec, Fitness: ', num2str(globalBestFitness)]);

    if mod(time, 2) == 0
        energy =0;
        for i = 1:numParticles
            energy= energy+ initialEnergyLevels(i);
        end
        coverage = coverageCalculatorFunction(globalBest, initialEnergyLevels, sensorLocations, areaSize, sensing_range);
        
        % Store values
        energySumOverTime = [energySumOverTime, energy];
        coverageOverTime = [coverageOverTime, coverage];
        timePoints = [timePoints, time];
    end
    

    activeSensorsIndices = [];
    standbySensorsIndices = [];
    for i = 1:numel(globalBest)
        if globalBest(i) == 1 && initialEnergyLevels(i) > 0
            activeSensorsIndices = [activeSensorsIndices , i];
        end
    end

    for i = 1:numel(globalBest)
        if globalBest(i) == 0 && initialEnergyLevels(i) > 0
            standbySensorsIndices = [standbySensorsIndices , i];
        end
    end

    figure(1);
    hold off;
    
    
    scatter(sensorLocations(activeSensorsIndices, 1), sensorLocations(activeSensorsIndices, 2), areaSize, 'red', 'filled');
    hold on;
    scatter(sensorLocations(standbySensorsIndices, 1), sensorLocations(standbySensorsIndices, 2), areaSize, 'blue', 'filled');
   
     
   
    
    xlabel('X');
    ylabel('Y');
    legend('Active Sensors', 'Stand by Sensors');
    title(['Sensor Network (Time: ', num2str(time), ' sec)']);
    axis([0 areaSize 0 areaSize]);
    drawnow;
end

finalTime = time;


% Plot energy sum and coverage over time
figure;
subplot(2, 1, 1);
plot(timePoints, energySumOverTime, '-o');
xlabel('Time (s)');
ylabel('Energy Sum of Active Sensors');
title('Energy Sum of Sensors Over Time');

subplot(2, 1, 2);
plot(timePoints, coverageOverTime, '-o');
xlabel('Time (s)');
ylabel('Coverage (%)');
title('Coverage of Over Time');


function minMaxValue = minMaxValueCalculate(particles,numParticles,energyLevels,sensorLocations,areaSize,sensing_range,numSensors,transmissionRange)
 
    coverage = [];
    for i = 1:numParticles
        coverageArea = coverageCalculatorFunction(particles(i,:), energyLevels, sensorLocations, areaSize, sensing_range);
        coverage = [coverage,coverageArea];
    end
    [ min_coverage , max_coverage ] = bounds( coverage );
    
    energyFitness= [];
    for i = 1:numParticles
        avgBatteryLifeFitness = calculateAvgBatteryLifeFitness(energyLevels,particles(i,:));
        energyFitness = [energyFitness,avgBatteryLifeFitness];
    end
    [min_avgBatteryLifeFitness,max_avgBatteryLifeFitness] = bounds(energyFitness);

    connectivity = [];
    for i = 1:numParticles
        [numSets, sensorSets] = calculateConnectivity(sensorLocations, numSensors, transmissionRange,particles(i,:),energyLevels);
        connectivity = [connectivity,numSets];
    end
    [min_connectivity, max_connectivity] = bounds(connectivity);

    minMaxValue = [min_coverage , max_coverage ,min_avgBatteryLifeFitness,max_avgBatteryLifeFitness,min_connectivity, max_connectivity ];
end

function fitness = fitnessFunction(particle, energyLevels, weights,sensorLocations,areaSize,sensing_range,numSensors,energyConsumptionRate,standbyEnergyConsumptionRate,transmissionRange,minMaxValue)
   
    min_coverage  = minMaxValue(1);
    max_coverage = minMaxValue(2);
    min_avgBatteryLifeFitness  = minMaxValue(3);
    max_avgBatteryLifeFitness = minMaxValue(4);
    min_connectivity = minMaxValue(5);
    max_connectivity  = minMaxValue(6);

    coverageArea = coverageCalculatorFunction(particle, energyLevels,sensorLocations,areaSize,sensing_range);
    scaledCoverageArea = ((coverageArea - min_coverage)/ (max_coverage - min_coverage)) * 100;


    numActiveSensor  = 0;
    totalAliveSensor = 0;
    for i = 1:numel(particle)
        if particle(i) == 1 && energyLevels(i) > 0
            numActiveSensor = numActiveSensor + 1;
        end
        if energyLevels(i) > 0 
            totalAliveSensor = totalAliveSensor +1;
        end
    end

    

    avgBatteryLifeFitness = calculateAvgBatteryLifeFitness(energyLevels,particle);
    scaled_avgBatteryLifeFitness = ((avgBatteryLifeFitness - min_avgBatteryLifeFitness) / (max_avgBatteryLifeFitness - min_avgBatteryLifeFitness)) * 100;
  
    
    [numSets, sensorSets] = calculateConnectivity(sensorLocations, numSensors, transmissionRange, particle, energyLevels);
    scaled_numSets = ((max_connectivity - numSets) / (max_connectivity - min_connectivity)) * 100;

    scaled_activeSensor = (numActiveSensor * 100)/totalAliveSensor;
    %disp([min_coverage,max_coverage]);
    fitness = (weights(1) * scaledCoverageArea) + (weights(2) * scaled_activeSensor) + (weights(3)* scaled_avgBatteryLifeFitness) + (weights(4) * scaled_numSets);
end


function coverageArea = coverageCalculatorFunction(particle, energyLevels,sensorLocations,areaSize,sensing_range)
    sensors = particle;
    range = sensing_range;
    gridLength = areaSize;

    % Initialize a grid to represent the coverage area
    coverageGrid = zeros(gridLength);

    % Iterate over each active sensor
    activeSensorIndices = find(sensors == 1);
    
    for i = 1:length(activeSensorIndices)
        sensorIndex = activeSensorIndices(i);
        sensorLocation = sensorLocations(sensorIndex, :);
        % Check if the sensor has sufficient battery level to be active
        if energyLevels(sensorIndex) > 0
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

    % Calculate the coverage area by counting the distinct 1*1 squares
    coverageArea = sum(coverageGrid(:)); 
end




function avgBatteryLifeFitness = calculateAvgBatteryLifeFitness(energyLevels,particle)

    initialEnergyLevels = energyLevels;
    activeSensorIndices = find(particle == 1);
    activeBatteryLevels = initialEnergyLevels(activeSensorIndices);
    avgBatteryLife = mean(activeBatteryLevels);
    avgBatteryLifeFitness = avgBatteryLife;

end



function [numSets, sensorSets] = calculateConnectivity(sensorLocations, numSensors, transmissionRange, particle, energyLevels)
    sensorSets = {};
    numSets = 0;
    
    visited = zeros(numSensors, 1);
    
    for i = 1:numSensors
        if visited(i) == 0 && particle(i) ~= 0 && energyLevels(i) > 0
            
            numSets = numSets + 1;
            currentSet = [i];
            visited(i) = 1;
            
            queue = [i];  % Initialize the queue with the current sensor
            
            while ~isempty(queue)
                currentSensor = queue(1);
                queue(1) = [];  % Dequeue
                
                % Check neighbors within transmission range
                for j = 1:numSensors
                    if i ~= j && particle(j) ~= 0 && energyLevels(j) > 0 && norm(sensorLocations(currentSensor, :) - sensorLocations(j, :)) <= transmissionRange && visited(j) == 0
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