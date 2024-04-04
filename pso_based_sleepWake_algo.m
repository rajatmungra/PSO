numSensors = 200;                   % Number of sensors
areaSize = 100;                     % Size of the area
numIterations = 5;                  % Number of iterations
numParticles = 20;                  % Number of particles
weights = [0.5, 0.2, 0.15, 0.15];   % Weights for fitness function
energyConsumptionRate = 5;         % Energy consumption rate
standbyEnergyConsumptionRate = 2;   % Standby energy consumption rate
sensing_range = 10;                 % Sensing range of sensors
energySumOverTime = [];             % Array to store energy over time
coverageOverTime = [];              % Array to store coverage over time
timePoints = [];                    % Array to store time points
transmissionRange = 7;              % Transmission range of sensors


% Randomly generate sensor locations
sensorLocations = rand(numSensors, 2) * areaSize;

% Randomly generate sensor energy
initialEnergyLevels = rand(numSensors, 1) * 100; 

% Initialize particles for optimization
particles = rand(numParticles, numSensors) > 0.5;

% cognition influence factor
c1 = 2; 
% social factor
c2 = 2; 
%inertia weight
w = 0.5; 

% velocity intialization
velocity = zeros(numParticles, numSensors);
personalBest = particles;
personalBestFitness = zeros(numParticles, 1);


% initial fitness evaluation
for i = 1:numParticles
    personalBestFitness(i) = fitnessFunction(particles(i, :), initialEnergyLevels, weights,sensorLocations,areaSize,sensing_range,numSensors,energyConsumptionRate,standbyEnergyConsumptionRate,transmissionRange);
end

% initializing global best
[globalBestFitness, globalBestIndex] = max(personalBestFitness);
globalBest = personalBest(globalBestIndex, :);


time = 0;
allSensorsDepleted = false;


% loop runs until network is alive
while ~allSensorsDepleted
    
    time = time + 2; 
    
    % after every time initialize the new personal and global best
    if(time > 2)
               particles = rand(numParticles, numSensors) > 0.5;
               personalBest = particles;
               for i = 1:numParticles
                    personalBestFitness(i) = fitnessFunction(particles(i, :), initialEnergyLevels, weights,sensorLocations,areaSize,sensing_range,numSensors,energyConsumptionRate,standbyEnergyConsumptionRate,transmissionRange);
               end 
               [globalBestFitness, globalBestIndex] = max(personalBestFitness);
               globalBest = personalBest(globalBestIndex, :);
    end
   
    %Main loop of PSO
    for iter = 1:numIterations
        
        % in every iteration update particles and velocities as per finess
        for i = 1:numParticles
            r1 = rand(1, numSensors);
            r2 = rand(1, numSensors);

            velocity(i, :) = w * velocity(i, :) + c1 * r1 .* (personalBest(i, :) - particles(i, :)) + ...
                c2 * r2 .* (globalBest - particles(i, :));

            
            % disp(velocity(i, :))
            max_velocity = 0.25; 
            % maximum velocity constraint
            velocity(i, :) = max(min(velocity(i, :), max_velocity), -max_velocity);
            
            % update particle position
            particles(i, :) = particles(i, :) + velocity(i, :);

            
            particles(i, particles(i, :) >= 0.5) = 1;
            particles(i, particles(i, :) < 0.5) = 0;

            % Evaluate fitness of particles
            fitness = fitnessFunction(particles(i, :), initialEnergyLevels, weights,sensorLocations,areaSize,sensing_range,numSensors,energyConsumptionRate,standbyEnergyConsumptionRate,transmissionRange);

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
    initialEnergyLevels(activeSensors,:) = max(0, initialEnergyLevels(activeSensors,:) - energyConsumptionRate); 
    % initialEnergyLevels(~activeSensors) = max(0, initialEnergyLevels(~activeSensors) - standbyEnergyConsumptionRate); 
    
   
    if all(initialEnergyLevels <= 0)
        allSensorsDepleted = true;
    end
    
    disp(['Time: ', num2str(time), ' sec, Fitness: ', num2str(globalBestFitness)]);

    % calculate total energy and coverage
    if mod(time, 2) == 0
        energy =0;
        for i = 1:numParticles
            energy= energy+ initialEnergyLevels(i);
        end
        coverage = coverageCalculatorFunction(globalBest, initialEnergyLevels, weights, sensorLocations, areaSize, sensing_range, numSensors, energyConsumptionRate, standbyEnergyConsumptionRate);
        
        % Store values
        energySumOverTime = [energySumOverTime, energy];
        coverageOverTime = [coverageOverTime, coverage];
        timePoints = [timePoints, time];
    end
    
    % identify active and stand by sensor for plot
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
    
    %plot
    scatter(sensorLocations(activeSensorsIndices, 1), sensorLocations(activeSensorsIndices, 2), 100, 'red', 'filled');
    hold on;
    scatter(sensorLocations(standbySensorsIndices, 1), sensorLocations(standbySensorsIndices, 2), 100, 'blue', 'filled');
   
     
   
    
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


% fitness function
function fitness = fitnessFunction(particle, energyLevels, weights,sensorLocations,areaSize,sensing_range,numSensors,energyConsumptionRate,standbyEnergyConsumptionRate,transmissionRange)
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
    coverage = (coverageArea * 100) / (gridLength^2);
    %disp(coverage);
    active_sensors =  logical(particle);
    active_sensor_positions = sensorLocations(active_sensors, :);
    
    num_active_sensor = size(active_sensor_positions(:,1));
    energy = (num_active_sensor(1,1) * 100)/ numSensors; 

    avgBatteryLifeFitness = calculateAvgBatteryLifeFitness(energyLevels,particle);
      %disp(avgBatteryLifeFitness);


    [numSets, sensorSets] = calculateConnectivity(sensorLocations, numSensors, transmissionRange,particle,energyLevels);
    
    % fitness calculation based on above parameters
    fitness = (weights(1) * coverage) - (weights(2) * energy) + (weights(3)* avgBatteryLifeFitness) + (weights(4) * numSets*100/numSensors)  ;
end


% coverage calculation
function coverageArea = coverageCalculatorFunction(particle, energyLevels, weights,sensorLocations,areaSize,sensing_range,numSensors,energyConsumptionRate,standbyEnergyConsumptionRate)
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




% average battery life of active sensors
function avgBatteryLifeFitness = calculateAvgBatteryLifeFitness(energyLevels,particle)

    initialEnergyLevels = energyLevels;
    activeSensorIndices = find(particle == 1);
    activeBatteryLevels = initialEnergyLevels(activeSensorIndices);
    avgBatteryLife = mean(activeBatteryLevels);
    avgBatteryLifeFitness = avgBatteryLife;

end


% connected sets calculation
function [numSets, sensorSets] = calculateConnectivity(sensorLocations, numSensors, transmissionRange, particle, energyLevels)
    sensorSets = {};
    numSets = 0;
    
    
    visited = zeros(numSensors, 1);
    % BFS based approach 
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