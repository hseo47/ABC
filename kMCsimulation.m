function runkMCsimulation()

    % Assumptions:                                                                           
    %   - Bacteria: E_bacteria = 100 MPa, nu_bacteria = 0.3.
    % 
    % Experimental data:
    %   - Rigid substrate: E_substrate = 84.5 MPa.
    %   - Hydro-softened substrate:  E_substrate = 1.8 MPa.
    % 
    % The composite effective moduli:
    %   1/E_eff = (1 - nu_b^2)/E_bacteria + (1 - nu_s^2)/E_substrate.
    %
    % Based on our fitted-paramteres from simulation
    %   - For rigid: delta_gamma = Jm^-2 gives U_threshold = -1.6×10e-12 J.
    %   - For hydro-softened:  delta_gamma = 1.49 Jm^-2 gives U_threshold = -1.25×10e-12 J.
    %
    % The simulation uses a logistic growth model (timePoints: 6, 12, 18, 24 hrs).
    % For each particle, random DLVO parameters are drawn, and adhesion is decided
    % based on whether the computed xDLVO energy is below the substrate-specific threshold.
    
    %% 1. Growth and Simulation Parameters
    N0 = 0.3;         % Initial optical density
    r_growth = 0.4;   % Growth rate (hr^-1)
    K = 2.0;          % Carrying capacity (OD)
    logisticFun = @(t) K ./ (1 + ((K-N0)/N0).*exp(-r_growth*t));
    
    timePoints = [6, 12, 18, 24];  
    totalParticles = 1000;         
    
    % Simulation domain (binary image mask)
    planeWidth = 1024;
    planeHeight = 1024;
    simMask_rigid = false(planeHeight, planeWidth);
    simMask_soft  = false(planeHeight, planeWidth);
    
    %% 2. Material Parameters
    E_bacteria = 100e6;   
    nu_bacteria = 0.3;
    
    % Substrate properties:
    % Rigid substrate:
    E_substrate_rigid = 84.5e6;   % 84.5 MPa
    % Hydro-softened substrate:
    E_substrate_soft = 1.8e6;   % 1.8 MPa
    nu_substrate = 0.3;         % assumed same for both
    
    % Effective contact radius (R) set to 1 µm:
    R_eff = 1e-6; % 1 µm
    
    DeltaGamma_rigid = 6.43; 
    DeltaGamma_soft  = 1.49; 
    
    % Compute composite effective moduli:
    E_eff_rigid = computeEffectiveModulus(E_bacteria, nu_bacteria, E_substrate_rigid, nu_substrate);
    E_eff_soft  = computeEffectiveModulus(E_bacteria, nu_bacteria, E_substrate_soft, nu_substrate);
    
    % Calculate Griffith thresholds:
    thresholdRigid = calculateGriffithThreshold(E_eff_rigid, R_eff, DeltaGamma_rigid);
    thresholdSoft  = calculateGriffithThreshold(E_eff_soft, R_eff, DeltaGamma_soft);
    fprintf('Computed Rigid Threshold: %.2e J\n', thresholdRigid);
    fprintf('Computed Soft  Threshold: %.2e J\n', thresholdSoft);
    
    %% 3. xDLVO Stochastic Parameter Ranges
    HamakerRange = [0.5e-20, 5e-20]; % in J
    Psi0Range    = [10e-3, 50e-3]; % in V
    IonicStrengthRange = [10e-3, 100e-3]; % in mol/L
    HydrationRange     = [1e-20, 1e-19]; % in J
    separation = 5e-9; % 5 nm
    
    %% 4. Figure Windows for Visualization
    figRigid = figure('Name','Rigid Substrate Adhesion');
    figSoft  = figure('Name','Soft Substrate Adhesion');
    
    backgroundColor = [209, 209, 209] / 255;  
    bacteriaColor   = [188, 124, 124] / 255;  
    noiseVariance = 0.005;  
    
    %% 5. KMC Simulation Loop Over Time Points
    totalAdheredRigid = 0;
    totalAdheredSoft  = 0;
    
    for t = timePoints
        fractionNow = logisticFun(t) / K;  
        numParticles = round(fractionNow * totalParticles);
        
        for p = 1:numParticles
            % --- Rigid Substrate ---
            [AH_r, Psi0_r, Ion_r, Hyd_r] = randomDLVOparams(HamakerRange, Psi0Range, IonicStrengthRange, HydrationRange);
            xE_rigid = calcXDLVOEnergy(AH_r, Psi0_r, Ion_r, Hyd_r, separation);
            if xE_rigid <= thresholdRigid
                totalAdheredRigid = totalAdheredRigid + 1;
                simMask_rigid = addParticleToMask(simMask_rigid);
            end
            
            % --- Hydro-softened Substrate ---
            [AH_s, Psi0_s, Ion_s, Hyd_s] = randomDLVOparams(HamakerRange, Psi0Range, IonicStrengthRange, HydrationRange);
            xE_soft = calcXDLVOEnergy(AH_s, Psi0_s, Ion_s, Hyd_s, separation);
            if xE_soft <= thresholdSoft
                totalAdheredSoft = totalAdheredSoft + 1;
                simMask_soft = addParticleToMask(simMask_soft);
            end
        end
        
        coloredRigid = convertMaskToColoredImage(simMask_rigid, backgroundColor, bacteriaColor);
        coloredRigidNoisy = imnoise(coloredRigid, 'gaussian', 0, noiseVariance);
        
        coloredSoft = convertMaskToColoredImage(simMask_soft, backgroundColor, bacteriaColor);
        coloredSoftNoisy = imnoise(coloredSoft, 'gaussian', 0, noiseVariance);
        
        figure(figRigid);
        imshow(coloredRigidNoisy);
        title(['Rigid Substrate Coverage at Time = ', num2str(t), ' hrs']);
        imwrite(coloredRigidNoisy, sprintf('Rigid_Adhesion_at_%02dhrs.png', t));
        
        figure(figSoft);
        imshow(coloredSoftNoisy);
        title(['Soft Substrate Coverage at Time = ', num2str(t), ' hrs']);
        imwrite(coloredSoftNoisy, sprintf('Soft_Adhesion_at_%02dhrs.png', t));
        
        fprintf('Time = %d hrs: Rigid adhered = %d, Soft adhered = %d\n', t, totalAdheredRigid, totalAdheredSoft);
        pause(1);
    end
    
    %% 6. Final Report
    ratio = (totalAdheredSoft + 1) / (totalAdheredRigid + 1);
    fprintf('\n=== Final Results ===\n');
    fprintf('Total Rigid Adhered: %d\n', totalAdheredRigid);
    fprintf('Total Soft  Adhered: %d\n', totalAdheredSoft);
    fprintf('Final Adhesion Ratio (Soft/Rigid) = %.2f\n', ratio);
end

%% computeEffectiveModulus
function E_effective = computeEffectiveModulus(E_bacteria, nu_bacteria, E_substrate, nu_substrate)
    E_effective = 1 / ((1 - nu_bacteria^2)/E_bacteria + (1 - nu_substrate^2)/E_substrate);
end

%% calculateGriffithThreshold
function U_threshold = calculateGriffithThreshold(E_eff, R, DeltaGamma)
    U_threshold = - (R^(4/3)) * (E_eff^(-2/3)) * ((DeltaGamma)^(5/3));
end

%% addParticleToMask
function simMask = addParticleToMask(simMask)
    [maskHeight, maskWidth] = size(simMask);
    cx = randi(maskWidth);
    cy = randi(maskHeight);
    radius = 5;  
    [xx, yy] = meshgrid(1:maskWidth, 1:maskHeight);
    circleMask = (xx - cx).^2 + (yy - cy).^2 <= radius^2;
    simMask = simMask | circleMask;
end

%% randomDLVOparams
function [AH, Psi0, IonicStrength, Hydration] = randomDLVOparams(HamakerRange, Psi0Range, IonicStrengthRange, HydrationRange)
    AH = HamakerRange(1) + rand()*(HamakerRange(2) - HamakerRange(1));
    Psi0 = Psi0Range(1) + rand()*(Psi0Range(2) - Psi0Range(1));
    IonicStrength = IonicStrengthRange(1) + rand()*(IonicStrengthRange(2) - IonicStrengthRange(1));
    Hydration = HydrationRange(1) + rand()*(HydrationRange(2) - HydrationRange(1));
end

%% calcXDLVOEnergy
function xDLVO_Energy = calcXDLVOEnergy(AH, Psi0, IonicStrength, Hydration, r)
    epsilon0 = 8.854e-12;  % F/m
    epsilonR = 78.5;       % dimensionless
    kB = 1.38064852e-23;   % J/K
    T = 298;               % K
    NA = 6.022e23;         % 1/mol
    e = 1.60217662e-19;    % C
    
    I_m3 = IonicStrength * 1000;
    
    kappa = sqrt((2 * NA * I_m3 * e^2) / (epsilon0 * epsilonR * kB * T));
    
    U_el = (64 * epsilon0 * epsilonR * kB * T / kappa) * ...
           tanh((e * Psi0)/(4*kB*T))^2 * exp(-kappa*r);
    
    U_vdw = -AH / (6 * r);
    
    decayLength = 1e-9;
    U_hydration = Hydration * exp(-r/decayLength);
    
    xDLVO_Energy = U_vdw + U_el + U_hydration;
end

%% convertMaskToColoredImage
function rgbImage = convertMaskToColoredImage(mask, backgroundColor, bacteriaColor)
    [rows, cols] = size(mask);
    rgbImage = repmat(reshape(backgroundColor, [1,1,3]), rows, cols);
    
    for channel = 1:3
        temp = rgbImage(:,:,channel);
        temp(mask) = bacteriaColor(channel);
        rgbImage(:,:,channel) = temp;
    end
end
