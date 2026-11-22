% =========================================================================
%% COMPARISON OF SOLUBLE CO2 FLUX BETWEEN ANALYTICAL AND MULTI-PHASE Kr UPSCALING
% =========================================================================

% DESCRIPTION:
%   This code evalautes percentage of CO2 flux for both exact analytical
%   solution, and immiscible multi-phase upscaling (conservation of phases flux)

% INPUTS
% This code takes results of analytical solution. The user is required to
% change directory and supply .mat file that contains the saved results

close all;
clear;
clc;
warning off;
% =========================================================================
%% PARSE INPUTS
% =========================================================================
% .mat file written by MultiLayer.m (variable results_file)
results_file = 'Otway_miscible_results.mat';
load(results_file, 'combined', 'solute_solvent_volume_ratio', 'pore_volume_at_breakthrough');


% =========================================================================
%% ESTIMATION OF SOLUBLE FLUX PERCENTAGE
% =========================================================================

% Estimating soluble flux before correction
q_soluble_no_alpha_correction = combined.F_w .* solute_solvent_volume_ratio;

% Estiamting soluble flux after correction
q_soluble_analytical = combined.alpha_factor .* q_soluble_no_alpha_correction;

% Estimating flux percentage
% Free flux is 1 - Fw
q_soluble_no_alpha_correction_ratio = q_soluble_no_alpha_correction ./ (1 - combined.F_w);
q_soluble_analytical_ratio = combined.alpha_factor .* q_soluble_no_alpha_correction_ratio;

% =========================================================================
% PLOTTING
% =========================================================================

figure('Color', 'w', 'Name', 'Dissolved to free CO_2 flux ratio', 'Position', [100, 100, 850, 600]);
hold on;
plot(combined.injected_pv, q_soluble_no_alpha_correction_ratio, 'b-', 'LineWidth', 3, 'DisplayName', 'No Correction');
plot(combined.injected_pv, q_soluble_analytical_ratio, 'r-', 'LineWidth', 3, 'DisplayName', '\alpha Correction');
xline(pore_volume_at_breakthrough, '--k');
xlabel('PVI', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('CO2 dissolved to free flux ratio', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 2);
grid on; box on; xlim([0,1.5]);
set(gca,'FontSize',20,'FontWeight','bold','linewidth',1.5); set(gcf,'color','w');