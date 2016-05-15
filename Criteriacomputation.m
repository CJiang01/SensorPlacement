function [MSE, WCEV, CondationalNO] = Criteriacomputation(sensorposition, MeasurementMatrix)
%% ************************  Introduction ********************************
% INPUT:
% MeasurementMatrix: a matrix from which we need to choose 'sensorNO' rows to
%                    construct a new matrix, and the chosen row index corresponds 
%                    to the sensor positions which are saved in 'sensorposition'.
%
% sensorposition:    a vector whose elements are the indices of the rows
%                    of 'MeasurementMatrix' chosen to place sensor
%
%
% OUTPUT:
% MSE:               The MSE of estiamted parameter when the sensing
%                    locations are determined by 'sensorposition'
%
%
%
% AUTHOR Information:
% Jiang Chaoyang, EEE, NTU 
% Email: cjiang003@e.ntu.edu.sg, chaoyangjiang@hotmail.com
% Finished at 10-sep-2014 
%% *****************************************************************
[N,n]=size(MeasurementMatrix);
if n>N
    error('Rows number is not enough')
end


[sensorNO, n1]=size(sensorposition);
if sensorNO<n1
    sensorposition=sensorposition';
end
[sensorNO, tilde_defaultParameter]=size(sensorposition);

Phi = zeros(sensorNO, n);
for i = 1:sensorNO
    Phi(i,:) = MeasurementMatrix(sensorposition(i),:);
end
eigenvalue = abs(eig(Phi'*Phi));
MSE = sum(1./eigenvalue);
WCEV = 1/min(eigenvalue);
CondationalNO = max(eigenvalue)/min(eigenvalue);