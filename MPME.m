function SensorPosition = MPME(MeasurementMatrix,sensorNO)

%% ************************  Introduction ********************************
% MPME: maximal projection on minimum eigenspace
%
% For detail about this algorithm, See:
%
%     C. Jiang, Y. C. Soh and H. Li, Sensor placement by maximal projection on
%     minimum eigenspace for linear inverse problems, 2015. (submitted to IEEE-
%     TSP).
%
% INPUT:
% MeasurementMatrix: a matrix from which we need to choose 'sensorNO' rows to
%                    construct a new matrix, and the chosen row index corresponds 
%                    to the sensor positions which are saved in 'SensorPosition'.
% sensorNO:          sensor number
%
%
% OUTPUT:
% SensorPosition:    a vector whose elements are the indices of the rows
%                    of 'MeasurementMatrix' chosen to place sensor
%

%
% AUTHOR Information:
% Jiang Chaoyang, EEE, NTU 
% Email: cjiang003@e.ntu.edu.sg, chaoyangjiang@hotmail.com
% Finished at 10-sep-zo14 
% Modified at 15-Nov-2015
%% ********************* Preparation *****************************
[N,n] = size(MeasurementMatrix);  
if (n > N); 
 error('Rows is not enough'); 
end
if (sensorNO < n); 
 error('More sensors are needed'); 
end

%% ******************* Initilization *******************************
SensorPosition = zeros(sensorNO,1);

%% ************* Determine the first to n-th sensor position ********************
    
%************* Determine the 1st sensor position ********************
RowNorm = sum(MeasurementMatrix.^2,2);
[Maxvalue, SensorPosition(1)]=max(RowNorm);
%*********** determine 2nd to n-th sensor position ****************
Phi= MeasurementMatrix(SensorPosition(1),:);
OrthBasis = orth(Phi');
P = eye(n,n)-OrthBasis*OrthBasis';
MeasurementMatrix(SensorPosition(1),:)=zeros(1,n); 
for i = 2:n
    ProjectionTemp = P*MeasurementMatrix';
    ProjectionTempNorm = sum(ProjectionTemp.^2);
    [Maxvalue, SensorPosition(i)]=max(ProjectionTempNorm);  
    
    Phi=[Phi;MeasurementMatrix(SensorPosition(i),:)];
    OrthBasis = orth(Phi');
    P = eye(n,n)-OrthBasis*OrthBasis';    
    MeasurementMatrix(SensorPosition(i),:)=zeros(1,n);
end
%% ************* Determine the reminder sensor position ********************
if sensorNO>n    
    for i = n+1:sensorNO
    [S,V,FullOrthBasis]=svd(Phi'*Phi);
    P = FullOrthBasis(:,n)*FullOrthBasis(:,n)';
        
    ProjectionTemp = P*MeasurementMatrix';
    ProjectionTempNorm = sum(ProjectionTemp.^2);
    [Maxvalue, SensorPosition(i)]=max(ProjectionTempNorm);  
    
    Phi=[Phi;MeasurementMatrix(SensorPosition(i),:)];
    MeasurementMatrix(SensorPosition(i),:)=zeros(1,n);    
    end
end    
