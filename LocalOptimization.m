function MSE = LocalOptimization(MeasurementMatrix,SubOptSPosition,iterationNO)

%% ************************  Introduction ********************************
% INPUT:
% MeasurementMatrix: a matrix from which we need to choose 'sensorNO' columns to
%                    construct a new matrix, and the chosen column index corresponds 
%                    to the sensor positions which are saved in 'SensorPosition'.
%
% SubOptSPosition:   a vector whose elements are the indices of the columns
%                    of 'MeasurementMatrix' chosen to place sensor, which
%                    obtained from one suboptimal algorithm.
%
% iterationNO:      The proposed algorithm is iterative algorithm. 
%                   'iterationNO' is the maximum iteration number. 
%
%
% OUTPUT:
% MSE:               MSE of the estimated paramter of the linear inverse problem.
%                    Here, MSE is used to assess the locally optimized solution 
%                    of sensor placement problem.
%
%
%Additionaly, the function can provide the other two parameters, which are
%             not used.
% SensorPosition:    The locally optimized solution of the sensor placement
%                    problem
% updataNumber:      The number of rows are updated.
%
%
%
% AUTHOR Information:
% Jiang Chaoyang, EEE, NTU 
% Email: cjiang003@e.ntu.edu.sg, chaoyangjiang@hotmail.com
% Finished at 10-sep-zo14 

%% ********************* Preparation *****************************
[n,N] = size(MeasurementMatrix);  
if (n > N); 
 error('size(MeasurementMatrix,1) > size(MeasurementMatrix,2)'); 
end

[sensorNO,default] = size(SubOptSPosition);
if (sensorNO < n); 
 error(' More sensors are needed'); 
end
updataNumber=1;
%% ***************** Initialization *****************************
SensorPosition = SubOptSPosition;
Phi = zeros(sensorNO,n);
for i = 1:sensorNO 
    Phi(i,:) = MeasurementMatrix(:,SensorPosition(i))';    
end
eigenvalue = abs(eig(Phi'*Phi));
MSE = sum(1./eigenvalue);


%% **************** Iteration algorithm ***************************
for iter = 1:iterationNO
    flag_stop_iteration=1;
    for k = 1: sensorNO  
        MSEvector = 10000000*ones(N,1);
        for i = 1:N            
            flag =1;
            for j = 1:sensorNO
                if i==SensorPosition(j) %****
                    flag = 0;
                    break;
                end
            end
            
            if flag
                Phi(k,:) = MeasurementMatrix(:,i)';
                eigenvalue = abs(eig(Phi'*Phi));
                MSEvector(i) = sum(1./eigenvalue); 
            end
        end
        [Minvalue, MinvalueIndex]=min(MSEvector);
        if Minvalue<MSE        
            MSE = Minvalue;
            SensorPosition(k) = MinvalueIndex;        
            Phi(k,:)= MeasurementMatrix(:,SensorPosition(k))';
            flag_stop_iteration =0;
        end        
    end
    if flag_stop_iteration==0
        updataNumber=updataNumber+1;
    else
        break
    end    
end