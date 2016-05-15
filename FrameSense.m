function SensorPosition = FrameSense(MeasurementMatrix,sensorNO)

%% ************************  Introduction ********************************
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
% Reference:         J. Ranieri et al.(2014) Near-Optimal sensor placement for
%                    linear inverse problems, IEEE-TSP.
%
% AUTHOR Information:
% Jiang Chaoyang, EEE, NTU 
% Email: cjiang003@e.ntu.edu.sg, chaoyangjiang@hotmail.com
% Finished at 10-sep-2014 
% Latest version at 15-OCT-2015

%% ********************* Preparation *****************************
[N,n] = size(MeasurementMatrix);  
if (n > N); 
  error('Rows is not enough'); 
end
if (sensorNO < n); 
 error(' More sensors are needed'); 
end


%% ************* Initialization ***********************************
Nindex = 1:1:N;
Lindex = zeros(1,N-sensorNO); % The index of removed rows

%% ****************** Calculate the square Matrix ***************** 

Psi2 = (MeasurementMatrix*MeasurementMatrix').^2; % Psi2 = (Phi*Phi^T).^2

%% ******************** First two deleted sensor locations *************
[C,IndexTemp] = max(Psi2); 
[Maxvalue, Lindex(1)] = max(C);
Lindex(2) = IndexTemp (Lindex(1));

Psi2(Lindex(1),:)=zeros(1,N);
Psi2(:,Lindex(1))=zeros(N,1);
Psi2(Lindex(2),:)=zeros(1,N);
Psi2(:,Lindex(2))=zeros(N,1);

%% ** Determine Lindex (Remove the other rows one-by-one)*****************
for i = 3:N-sensorNO
    OneRowInfluence = 2*sum(Psi2,2)-diag(Psi2);
    [Maxvalue, maxindex] = max(OneRowInfluence);
    Lindex(i)=maxindex;    

    Psi2(Lindex(i),:)=zeros(1,N);
    Psi2(:,Lindex(i))=zeros(N,1);
end
SensorPosition = setdiff(Nindex,Lindex);
