function w = ConvexOpt(MeasurementMatrix,sensorNO) 
%% ************************  Introduction ********************************
%% Toolbox <SDPT3-4.0> is needed.
%% <SDPT3-4.0>:
%% Step one:          run Installmex
%% Step two:          run startup (run both in command window)
%% The toolbox can be freely downloaded from:
%% www.math.nus.edu.sg/~mattohkc//sdpt3.html
%
% INPUT:
% MeasurementMatrix: a matrix from which we need to choose 'sensorNO' columns to
%                    construct a new matrix, and the chosen column index corresponds 
%                    to the sensor positions which are saved in 'SensorPosition'.
% sensorNO:          sensor number
%
% OUTPUT:
% w:                 a vector, the indices of whose elements correspond to the 
%                    sensing indices, and the magnitude of its element implies 
%                    the significance of the correspinding sensing location.
%
%
%Reference:         [1] S. Joshi and S. Boyd (2009) Sensor selection via convex 
%                       optimization, IEEE-TSP.
%                   [2] K.-C. Toh, et al.(1999) SDPT3-a Matlab software package
%                       for semidefinite programming, version 1.3, Optimization
%                       methods and software.
%
% AUTHOR Information:
% Jiang Chaoyang, EEE, NTU, Singapore 
% Email: cjiang003@e.ntu.edu.sg, chaoyangjiang@hotmail.com
% Finished at 10-sep-zo14  
%% ******************* Initalization ********************************
   [n,N] = size(MeasurementMatrix);  
   if (n > N); 
     error(' size(MeasurementMatrix,1) > size(MeasurementMatrix,2)'); 
   end
   
   b = zeros(N,1);
   blk{1,1} = 's'; blk{1,2} = n; 
   F = cell(1,N); 
   for k = 1:N
      F{1,k} = -MeasurementMatrix(:,k)*MeasurementMatrix(:,k)';
   end
   At(1) = svec(blk(1,:),F,1);
   C{1,1} = sparse(n,n);
 
   blk{2,1} = 'l'; blk{2,2} = N; 
   At{2,1} = -speye(N,N); 
   C{2,1} = zeros(N,1); 
   
   blk{3,1} = 'l'; blk{3,2} = N; 
   At{3,1} = speye(N,N); 
   C{3,1} = ones(N,1); 
 
   blk{4,1} = 'u'; blk{4,2} = 1;
   At{4,1} = ones(1,N); 
   C{4,1} = sensorNO; 
 
   OPTIONS.parbarrier{1,1} = 1; 
   OPTIONS.parbarrier{2,1} = 0; 
   OPTIONS.parbarrier{3,1} = 0; 
   
%% ***************** solve the convex optimization Problem ***************    
[defaultParamter1,defaultParameter2,w,defaultParamter3] = sqlp(blk,At,C,b,OPTIONS);


