function w = SparSenSe_CVX(MeasurementMatrix,gamma) 
%% ************************  Introduction ********************************
%% CVX toolbox is needed.
%% The CVX toolbox can be freely downloaded from
%% cvxr.com/cvx
%
%% The main part of the code is contributed by Dr. H. Jamali-Rad and Dr. A. Simonetto. 
%% The SparSenSe method is presented in the follwoing paper:
%
% H. Jamali-Rad, A. Simonetto, and G. Leus, Sparsity-aware sensor selection: 
% Centralized and distributed algorithms, IEEE Signal Processing Letters, 22(12),
% pp. 217-220, 2014.
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
% AUTHOR Information:
% Jiang Chaoyang, EEE, NTU, Singapore 
% Email: cjiang003@e.ntu.edu.sg, chaoyangjiang@hotmail.com
% Finished at 15-Nov-zo15  
%% *****************************************************************
[N,n] = size(MeasurementMatrix);  
if (n > N); 
 error('Rows is not enough'); 
end
%% SparSenSe (Centralized Algorithm) 
tic
cvx_begin sdp quiet
    variable w(N)
    variable u(n)
    dual variable y{n} % Indexed dual variable
    expression Z(n,n)
    %Z = zeros(n,n);
    for j = 1:N
        Z = Z + w(j)*(MeasurementMatrix(j,:)'*MeasurementMatrix(j,:));
    end        
    minimize(sum(w) )
    % Constraints
    subject to    
        for i = 1:n
            e_i = zeros(n,1);
            e_i(i) = 1;
            y{i} : [Z e_i; e_i.' u(i)] >= 0;
        end    
        sum(u) <= gamma;
        0 <= u;

        0 <= w;
        w <= 1;         
cvx_end
