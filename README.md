# Sensor_Placement


The Matlab Code is used in the paper: 
    
    Sensor placement by maximal projection on minimum eigenspace for linear inverse problems, IEEE-TSP, 2016


All the code have been tested using MATLAB 7.8.0 (2009a)



Before using the code, please make sure that you have installed the SDPT3‐4.0 solver and CVX toolbox.


1. The CovnexOpt.m requires the SDPT3‐4.0.
    
        The users can freely download the SDPT3‐4.0 solver from: www.math.nus.edu.sg/~mattohkc//sdpt3.html
        
        Every time before using the code, do the following two steps to start up the SDPT3‐4.0 solver:
        
                Step one: run Installmex
            
                Step two: run startup (run both in command window)
        
        For more detail about the SDPT3‐4.0, please visit: www.math.nus.edu.sg/~mattohkc//sdpt3.html


2. The SparSenSe_CVX.m requires the CVX toolbox.

        To acquire the CVX toolbox, please Visit: http://cvxr.com/cvx/



Package List:
  
  Three main functions:
      
      MainAccuracyComparision.m is used to generate Fig. 1, Fig. 2, Fig. 5, and Fig. 6.
      
      MainSparSenSe1500.m is used to generate Fig. 3 and Fig. 4.
      
      MainTimeCostComparison.m is used to generate Fig. 7.
  
  Five Sensor placement algorithms:
      
      Convexopt.m
      
      SparSenSe_CVX.m
      
      FrameSense.m
      
      MNEP.m
      
      MPME. m
      
  Other functions:
    
        LocalOptimization.m: local optimization in terms of MSE index
        
        LocalOptimizationWCEV.m: local optimization in terms of WCEV index
        
        Criteriacomputation.m: compute the MSE index, the WCEV index and the condition number of the dual observation matrix for a given sensor configuration.


Author Information: 

                    Chaoyang Jiang,
                    
                    School of Electrical and Electronic Engineering,
                    
                    Nanyang Technological University, Singapore
                    
                    E‐mail: chaoyangjiang@hotmail.com, cjiang003@e.ntu.edu.sg

If you meet any problem in using the code, please feel free to contact with the author
