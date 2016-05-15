close all
clear all
clc
%% ************************  Introduction ********************************
%     This is a main function to estimate the mean MSE index, mean WCEV index,
%     and mean condition number of the solutions and the locally optimizated 
%     solutions of five sensor placement algorithms: the convex relaxation 
%     method (Joshi & Boyd 2009), the SparSenSe method (Jamali-Rad et.al. 2014), 
%     the FrameSense method (Ranieri et.at. 2014), Minimum Nonzero Eigenvalue 
%     Puisuit and Maximal Projection on Minimum Eigenspace. 
% 
%     The results of this main function is used as Fig. 1, Fig. 2, Fig.5
%     and Fig. 6 in the following paper:
%
%     C. Jiang, Y. C. Soh and H. Li, Sensor placement by maximal projection on
%     minimum eigenspace for linear inverse problems, 2015. (submitted to IEEE-
%     TSP).
%
% 
%      ATTENTION: For different examples, please properly choose the database
%
%
%     AUTHOR Information:
%     Jiang Chaoyang, EEE, NTU, Singapore
%     Email: cjiang003@e.ntu.edu.sg, chaoyangjiang@hotmail.com
%
%     Finished at 10-sep-2014 
%     Latest version at 15-OCT-2015
%% ************************************************************************

%% *************************** Initialization *****************************

% ****************** Initialize some parameters ***************************
MonteCarlNO = 200;
rowNO =100;
n = 20; 
MaxsensorNO = 40;
DatabaseRandom = cell(MonteCarlNO,1);

%********************** Generate the database ****************************
% for each simulation, choose one of the following three types of database.

%*********** generate 100x20 Gaussian random matrices ***********
gamma = 1.5; % for SparSenSe, The MSE threshold
WCEV_Threshold =0.3; % for SparSenSe
for i = 1:MonteCarlNO
    V = randn(rowNO,n);
    DatabaseRandom{i,1}= V;    
end

%********** generate 100x20 Bernoulli random matrices ***********
% gamma = 8; % for SparSenSe, The MSE threshold
% WCEV_Threshold =1.5; % for SparSenSe
% for i = 1:MonteCarlNO   
%     V = binornd(1,0.5,[rowNO,n]);   
%     DatabaseRandom{i,1}= V;    
% end

%********** generate 100x20 normalized Gaussian random matrices ***********
% gamma = 35; % for SparSenSe, The MSE threshold
% WCEV_Threshold =5; % for SparSenSe
% for i = 1:MonteCarlNO
%     V = randn(rowNO,n);
%     for j = 1:rowNO
%         V(j,:)=V(j,:)/sqrt(V(j,:)*V(j,:)');
%     end
%     DatabaseRandom{i,1}= V;    
% end

%***************** Initilize some used matrices and vectors **********
CriteriaMatrixConvex = zeros(3,MaxsensorNO-n+1);
CriteriaMatrixFrameSense = zeros(3,MaxsensorNO-n+1);
CriteriaMatrix_MNEP = zeros(3,MaxsensorNO-n+1);
CriteriaMatrix_MPME = zeros(3,MaxsensorNO-n+1);
CriteriaMatrix_Sparse = zeros(3,MaxsensorNO-n+1);

%***** used for local optimization
SensorPositionMatrixCov = cell(MonteCarlNO,MaxsensorNO-n+1);
SensorPositionMatrixSpar = zeros(MonteCarlNO,MaxsensorNO);
SensorPositionMatrixMNEP = zeros(MonteCarlNO,MaxsensorNO);
SensorPositionMatrixMPME = zeros(MonteCarlNO,MaxsensorNO);
SensorPositionMatrixFrame = cell(MonteCarlNO,MaxsensorNO-n+1);

%% ********** Five methods:the results **************************

for MonteCarlNumber = 1:MonteCarlNO     
    V = DatabaseRandom{MonteCarlNumber,1};    
    %************ results of convex relaxation & FrameSense *************
    for sensorNumber = n:1:MaxsensorNO 
        %****************** convex relaxation ****************
        w = ConvexOpt(V',sensorNumber);
        [default_wRank,SensorPositionCOV]=sort(w,'descend');        
        [MSE, WCEV, ConditionalNO] = Criteriacomputation(SensorPositionCOV(1:sensorNumber), V);
        CriteriaMatrixConvex(:,sensorNumber-n+1)=CriteriaMatrixConvex(:,sensorNumber-n+1)+[MSE;WCEV;ConditionalNO];  

        SensorPositionMatrixCov{MonteCarlNumber,sensorNumber-n+1}=SensorPositionCOV(1:sensorNumber);% used for local optimization
        %**** FrameSense *******************************************
        SensorPosition = FrameSense(V,sensorNumber);
        [MSE, WCEV, ConditionalNO] = Criteriacomputation(SensorPosition, V);
        CriteriaMatrixFrameSense(:,sensorNumber-n+1)=CriteriaMatrixFrameSense(:,sensorNumber-n+1)+[MSE;WCEV;ConditionalNO];
        
        SensorPositionMatrixFrame{MonteCarlNumber,sensorNumber-n+1}=SensorPosition;% used for local optimization
    end
    % ************* End of  convex relaxation & FrameSense ****************
    
    
    %********************* results of SparSenSe, MNEP & MPME **************  
    %******************* Determine the sensor placement *******************
    %****************** SparSenSe *************************
    w = SparSenSe_CVX(V,gamma);
    [wRank,SensorPositionSparSenSe]=sort(w,'descend');     
    SensorPositionMatrixSpar(MonteCarlNumber,:) = SensorPositionSparSenSe(1:MaxsensorNO); % used for local optimization  
    
    %**** Minimum nonzero eigenvalue pursuit *****************
    SensorPositionMNEP = MNEP(V,MaxsensorNO);    
    SensorPositionMatrixMNEP(MonteCarlNumber,:) = SensorPositionMNEP; % used for local optimization    
    
    %**** Maximal projection on minimum eigenspace *******
    SensorPositionMPME = MPME(V,MaxsensorNO);
    SensorPositionMatrixMPME(MonteCarlNumber,:) = SensorPositionMPME; % used for local optimization    
    %*************End of the determination of sensor placement*************
    
    %************** compute the criteria: MSE, WCEV & Cond ****************    
    for sensorNumber = n:1:MaxsensorNO
        %****************** SparSenSe *************************
        [MSE, WCEV, ConditionalNO] = Criteriacomputation(SensorPositionSparSenSe(1:sensorNumber), V);
        CriteriaMatrix_Sparse(:,sensorNumber-n+1)=CriteriaMatrix_Sparse(:,sensorNumber-n+1)+[MSE;WCEV;ConditionalNO]; 
        
        %**** Minimum nonzero eigenvalue pursuit *****************
        [MSE, WCEV, ConditionalNO] = Criteriacomputation(SensorPositionMNEP(1:sensorNumber), V);
        CriteriaMatrix_MNEP(:,sensorNumber-n+1)=CriteriaMatrix_MNEP(:,sensorNumber-n+1)+[MSE;WCEV;ConditionalNO]; 
        
        %**** Maximal projection on minimum eigenspace *******
        [MSE, WCEV, ConditionalNO] = Criteriacomputation(SensorPositionMPME(1:sensorNumber), V);    
        CriteriaMatrix_MPME(:,sensorNumber-n+1)=CriteriaMatrix_MPME(:,sensorNumber-n+1)+[MSE;WCEV;ConditionalNO];       
    end
    %*****************End of Criteria computation *************************
end %end MonteCarlNumber
%**************** End of SparSenSe, MNEP & MPME ***************************

%***************** find the mean value **********************************
CriteriaMatrixConvex = CriteriaMatrixConvex./MonteCarlNO;
CriteriaMatrixFrameSense = CriteriaMatrixFrameSense./MonteCarlNO;
CriteriaMatrix_MNEP = CriteriaMatrix_MNEP./MonteCarlNO;
CriteriaMatrix_MPME = CriteriaMatrix_MPME./MonteCarlNO;
CriteriaMatrix_Sparse =CriteriaMatrix_Sparse./MonteCarlNO;

%% ********* Plot the comparison of five methods ******************
SensorNumber = n:MaxsensorNO;
OneMatrix = ones(size(SensorNumber));
figure(1)
plot(SensorNumber,CriteriaMatrixConvex(1,:)','-ko',...
     SensorNumber,CriteriaMatrix_Sparse(1,:)','-m*',...
     SensorNumber,CriteriaMatrixFrameSense(1,:)','-gs',...
     SensorNumber,CriteriaMatrix_MNEP(1,:)','-b+',...
     SensorNumber,CriteriaMatrix_MPME(1,:)','-rx',...
     SensorNumber,gamma*OneMatrix,'g'); 
  xlabel('number of sensor nodes');  ylabel('mean MSE index');
  legend('convex relaxation',...
       'SparSenSe',...
       'FrameSense',...
       'MNEP',...
       'MPME',...
       ['MSE index threshold=',num2str(gamma)])
figure(2)
plot(SensorNumber,CriteriaMatrixConvex(2,:)','-ko',...
     SensorNumber,CriteriaMatrix_Sparse(2,:)','-m*',...
     SensorNumber,CriteriaMatrixFrameSense(2,:)','-gs',...
     SensorNumber,CriteriaMatrix_MNEP(2,:)','-b+',...
     SensorNumber,CriteriaMatrix_MPME(2,:)','-rx',...
     SensorNumber,WCEV_Threshold*OneMatrix,'g'); 
  xlabel('number of sensor nodes');  ylabel('mean WCEV index');
   legend('convex relaxation',...
       'SparSenSe',...
       'FrameSense',...
       'MNEP',...
       'MPME',...
       ['WCEV index threshold=',num2str(WCEV_Threshold)])  
figure(3)
plot(SensorNumber,CriteriaMatrixConvex(3,:)','-ko',...
     SensorNumber,CriteriaMatrix_Sparse(3,:)','-m*',...
     SensorNumber,CriteriaMatrixFrameSense(3,:)','-gs',...
     SensorNumber,CriteriaMatrix_MNEP(3,:)','-b+',...
     SensorNumber,CriteriaMatrix_MPME(3,:)','-rx'); 
  xlabel('number of sensor nodes (k) ');  ylabel('mean condition number of \Psi_k');
   legend('convex relaxation',...
       'SparSenSe',...
       'FrameSense',...
       'MNEP',...
       'MPME') 
%**************************************************************************

%% ******** Five methods: Local Optimization ***********************
%**************************************************************************  
MaxiterationNO = 100; % for local optimization
MSEMatrix_loc = zeros(5,MaxsensorNO-n+1);
WCEVMatrix_loc = zeros(5,MaxsensorNO-n+1);
%***** result of different sensor number ********
for MonteCarlNumber = 1:MonteCarlNO 
    V = DatabaseRandom{MonteCarlNumber,1};
    for sensorNumber = n:1:MaxsensorNO
        % ******  convex relaxation *********        
        MSE = LocalOptimization(V',SensorPositionMatrixCov{MonteCarlNumber,sensorNumber-n+1},MaxiterationNO);
        MSEMatrix_loc(1,sensorNumber-n+1) = MSEMatrix_loc(1,sensorNumber-n+1)+MSE;
        
        WCEV= LocalOptimizationWCEV(V',SensorPositionMatrixCov{MonteCarlNumber,sensorNumber-n+1},MaxiterationNO);
        WCEVMatrix_loc(1,sensorNumber-n+1) = WCEVMatrix_loc(1,sensorNumber-n+1)+WCEV;        
        % ******* SparSenSe **********        
        MSE = LocalOptimization(V',SensorPositionMatrixSpar(MonteCarlNumber,1:sensorNumber)',MaxiterationNO);
        MSEMatrix_loc(2,sensorNumber-n+1) = MSEMatrix_loc(2,sensorNumber-n+1)+MSE;

        WCEV= LocalOptimizationWCEV(V',SensorPositionMatrixSpar(MonteCarlNumber,1:sensorNumber)',MaxiterationNO);
        WCEVMatrix_loc(2,sensorNumber-n+1) = WCEVMatrix_loc(2,sensorNumber-n+1)+WCEV;        
        %**** FrameSense ******************************************* 
        MSE = LocalOptimization(V', SensorPositionMatrixFrame{MonteCarlNumber,sensorNumber-n+1}',MaxiterationNO);
        MSEMatrix_loc(3,sensorNumber-n+1) = MSEMatrix_loc(3,sensorNumber-n+1)+MSE;
    
        WCEV = LocalOptimizationWCEV(V',SensorPositionMatrixFrame{MonteCarlNumber,sensorNumber-n+1}',MaxiterationNO);
        WCEVMatrix_loc(3,sensorNumber-n+1) = WCEVMatrix_loc(3,sensorNumber-n+1)+WCEV;       
        %**** Minimum nonzero eigenvalue pursuit *****************
        MSE = LocalOptimization(V',SensorPositionMatrixMNEP(MonteCarlNumber,1:sensorNumber)',MaxiterationNO);
        MSEMatrix_loc(4,sensorNumber-n+1) = MSEMatrix_loc(4,sensorNumber-n+1)+MSE;

        WCEV = LocalOptimizationWCEV(V',SensorPositionMatrixMNEP(MonteCarlNumber,1:sensorNumber)',MaxiterationNO);
        WCEVMatrix_loc(4,sensorNumber-n+1) = WCEVMatrix_loc(4,sensorNumber-n+1)+WCEV;        
        %**** Maximal projection on minimum eigenspace *******        
        MSE = LocalOptimization(V',SensorPositionMatrixMPME(MonteCarlNumber,1:sensorNumber)',MaxiterationNO);
        MSEMatrix_loc(5,sensorNumber-n+1) = MSEMatrix_loc(5,sensorNumber-n+1)+MSE;        

        WCEV = LocalOptimizationWCEV(V',SensorPositionMatrixMPME(MonteCarlNumber,1:sensorNumber)',MaxiterationNO);    
        WCEVMatrix_loc(5,sensorNumber-n+1) = WCEVMatrix_loc(5,sensorNumber-n+1)+WCEV;
    end
end
%%****find the mean value
MSEMatrix_loc = MSEMatrix_loc./MonteCarlNO;
WCEVMatrix_loc = WCEVMatrix_loc./MonteCarlNO;

%% ******** Plot the results: Local optimization *******************
figure(4)
plot(SensorNumber,MSEMatrix_loc(1,:)','-ko',...    
     SensorNumber,MSEMatrix_loc(2,:)','-m*',... 
     SensorNumber,MSEMatrix_loc(3,:)','-gs',... 
     SensorNumber,MSEMatrix_loc(4,:)','-b+',...
     SensorNumber,MSEMatrix_loc(5,:)', '-rx',...
     SensorNumber,CriteriaMatrix_MPME(1,:)','-kd',...
     SensorNumber,gamma*OneMatrix,'g'); 
  xlabel('number of sensor nodes');  ylabel('mean MSE index '); 
  title('local minimization of MSE                100\times20 Gaussian random matrices');
   legend('convex relaxation + local optimization',...
       'SparSenSe + local optimization',...
       'FrameSense + local optimization',...
       'MNEP + local optimization',...
       'MPME + local optimization',...
       'MPME',...
       ['MSE index threshold=',num2str(gamma)]) 
  
 figure(5)
plot(SensorNumber,WCEVMatrix_loc(1,:)','-ko',...    
     SensorNumber,WCEVMatrix_loc(2,:)','-m*',... 
     SensorNumber,WCEVMatrix_loc(3,:)','-gs',... 
     SensorNumber,WCEVMatrix_loc(4,:)','-b+',...
     SensorNumber,WCEVMatrix_loc(5,:)', '-rx',...
     SensorNumber,CriteriaMatrix_MPME(2,:)','-kd',...
     SensorNumber,WCEV_Threshold*OneMatrix,'g'); 
  xlabel('number of sensor nodes');  ylabel('mean WCEV index');  
  title('local minimization of WCEV                100\times20 Gaussian random matrices')
   legend('convex relaxation + local optimization',...
       'SparSenSe + local optimization',...
       'FrameSense + local optimization',...
       'MNEP + local optimization',...
       'MPME + local optimization',...
       'MPME',...
       ['WCEV index threshold=',num2str(WCEV_Threshold)]) 
  
 