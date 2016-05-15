close all
clear all
clc
%% ************************  Introduction ********************************
%     This is a main function to estimate the selected sensor indices used SparSenSe
%     and mean MSE index of the solutions of five sensor placement algorithms: 
%     the convex relaxation method (Joshi & Boyd 2009), SparSenSe(Jamali-Rad 
%     et.al. 2014), FrameSense (Ranieri et.at.2014), Minimum Nonzero Eigenvalue 
%     Puisuit and Maximal Projection on Minimum Eigenspace. 
% 
%     The results of this main function is used as Fig. 3, and Fig. 4 
%     in the following paper:
%
%     C. Jiang, Y. C. Soh and H. Li, Sensor placement by maximal projection on
%     minimum eigenspace for linear inverse problems, 2015. (submitted to IEEE-
%     TSP).
%
%
%
%     AUTHOR Information:
%     Jiang Chaoyang, EEE, NTU, Singapore
%     Email: cjiang003@e.ntu.edu.sg, chaoyangjiang@hotmail.com
%
%     Finished at 15-OCT-2015
%% ************************************************************************

%% *************************** Initialization *****************************
rowNO =1500;
n = 20; 
MaxsensorNO = 40;
%*********** generate 1500x20 Gaussian random matrices ***********
gamma = 1.5; % for SparSenSe, The MSE threshold
WCEV_Threshold =0.3; % for SparSenSe
V = randn(rowNO,n);

CriteriaMatrixConvex = zeros(3,MaxsensorNO-n+1);
CriteriaMatrixFrameSense = zeros(3,MaxsensorNO-n+1);
CriteriaMatrix_MNEP = zeros(3,MaxsensorNO-n+1);
CriteriaMatrix_MPME = zeros(3,MaxsensorNO-n+1);
CriteriaMatrix_Sparse = zeros(3,MaxsensorNO-n+1);

%% ********** Five methods:the results **************************

%********************************* SparSenSe ******************************
w = SparSenSe_CVX(V,gamma);
[wRank,SensorPositionSparSenSe]=sort(w,'descend');
%*********** compute the criteria: MSE, WCEV & Cond***********
for sensorNumber = n:1:MaxsensorNO
[MSE, WCEV, ConditionalNO] = Criteriacomputation(SensorPositionSparSenSe(1:sensorNumber), V);
CriteriaMatrix_Sparse(:,sensorNumber-n+1)=[MSE;WCEV;ConditionalNO];         
end
    %***************** SparSenSe: show the sensing indices ****************
    selThresh=0.02;
    Storagew = w;
    wCentr=w;
    wCentr(wCentr < selThresh) = 0;
    indCentr = find(wCentr);
    SensorNO_SparSenSe = length(indCentr);
    valCentr = wCentr(indCentr);     
    figure(1)
    h1 = stem(indCentr,valCentr,'--','markersize',14,'color','k','linewidth',2);     
    axis([1 rowNO 0 max(wCentr)+0.2])  % Axis (Play with 0.2 if necessary!) 
    hold on;
    SensorIndex = 1:rowNO;
    OneMatrixSensorIndex = ones(size(SensorIndex));
    plot(SensorIndex,0.05*OneMatrixSensorIndex,'g',...
         SensorIndex,0.1*OneMatrixSensorIndex,'r'); 
    xlabel('sensor index')
    ylabel('w^*')
    title('1500\times20 Gaussian random matrices')
    legend('selected sensor indices when \tau=0.02',...
         '\tau=0.05',...
         '\tau=0.1')
%**************************************************************************


%************** convex relaxation, FrameSense, MNEP & MPME ****************
for sensorNumber = n:1:MaxsensorNO   
    %****************** convex relaxation ****************
    w = ConvexOpt(V',sensorNumber);
    [wRank,SensorPositionCOV]=sort(w,'descend');        
    [MSE, WCEV, ConditionalNO] = Criteriacomputation(SensorPositionCOV(1:sensorNumber), V);
    CriteriaMatrixConvex(:,sensorNumber-n+1)=[MSE;WCEV;ConditionalNO];  

    %**** FrameSense *******************************************
    SensorPosition = FrameSense(V,sensorNumber);
    [MSE, WCEV, ConditionalNO] = Criteriacomputation(SensorPosition, V);
    CriteriaMatrixFrameSense(:,sensorNumber-n+1)=[MSE;WCEV;ConditionalNO];

    %**** Minimum nonzero eigenvalue pursuit *****************
    SensorPosition = MNEP(V,sensorNumber);
   [MSE, WCEV, ConditionalNO] = Criteriacomputation(SensorPosition, V);
    CriteriaMatrix_MNEP(:,sensorNumber-n+1)=[MSE;WCEV;ConditionalNO];

    %**** Maximal projection on minimum eigenspace *******
    SensorPosition = MPME(V,sensorNumber);
    [MSE, WCEV, ConditionalNO] = Criteriacomputation(SensorPosition, V);
    CriteriaMatrix_MPME(:,sensorNumber-n+1)=[MSE;WCEV;ConditionalNO]; 
end
%********************** Plot the result *****************************
SensorNumber = n:MaxsensorNO;
OneMatrix = ones(size(SensorNumber));


figure(2)
plot(SensorNumber,CriteriaMatrixConvex(1,:)','-ko',...
     SensorNumber,CriteriaMatrix_Sparse(1,:)','-m*',...
     SensorNumber,CriteriaMatrixFrameSense(1,:)','-gs',...
     SensorNumber,CriteriaMatrix_MNEP(1,:)','-b+',...
     SensorNumber,CriteriaMatrix_MPME(1,:)','-rx',...
     SensorNumber,gamma*OneMatrix,'g'); 
  xlabel('number of sensor nodes');  ylabel('MSE index');
  title('1500\times20 Gaussian random matrices')
  legend('convex relaxation',...
       'SparSenSe',...
       'FrameSense',...
       'MNEP',...
       'MPME',...
       ['MSE index threshold=',num2str(gamma)])
%**************************************************************************