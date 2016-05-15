close all
clear all
clc
%% ************************  Introduction ********************************
% This is a main function to estimate the mean computation time of five
% sensor placement algorithms: the convex relaxation method (Joshi & Boyd
% 2009), SparSenSe(Jamali-Rad et.al. 2014), FrameSense (Ranieri et.at.
% 2014), Minimum Nonzero Eigenvalue Puisuit and Maximal Projection on
% Minimum Eigenspace. 
%
%The result of this main function is the Fig. 7 in the following paper:
%
% C. Jiang, Y. C. Soh and H. Li, Sensor placement by maximal projection on
% minimum eigenspace for linear inverse problems, 2015. (submitted to IEEE-
% TSP).
%
%SparSenSe is very time consuming if the Rowsize (the number of posential
%sensing locations) is large. We set the repeatTimes = 50 and this code is
%run in laptop with 2.4GHz Inter i3-3110M processor. The computation time
%is around 23 hours.
%
%ATTENTION: If the repeatTime is small, the estimation will be not very
%           accurate.
%
% AUTHOR Information:
% Jiang Chaoyang, EEE, NTU 
% Email: cjiang003@e.ntu.edu.sg, chaoyangjiang@hotmail.com
% Finished at 15-OCT-2015
%% ************************************************************************
%% *************************** Initialization *****************************
n = 20;
sensorNumber = 20; % M
gamma = 1.5; % for SparSenSe, The MSE threshold
MaxRowSize = 1000; % the maximum N
caseNumber = MaxRowSize/100;
T = zeros(5,caseNumber);

repeatTimes =50;
% ****************** Estimation of the computaion time ********************
for  i= 1:repeatTimes
    
    for RowSize = 100:100:MaxRowSize % N
        Sequence = RowSize/100;
        V = randn(RowSize,n); % the signal representation matrix \tilde{\Phi}

        %****************** convex relaxation ****************
        tic; 
        defaultParameter1 = ConvexOpt(V',sensorNumber);
        T(1,Sequence)= T(1,Sequence)+ toc;

        %********************* SparSenSe *********************  
        tic;
        defaultParameter2 = SparSenSe_CVX(V,gamma);
        T(2,Sequence)= T(2,Sequence)+ toc;

        %**** FrameSense *******************************************
        tic;
        defaultParameter3 = FrameSense(V,sensorNumber);
        T(3,Sequence)= T(3,Sequence)+toc;

        %**** Minimum nonzero eigenvalue pursuit *****************  
        tic;
        defaultParameter4 = MNEP(V,sensorNumber);
        T(4,Sequence)= T(4,Sequence)+toc;

        %**** Maximal projection on minimum eigenspace *******
        tic;
        defaultParameter5 = MPME(V,sensorNumber);
        T(5,Sequence)= T(5,Sequence)+toc;
    end
end
%********************** Plot the result *****************************
T=T/repeatTimes;
RowSizeSequence = 100:100:MaxRowSize;

plot(RowSizeSequence,T(1,:)','ko',...    
     RowSizeSequence,T(2,:)','m*',... 
     RowSizeSequence,T(3,:)','gs',... 
     RowSizeSequence,T(4,:)','b+',...
     RowSizeSequence,T(5,:)','rx','linewidth',2);
  xlabel('number of potential sensor locations');  ylabel('mean computation time(second)');
  title('100\times20 Gaussian random matrices');
  legend('convex relaxation',...
       'SparSenSe',...
       'FrameSense',...
       'MNEP',...
       'MPME')
