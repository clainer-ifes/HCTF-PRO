clear;
clc;
%format bank;

% Internal settings of the probabilistic analysis
ProbabilityCurvePoints = 500;
NumberOfInputVariableClasses = 5;
NumberOfOutputVariableClasses = 15;
MultiplicativeFactorForExpandedResultCurvePoints = 1e4;

% Prompt the user to specify the desired PDF format.
Chosen_PDF_Default = 1;
Chosen_PDF = input('Select the desired PDF [default: pdf #1]:\n1 - Lognormal\n2 - Weibull\n3 - Gamma\n4 - Normal\n\n');
if(length(Chosen_PDF)==0)
    Chosen_PDF = Chosen_PDF_Default;
end

% Prompt the user to specify the desired output quantity.
Output_Quantity_Default = 1;
Output_Quantity = input('Select the desired output quantity [default: output quantity #1]:\n1 - Turbidity Removal Efficiency\n\n');
if(length(Output_Quantity)==0)
    Output_Quantity = Output_Quantity_Default;
end

% Prompt the user to specify the confidence interval - Considering both sides of the probability curve combined.
ConfidenceInterval_Default = 90;
ConfidenceInterval = input('Select the confidence interval [default: 90%]: ');
if(length(ConfidenceInterval)==0)
    ConfidenceInterval = ConfidenceInterval_Default;
end

% Prompt the user to specify the Relative Standard Deviation quantity - RSD
RelativeStandardDeviation_Default = 10;
RelativeStandardDeviation_Chosen = input('Enter the relative standard deviation - RSD [default: 10%]: ');
if(length(RelativeStandardDeviation_Chosen)==0)
    RelativeStandardDeviation_Chosen = RelativeStandardDeviation_Default;
end

% Model Coefficient Vector
ModelCoefficients = [0.8731 1.8151e-6 1.5153e-5 28.63773 0.0035];

% Load the HTF data
test_data;

% Ask the user to identify the HCTF they wish to analyze.
ID_HCTF_Chosen_Default = 1;
ID_HCTF_Chosen = input(['Enter the ID of the HCTF you want to analyze [between 1 and ' num2str(size(Model_Data,1),'%02d') ' - default: HCTF #1]: ']);
if(length(ID_HCTF_Chosen)==0)
    ID_HCTF_Chosen = ID_HCTF_Chosen_Default;
end
HCTF_Data_Chosen = HCTFs_Data(ID_HCTF_Chosen, :);
Model_Data_Chosen = Model_Data(ID_HCTF_Chosen:(ID_HCTF_Chosen+7),:);

% Define the probabilistic variables.
MeanValue_G = HCTF_Data_Chosen(2);
StandardDeviation_G = RelativeStandardDeviation_Chosen / 100 * MeanValue_G;

MeanValue_Q = HCTF_Data_Chosen(8);
StandardDeviation_Q = RelativeStandardDeviation_Chosen / 100 * MeanValue_Q;

% Define the constants.
d = HCTF_Data_Chosen(3);
p = HCTF_Data_Chosen(4);
D = HCTF_Data_Chosen(5);

ro = 997; %Density 
mi = 0.00089; % Viscosity

c10_deterministic = (ModelCoefficients(2)*MeanValue_G*pi*(d^2)/(4*MeanValue_Q));

L_optimal = sqrt(ModelCoefficients(4)*p/c10_deterministic);

disp('Optimal length [m]:');
disp(L_optimal);

% **************************************************************

if(Chosen_PDF==1)
	MU = log(MeanValue_G^2 / sqrt(StandardDeviation_G^2+MeanValue_G^2));
    SIGMA = sqrt(log(StandardDeviation_G^2/MeanValue_G^2 + 1));
    ValueVector_G = lognrnd(MU,SIGMA,ProbabilityCurvePoints,1);
    
	MU = log(MeanValue_Q^2 / sqrt(StandardDeviation_Q^2+MeanValue_Q^2));
    SIGMA = sqrt(log(StandardDeviation_Q^2/MeanValue_Q^2 + 1));
    ValueVector_Q = lognrnd(MU,SIGMA,ProbabilityCurvePoints,1);
elseif(Chosen_PDF==2)    
    k = (StandardDeviation_G/MeanValue_G)^(-1.086);   % Intrinsic shape factor of the curve
    c = MeanValue_G / gamma(1+1/k);               % Intrinsic scale factor of the curve
    ValueVector_G = wblrnd(c, k, ProbabilityCurvePoints, 1);
    
    k = (StandardDeviation_Q/MeanValue_Q)^(-1.086);   % Intrinsic shape factor of the curve
    c = MeanValue_Q / gamma(1+1/k);               % Intrinsic scale factor of the curve
    ValueVector_Q = wblrnd(c, k, ProbabilityCurvePoints, 1);
elseif(Chosen_PDF==3)    
    teta = (StandardDeviation_G^2) / MeanValue_G;
    k = MeanValue_G/ teta;
    ValueVector_G = gamrnd(k, teta, ProbabilityCurvePoints, 1);
    
    teta = (StandardDeviation_Q^2) / MeanValue_Q;
    k = MeanValue_Q/ teta;
    ValueVector_Q = gamrnd(k, teta, ProbabilityCurvePoints, 1);
elseif(Chosen_PDF==4)    
    ValueVector_G = normrnd(MeanValue_G, StandardDeviation_G, ProbabilityCurvePoints, 1);
    
    ValueVector_Q = normrnd(MeanValue_Q, StandardDeviation_Q, ProbabilityCurvePoints, 1);
end

[P_ValueVector_G, ValueVector_G] = hist(ValueVector_G, NumberOfInputVariableClasses);
P_ValueVector_G = P_ValueVector_G / sum(P_ValueVector_G);

[P_ValueVector_Q, ValueVector_Q] = hist(ValueVector_Q, NumberOfInputVariableClasses);
P_ValueVector_Q = P_ValueVector_Q / sum(P_ValueVector_Q);

Results = [];

for index_G=1:size(ValueVector_G,2)
    for index_Q=1:size(ValueVector_Q,2)
        % **************************************************************
        % **************************************************************
        % Here, the calculation function for each combination of values is established.

        if(Output_Quantity==1)
            Part1 = ModelCoefficients(1);
            Part2 = ((-ModelCoefficients(2)*ValueVector_G(index_G)*L_optimal*pi*(d^2)) / (4*ValueVector_Q(index_Q)));
            Part3 = ((-ModelCoefficients(3)*4*ro*ValueVector_Q(index_Q)) / (mi*pi*d));
            Part4 = ((-ModelCoefficients(4)*p) / (L_optimal));
            Part5 = ((+ModelCoefficients(5)*D) / (d));

            ResultFunction = (Part1 + Part2 + Part3 + Part4 + Part5)*100;
        end

        % **************************************************************
        % **************************************************************

        JointProbability = P_ValueVector_G(index_G)*P_ValueVector_Q(index_Q);

        Results = [Results; ResultFunction JointProbability;];
    end
end

% Capture the results WITHOUT repetitions
Concatenated_result = unique(Results(:,1));

% Add a column of zeros to the results without repetitions
zero_vector=zeros(size(Concatenated_result, 1),1);
Concatenated_result = [Concatenated_result zero_vector];

% Determine the probability for each item in the results without repetitions.
for index_final=1:size(Concatenated_result,1)
    Mask = Results(:,1)==Concatenated_result(index_final);
    
    Probability = sum(Results(:,2).*Mask);
    
    Concatenated_result(index_final,2) = Probability;
end

% Expand the results with thousands of points according to the calculated probability.
Expanded_result = [];
Concatenated_result = [Concatenated_result Concatenated_result(:,2)*MultiplicativeFactorForExpandedResultCurvePoints];
for index_1=1:size(Concatenated_result,1)
    index_2=1;
    while(index_2<=Concatenated_result(index_1,3))
        Expanded_result=[Expanded_result; Concatenated_result(index_1,1)];
        
        index_2=index_2+1;
    end
end

% Generate the result graphs
[Y,X] = hist(Expanded_result, NumberOfOutputVariableClasses);
Y=Y/sum(Y)*100;
Y_accumulated = cumulative_distribution(Y);
figure('Color',[1 1 1]);
hold on;
bar(X, Y, 1, 'FaceColor', [0.8 0.8 0.8]);

% Add labels
if(Output_Quantity==1)
    xlabel('Turbidity Removal Efficiency [%]', 'FontSize',14);
end
ylabel('Probability of Occurrence [%]', 'FontSize',14);
set(gca,'FontSize', 12);

% Set the graph limits
UpperLimitX = (floor(max(X))+1);
LowerLimitX = (floor(min(X))-0.5);

if(LowerLimitX==UpperLimitX)
    UpperLimitX=UpperLimitX+1;
    LowerLimitX=LowerLimitX-1;
end
UpperLimitY = (floor(max(Y))+1);
grid on;

axis([LowerLimitX, UpperLimitX, 0, UpperLimitY]);

% Generate the highlighted area based on the selected quantity
if(Output_Quantity==1)
    Part1 = ModelCoefficients(1);
    Part2 = ((-ModelCoefficients(2)*MeanValue_G*L_optimal*pi*(d^2)) / (4*MeanValue_Q));
    Part3 = ((-ModelCoefficients(3)*4*ro*MeanValue_Q) / (mi*pi*d));
    Part4 = ((-ModelCoefficients(4)*p) / (L_optimal));
    Part5 = ((+ModelCoefficients(5)*D) / (d));

    DETERMINISTICValueOfSelectedQuantity = (Part1 + Part2 + Part3 + Part4 + Part5)*100;
end

Vector = DETERMINISTICValueOfSelectedQuantity<X;

index_1=1;
while(index_1<=length(Vector))
   if(Vector(index_1)==1)
       IDClass = index_1;
       index_1 = length(Vector) + 1;
   end

   index_1 = index_1 + 1;
end

VisualLimitWithPreviousClass = ((X(IDClass) - X(IDClass-1))/2+X(IDClass-1));
if(DETERMINISTICValueOfSelectedQuantity<VisualLimitWithPreviousClass)
    VisualClassID = IDClass - 1;
else
    VisualClassID = IDClass;
end

if(Output_Quantity==1)
    bar(X(1,(VisualClassID+1):length(X)), Y(1,(VisualClassID+1):length(Y)), 1, 'b');
end

if(Output_Quantity==1)
    PositionX = DETERMINISTICValueOfSelectedQuantity;
    PositionY = 0;
    Width = (((X(VisualClassID+1) - X(VisualClassID))/2+X(VisualClassID)) - DETERMINISTICValueOfSelectedQuantity);
    Height = Y(VisualClassID);
end

rectangle('Position',[PositionX PositionY Width Height], 'FaceColor', 'b');

% Plot the lower and upper limits of the confidence interval and inform the user of the values.
LowerLimitValueConfidenceInterval=prctile(Expanded_result, (100-ConfidenceInterval)/2);
plot([LowerLimitValueConfidenceInterval LowerLimitValueConfidenceInterval], [0 UpperLimitY], 'Color', 'r', 'LineWidth', 2);

disp('LOWER limit:');
disp(LowerLimitValueConfidenceInterval);

UpperLimitValueConfidenceInterval=prctile(Expanded_result, 100-(100-ConfidenceInterval)/2);
plot([UpperLimitValueConfidenceInterval UpperLimitValueConfidenceInterval], [0 UpperLimitY], 'Color', 'r', 'LineWidth', 2);

disp('UPPER limit:');
disp(UpperLimitValueConfidenceInterval);

disp('DETERMINISTIC Value:');
disp(DETERMINISTICValueOfSelectedQuantity);

% Draw a line corresponding to the value of the quantity being analyzed
if(Output_Quantity==1)
    plot([DETERMINISTICValueOfSelectedQuantity DETERMINISTICValueOfSelectedQuantity], [0 UpperLimitY], 'y', 'LineWidth', 3);
end

% Calculate the highlighted area in the figure and inform the user.
AreaCalculationVector = DETERMINISTICValueOfSelectedQuantity>X;
LowerArea = sum(AreaCalculationVector.*Y);

DeltaY = Y(IDClass);
DeltaX = X(IDClass) - X(IDClass-1);
XVariation = DETERMINISTICValueOfSelectedQuantity - X(IDClass-1);
LowerArea = LowerArea + XVariation/DeltaX*DeltaY;

if(Output_Quantity==1)
    QuantityArea = 100 - LowerArea;
end

disp('Highlighted area:');
disp(QuantityArea);


% ***********************************
% ***********************************
% Comparative graph of the optimal L with the model
% ***********************************
% ***********************************

figure('Color',[1 1 1]);
hold on;

% Add labels
xlabel('Length [m]', 'FontSize',14);
ylabel('Turbidity Removal Efficiency [%]', 'FontSize',14);
set(gca,'FontSize', 12);
grid on;

% Set the graph limits
UpperLimitX = (floor(max(Model_Data_Chosen(:,2)))+1);
LowerLimitX = (floor(min(Model_Data_Chosen(:,2)))-1);

if(LowerLimitX==UpperLimitX)
    UpperLimitX=UpperLimitX+1;
    LowerLimitX=LowerLimitX-1;
end
UpperLimitY = (floor(max(Model_Data_Chosen(:,3)))+1);
LowerLimitY = (floor(min(Model_Data_Chosen(:,3)))-1);

axis([LowerLimitX, UpperLimitX, LowerLimitY, UpperLimitY]);

[MaximumY, MaximumYIndex] = max(Model_Data_Chosen(:,3));

if(MaximumYIndex==1)
    PositionX = 0;
else
    PositionX = Model_Data_Chosen(MaximumYIndex-1, 2);
end
PositionY = LowerLimitY;
Width = Model_Data_Chosen(MaximumYIndex+1, 2) - PositionX;
Height = (UpperLimitY-LowerLimitY);

rectangle('Position',[PositionX PositionY Width Height], 'FaceColor',[0.9 0.9 0.9], 'LineStyle', 'none');

scatter(Model_Data_Chosen(:,2), Model_Data_Chosen(:,3), 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');

scatter(L_optimal, DETERMINISTICValueOfSelectedQuantity, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
legend('Turbidity removal efficiency from Oliveira et al. (2019)', 'Expected turbidity removal efficiency at the optimal length', 'location', 'southwest');

disp('LOWER limit of the length:');
disp(PositionX);

disp('UPPER limit of the length:');
disp(PositionX+Width);


