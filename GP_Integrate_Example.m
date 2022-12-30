%BSS-ANOVA GP Integration Example Script. Kyle Hayes, using code from Dave Mebane.
%This case assumes you are starting your BSS-ANOVA model building from scratch, 
%and have fairly minimal knowledge of the code found in the 
%ESMS-Group-Public/FoKL-GP repository. If your models are already built, 
%or you are comfortable builing these models, the integration of models
%begins on line 94.

%Our first step is to load and format the data correctly in order to create 
%a BSS-ANOVA model using emulator.
%For the purposes of this example we will use data concerning a set of
%cascaded tanks (2 outputs, 3 inputs)

dataset1 = load('Tank2.mat');%Replace with your data

%In this case, the training output data for our models comes from differentiating 
%the two columns of the 'y' portion of the dataset1 structure. The 'u'
%structure contains one of the inputs to the model (the forcing function).

y = dataset1.y';
u = dataset1.u';

%Now that the data has been loaded, we will create matrices of training
%inputs and outputs. For the purposes of this example we will use half the
%data for creating the model and half for testing
%Differentiation is carried out here by a center difference method since
%the data is smooth. If the data is noisy, using a smoothing method or
%alternate differentiation method. If derivatives are already available,
%skip this step.

%In order to differentiate, we note that the data is sampled at a uniform 4
%second rate
traindata1 = gradient(y(1:(end),1),4);
traindata2 = gradient(y(1:(end),2),4);

%Now we create matricies of the training inputs. All inputs need to be
%normalized on a 0 to 1 range. In this case both models use the same sets
%of inputs, so only one input matrix will be created, but if the inputs
%differed for each model seperate matrices would need to be built for each.

%The order of the inputs in the input matrices matters for the integrator.
%Although inconsistantly ordered matrices can be corrected later, it will be
%fastest if input matrices are all set up in the following manner:
%The first columns should be the undifferentiated verisons of the output
%itself (the data for the 1st model, then the 2nd and so on)
%The following columns should contain any other inputs, which can be
%ordered in any manner, so long as that order is consistant for all models.

traininputs = normalize([y(1:(end),:),u(1:end)],'range');


%Creating BSS-ANOVA models

%The BSS-ANOVA model requires the selection hyperparameters (a, b, atau,
%and btau). 4 is a good starting value for a and atau. Increase the value of
%a to decrease the 'spread' of the predictions made by the model. b and btau can be
%calculated from the chosen values of a and atau, as well as an initial
%guess concerning the varience of the noise in the datasets (sigma squared), 
%and for btau the scale of the data

a=1000;
atau=4;
sigmasq = 10^(-3); %Due to the smoothness of the data
scale1 = abs(mean(traindata1));
scale2 = abs(mean(traindata2));
b = sigmasq*(a+1);
btau1 = (scale1/sigmasq)*(atau+1);
btau2 = (scale2/sigmasq)*(atau+1);

%Now we load the splines containing our basis set
x = dlmread('spline_coefficient_500.txt');
phis = splineconvert500(x);

%We can now create our BSS-ANOVA models. Default values are used for
%several inputs to emulator, better results may be obatained with fine
%tuning these variables
[betas1, mtx1, evs1] = emulator(traininputs, traindata1, phis, [1;1;1;1;1;1], a, b, atau, btau1, 3, 2000, 0, 1, 0, 0, 100, 0);
[betas2, mtx2, evs2] = emulator(traininputs, traindata2, phis, [1;1;1;1;1;1], a, b, atau, btau2, 3, 2000, 0, 1, 0, 0, 100, 0);

%It's a good idea to throw out the draws from the 'burn in' period; the 1st
%1000 draws should be a safe estimate

betas1 = betas1(1000:end,:);
betas2 = betas2(1000:end,:);

%Models are now created for each of the derivatives we seeked to model. We
%can plot the predictions of the derivatives using coverage, or evaluate
%these predictions for test data using bss_eval. An example of using coverage
%for the first model is shown below:
coverage(betas1, traininputs, traindata1, phis, mtx1, 50, 1)


% For this example, we will now move into integrating these models.

%Integration

%Integration requires specifying several basics: initial conditions, step
%sizes, start and stop points, and the values of any additional inputs.
start = 4;
stop = 3750*4;
stepsize = 4;
%Inputs to the integrator should be placed in the same order as they were
%in the integrator. Since only one set of inputs is provided to the
%integrator (since all models need to be evaluated iteratively) if models
%were built with their inputs in different orders, then the used_inputs
%input to GP_Integrate can be provided additional information to rectify
%the issue (see the documentation of GP_Integrate for more information).

%In the case where all inputs are used for all models (in the correct order),
%used_inputs should be a cell matrix constisting of vectors of ones
used_inputs = {[1,1,1],[1,1,1]};

%Our initial condition will be the 1st value of the 2nd half of the
%dataset:
ic = y(end/2,:)';

%'u' is our forcing function: an independent variable that varies with
%time. The integrator can handle an arbitrary number of these inputs so
%long as their value is known (or can be approximated) at each time step
%within the integrated range. Here u is our only additional input, but
%for any cases with additional inputs, GP_Integrate would be provided a
%matrix of all these inputs rather than the vector shown here.
%All additional inputs still need to be normalized
utest = normalize(u(end/2+1:end),'range');

%All derivative models need to be concatenated into a cell array (in the
%same order as they are contained within the inputs). If the mean
%prediction is all that is required, then the mean beta model should be
%provided to GP_Integrate (accomplished by averaging the columns of the
%beta matricies)
betas = {mean(betas1);mean(betas2)};
matricies = {mtx1;mtx2};

%norms is the matrix of min and max values of the outputs in the test set
%(used for denormalizing the predictions)
norms = [min(y(1:end/2,:));max(y(1:end/2,:))];

%We now run GP_Integrate
[T,Y] = GP_Integrate(betas,matricies,utest,norms,phis,start,stop,ic,stepsize,used_inputs);


%Simple plots comparing the predictions to the actual data
figure
hold on
plot(T,Y(1,:),'b')
hold on
plot(T,y(end/2+1:end,1),'r')


figure
hold on
plot(T,Y(2,:),'b')
hold on
plot(T,y(end/2+1:end,2),'r')

