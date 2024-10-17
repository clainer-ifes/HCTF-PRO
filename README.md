HCTF-PRO: A Code for Probabilistic Analysis of Process Efficiency in Helically Coiled Tube Flocculators

Authors:
- Danieli Soares de Oliveira (Federal Institute of Espírito Santo - campus Cariacica)
- Clainer Bravin Donadel (Federal Institute of Espírito Santo - campus Vitória)

Contact Information:
- danieli@ifes.edu.br
- cdonadel@ifes.edu.br

Introduction:
HCTF-PRO is a computational tool designed for the probabilistic analysis of helically coiled tube flocculators (HCTFs) used in water treatment. The tool helps determine the optimal length of the HCTF and conducts a probabilistic performance analysis considering uncertainties in hydraulic parameters.

Code Metadata:
- Current code version: v1
- Permanent link to code repository: https://github.com/clainer-ifes/HCTF-PRO
- Legal code license: MIT License
- Software code language: MATLAB®

Software Usage Instructions:

1. Download or clone the repository from GitHub:
   git clone https://github.com/clainer-ifes/HCTF-PRO

2. Ensure MATLAB® is installed on your system.

3. Open MATLAB® and set the working directory to the location of the cloned repository.

4. Input Data:
   - Probability density function: Select from Lognormal, Weibull, Gamma, or Normal.
   - Output variable: Turbidity removal efficiency (default variable analyzed).
   - Confidence interval: Default is 90%.
   - Relative standard deviation (RSD): Default is 10%.
   - HCTF Group: Indicate the ID of the first HCTF to analyze from the provided dataset.

5. Running the Code:
   Execute the MATLAB® script to begin the analysis. The software will prompt the user for input parameters and will perform the probabilistic analysis based on the dataset 'test_data.m'.

6. Output:
   - Optimal flocculator length.
   - Probabilistic analysis results, including graphical histograms and text output of confidence intervals and turbidity removal efficiency.

License:
This software is provided under the MIT License.
