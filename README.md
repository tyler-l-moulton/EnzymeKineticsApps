# EnzymeKineticsApps
Apps, readme, and sample data files for sample enzyme kinetics lab simulator 

These apps can be used to explore Michaelis-Menten enzyme kinetics. The "tutorial" app is very straightforward. Users toggle Km and Vmax as well as axis min/max values to see how changing these parameters changes Michaelis-Menten and Lineweaver-Burk plots.
The """ app allows users to run a full experiment to determine the Km of alpha amylase (in the app, it is specified to be 4). This is done using the conversion of DNSA->ANSA (yellow to orange) by reducing sugars as a proxy for enzyme digestion rate. Using this app students can:
1) Create a standard curve to quantify the amount of reducing ends in solution by reacting maltose with DNSA. 
2) Determine optimal pH for enzyme activity of alpha amylase (based on Gangadharan et al 2009 Appl Biochem Biotech)
3) Quantify reducing ends in starch solutions at various concentrations WITHOUT digestion by alpha amylase
4) Quantify reducing ends in starch solutions at various concentrations AFTER digestion by alpha amylase

All concentrations are set by uploading a 96 well plate template csv file (which MUST meet the same dimensions and labels as the template that is downloadable from the app). PH can be specified as constant using a slider, or individually for each well by uploading a separate 96 well plate template csv file. All other parameters are set with dropdown menus.Each time the user clicks "Run Simulation", the amylase digestion is run (if amylase is added to the solution) followed by the reaction with DNSA. An image of the resulting 96 well plate appears, and users can download the corresponding OD readings (at 540nm). They are once again in the format of a 96 well plate. 

The time of the amylase digestion is specified to be 3 minutes. Solutions are reacted with DNSA at 90 C for 10 minues. Enzyme concentration is not specified, but can be specified by the instructor if you wish students to calculate Vmax. All optical density readings are taken at 540nm.

A sample file for pH and for a solute (either maltose or starch) is included in this repository. Blank templates can be downloaded from the app itself.

By pulling these script files, you will be able to edit the code in R publish the app under your own shiny account. You must change the app name to "app.R" before publishing, otherwise, you will get an error.
