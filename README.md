# PZS_analysis
An alternative FFCF extraction method - Matlab Code

This Matlab code evaluates FFCF from 2DIR spectra of a single or multi-component system (overlapping spectral bands - resolved or unresolved).

Datasets : "Single_Transition_Simulated_2DIR", "Two_Transition_Resolved_Simulated_2DIR" and "Two_Transition_Unresolved_Simulated_2DIR"

# INPUT PARAMETERS
  folder (required) -> Foldername containing 2DIR spectra
  
                        Matlab Figure Format
                        
                        FigureName -'Tw_()' where () is filled with waiting time in femtoseconds
                        
                        FigureAxes - wtau (y-direction) and wt (x-direction)
  
  N (required)-> Different analysis
  
                'single' or 'overlapping'

  freq_r (optional) -> Frequency window parameters
  
                     scalar input - [radius] (single component anaysis)
                     
                     vector input - [wtL wtU wtauL wtauU] (multi-component analysis) 
                     
                        wtL - lower bound for wt frequency 
                        
                        wtU - upper bound for wt frequency
                        
                        wtauL - lower bound for wtau frequency
                        
                        wtauU - upper bound for wtau frequency
                        
                      'Default' = []
  
  norder (optional) -> nth order pseudo-Zernike moment (larger than 4 is recommended)
  
                        'Default' = 15
                        
  nExpDec (optional) -> type of exponential decay fit
  
                          single exponential decay fit (1)
                          
                          area under the curve for multi-exponential decay fit (0)
                          
                          'Default' = 1
 
  datatype (optional)-> center peak adjustment for multi-component system
  
                          'psc' - peak shift correction
                          
                            ''  - no correction for peak shift
                            
                            'Default' = ''
                            
  # OUTPUT
       data_2D -> spectral data extracted from 2DIR matlab figures
       
       output -> struct file ('spectral_region','analysis','results','figures')
       
                'single'
                  
                    Selected radius of frequency window 

                    Selected TPZS vs. waiting time

                    CLS vs. TPZS

                    CL-PZS
                  
                'overlapping'
                  
                    FR-PZS

                    Selected FR-PZS vs. waiting time for each component
                
                
 # DEMO 
    Default parameters are used for following examples:
    
      First, copy both .m file (matlab code) into the folder containg .txt file (2DIR data)

      Second, txt2fig_2DIR.m file is executed to convert 2DIR .txt file to matlab figures and save in a new folder.

      Then, PZS.m file is executed for PZS analysis.
    
            Single component analysis - Download "Single_Transition_Simulated_2DIR"

              txt2fig_2DIR('Simulation1')

              [data_2D,output]= PZS('Simulation1','single')

            Two component analysis (resolved) - Download "Two_Transition_Resolved_Simulated_2DIR"

              txt2fig_2DIR('Simulation2')

              [data_2D,output]= PZS('Simulation2','overlapping')

            Two component anaylsis (unresolved) - Download "Two_Transition_Unresolved_Simulated_2DIR"

              txt2fig_2DIR('Simulation3')

              [data_2D,output]= PZS('Simulation3','overlapping')

      
