# PZS_analysis
An alternative FFCF extraction method - Matlab Code



This Matlab code evaluates FFCF from 2DIR spectra of a single or multi-component system (overlapping spectral bands - resolved or unresolved).

Datasets : "Single_Transition_Simulated_2DIR", "Two_Transition_Resolved_Simulated_2DIR" and "Two_Transition_Unresolved_Simulated_2DIR"

# INPUT PARAMETERS
  [data_2D,output]=PZS(N,freq_r,norder,nExpDec,datatype) 
  
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
                          
                          multi-exponential decay fit (>1)
                          
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
    
      First, txt2fig_2DIR.m file is executed to convert 2DIR .txt file to matlab figures and save in a new folder.

      Second, PZS.m file is executed for PZS analysis.                   
    
                        Single component analysis - Download "Single_Transition_Simulated_2DIR"

                          txt2fig_2DIR('Simulation1')

                          [data_2D,output]= PZS('single') - Select the 'Simulation1' folder when the dialog box pops up.

                        Two component analysis (resolved) - Download "Two_Transition_Resolved_Simulated_2DIR"

                          txt2fig_2DIR('Simulation2')

                          [data_2D,output]= PZS('overlapping') - Select the 'Simulation2' folder when the dialog box pops up.

                        Two component anaylsis (unresolved) - Download "Two_Transition_Unresolved_Simulated_2DIR"

                          txt2fig_2DIR('Simulation3')

                          [data_2D,output]= PZS('overlapping') - Select the 'Simulation3' folder when the dialog box pops up.
                          
      FR-PZS analysis with tunable parameters (One Example)
      
            [data_2D,output]= PZS('overlapping',[2200 2240 2215 2235],15,1,'')
              
              N = 'overlapping'
              
              freq_r = [2200 2240 2215 2235]
              
              norder = 15
              
              nExpDec = 1 (single exponential decay fit)
              
              datatype = '' (no peak shift corretion)
      

      
