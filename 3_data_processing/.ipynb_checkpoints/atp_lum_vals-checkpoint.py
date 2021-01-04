import numpy as np

def get_atp_conc(lum):
    '''Input lum measurement value with gain 135 from a biotek
    Representative of ATP from Cayman Chem firefly luminescence assay 
    fluor = fluorescence values, input as list in mM units,
    Output: ATP Concentration, list'''
    try:
        # Define variables
        slope = 0.6177860271592077
        intercept = 10.283052798607272
        avg_blank = 11.5
        atp_list = []
        
        for val in lum:
            # Check inputs
            if val <= 11.5:
                print('Inputs must be greater than the lum value for the blanks, 11.5')
                
            # Get log correct sample RLU
            correct_lum = val - avg_blank
            log_lum = np.log(correct_lum)
            
            # Get log atp conc
            log_atp = (log_lum - intercept) / slope
            
            # Get atp conc
            atp = np.exp(log_atp)
            atp_list.append(atp)
            
        return atp_list
        
    except:
        print('Input must be a list and above 11.5')
        
        