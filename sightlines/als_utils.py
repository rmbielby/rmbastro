import numpy as np

def curve_of_growth(ew):
    from scipy import interpolate
    CoG = np.array([[11.419   , 11.892   , 12.271   , 12.634   , 13.123   , 13.422   ,
        14.036   , 14.822   , 16.426   , 17.291   , 17.716   , 17.952   ,
        18.251   , 18.834   , 19.543   , 20.125   , 20.676   , 20.991   ],
       [-2.8301  , -2.3699  , -1.9876  , -1.6407  , -1.2655  , -1.0389  ,
        -0.81947 , -0.69204 , -0.5292  , -0.39469 , -0.2531  , -0.15398 ,
        -0.026549,  0.2354  ,  0.59646 ,  0.85133 ,  1.0425  ,  1.1416  ]])
    cogfunc = interpolate.interp1d(CoG[1,:],CoG[0,:])
    coldens = cogfunc(np.log10(ew))
    return coldens
