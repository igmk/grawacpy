def rh_to_qabs(rh,t):
    '''
    conversion of relative humidity[%] to absolute humidity[kg/m^3]:
    INPUT:
        relative humidity   (rh)      [%]
        temperature         (T)       [K]
    OUTPUT:
        absolute humitity   (qabs)    [kg/m^3]

    calculates saturation pressure according to Goff and Gratch
    '''
    
    rh = np.asarray(rh)
    t = np.asarray(t)
    
    psat=1013.246 * 10**( -7.90298*(373.16/t-1) \
      + 5.02808*np.log10(373.16/t) \
      - 1.3816e-7*(10.**(11.344*(1.-t/373.16))-1.) \
      + 8.1328e-3 * (10.**(-3.49149*(373.16/t-1.))-1.) )
    
    qabs = rh * psat / (461.5*t)
    
    return qabs
