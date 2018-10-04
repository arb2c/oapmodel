FUNCTION baumgardner_delay,intensity,res,tas,tau=tau,xpos=xpos,cleanstart=cleanstart
   ;FUNCTION to compute the time delay of OAP diodes based on 
   ;Baumgardner and Korolev JTECH 1997
   ;Res, tas, and tau all in seconds and meters
   ;Cleanstart sets the first element of the array to 1.0, this is to avoid noisy starts.
   ;Intensity can be in any unit.
   ;AB 5/2012
   
   IF n_elements(tau) eq 0 THEN tau=0.4e-6 ;Default time response
   IF n_elements(cleanstart) eq 0 THEN cleanstart=0
   IF (res gt 3e-3) or (tau[0] gt 1e-3) THEN BEGIN
      print, 'Res/tau too big, enter in meters/seconds', res, tau
      return,0
   ENDIF
   localtau=tau
   s=size(intensity)
   
   IF (n_elements(xpos) gt 0) and (n_elements(tau) gt 1) THEN BEGIN
      ;This allows for a flexible tau depending on the diode
      ;Will use xpos to decide which pixels get which tau
      xpos=round(xpos)
      localtau=fltarr(s[1])
      FOR i=0,n_elements(xpos)-1 DO localtau[xpos[i]>0<(s[1]-1):*]=tau[i]
   ENDIF

   deltat=float(res)/tas
   term1=(1-exp(-deltat/localtau))
   r=intensity  ;Start with unaltered array 
   
   ;Check for a 2-D image, otherwise assume linear pattern
   IF s[0] eq 1 THEN BEGIN      
      ;Linear version
      IF cleanstart eq 1 THEN r[0]=1
      FOR i=1,n_elements(r)-1 DO BEGIN
         ;Time response based on difference in intensity over difference in time
         r[i]=r[i-1]-(r[i-1]-intensity[i])*term1
      ENDFOR
   ENDIF ELSE BEGIN
      ;2-D version
      ;Faster version, do entire slices at a time
      IF cleanstart eq 1 THEN r[*,0]=1
      FOR i=1,s[2]-1 DO BEGIN
         ;Time response based on difference in intensity over difference in time
         r[*,i]=r[*,i-1]-(r[*,i-1]-intensity[*,i])*term1
      ENDFOR
   ENDELSE
   
   return,r
END
