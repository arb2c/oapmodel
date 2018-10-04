FUNCTION korolev_fresnel, diam, z, wavelength=wavelength, makeplot=makeplot, littler=littler, microns=microns
   ;Computes Fresnel diffraction pattern following Korolev 1991 and 1998.
   ;All input units in meters
   ;'z' is distance from center of focus
   ;littler specifies radius points to compute, otherwise makes a cross-section
   
   IF keyword_set(microns) THEN BEGIN
      ;Units entered in microns, make conversions to meters
      diam=diam/1.0e6
      z=z/1.0e6
      IF keyword_set(wavelength) THEN wavelength=wavelength/1.0e6
   ENDIF
   
   IF n_elements(wavelength) eq 0 THEN wavelength=0.65e-6 ;meters
   IF n_elements(makeplot) eq 0 THEN makeplot=0 
   IF (diam gt 1) or (z gt 1) THEN BEGIN
      print,'All units in meters'
      return,0
   ENDIF
   
   bigr=diam/2.0             ;radius
   IF n_elements(littler) eq 0 THEN BEGIN
      n=500  ;Number of zd points to compute, range -5 to +5
      roverR_forplots=findgen(n)/(n/10) -5 ;Normalized distance from center of particle
      roverR=abs(roverR_forplots)  ;ignore negative r/R, just for plotting convenience
      littler=abs(roverR*bigr)  ;actual distance from center of particle
   ENDIF ELSE BEGIN
      n=n_elements(littler)
      roverR_forplots=littler/bigr
      roverR=abs(roverR_forplots)  ;ignore negative r/R, just for plotting convenience
      littler=abs(littler)
   ENDELSE
   
  
   k=2*!pi/wavelength        ;wavenumber
   zd=z*wavelength/bigr^2    ;dimensionless variable
   IF zd eq 0 THEN zd=1e-5   ;avoid divide by zero errors when z=0
   
   nalpha=500  ;Number of angles to integrate over
   alpha=findgen(nalpha)/nalpha*2*!pi
   dalpha=alpha[1]-alpha[0]
   
   ua=complexarr(n)
   ;ub98=complexarr(n) ;Using korolev 1998 formulation
   ub=complexarr(n)+1  ;Using korolev 1991 formulation
   coeff=-exp(complex(0,k*z))/(2*!pi)
   uadefault=exp(complex(0,k*z))
   w=where(roverr le 7,nw)   ;Only compute within roverr lower than 7
   FOR j=0L,nw-1 DO BEGIN   ;n-1 DO BEGIN
      i=w[j]
      ;Compute Ua
      IF (littler[i] gt bigr) THEN ua[i]=uadefault ELSE ua[i]=0
      
      ;Compute Ub
      ;s=sqrt(bigr^2 + z^2 + littler[i]^2 -2*bigr*littler[i]*cos(alpha))   
      ;ub98[i]=  -1.0/(4*!pi) * total(exp(complex(0,k*s)) * (bigr^2 - bigr*littler[i] * cos(alpha)) *dalpha / (z*(s-z)))
      
      ;Ub 1991 formula
      term1=exp(complex(0,!pi)/zd * (1+roverr[i]^2-2*roverr[i]*cos(alpha)))
      term2=1-roverr[i]*cos(alpha)
      term3=1+roverr[i]^2-2*roverr[i]*cos(alpha)         
      ub[i]= coeff*total(term1*term2*dalpha/term3)
   ENDFOR
   
   ;intensity98=abs(ub98-ua)^2
   intensity=abs(ub-ua)^2
   
   bad=where(abs(roverr-1) lt 0.01,nbad)   ;Discontinuties here, nan them
  
   IF nbad gt 0 THEN BEGIN
      good=where(abs(roverr-1) ge 0.01,nbad)
      intensity[bad]=interpol(intensity[good],good,bad) ;!values.f_nan
   ENDIF
   
   IF makeplot THEN BEGIN
      plot,roverR_forplots,intensity
      ;oplot,roverR,intensity98,color=100
   ENDIF
   return,{roverr:roverr, intensity:intensity, zd:zd}
END
