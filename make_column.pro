FUNCTION make_column, diam, aspectratio, phimajor, res=res, center=center, fixedn=fixedn, edge=edge, $
   transmission=transmission, areapercent=areapercent, holetransmission=holetransmission
   ;FUNCTION to return an image of a rectangle with given aspect ratio.
   ;Diam and res in microns
   ;phimajor is angle of major axis in radians
   ;AB 5/2012
   
   IF n_elements(res) eq 0 THEN res=1.0
   IF aspectratio lt 1 THEN stop,'Must have aspect ratio above 1'
   IF n_elements(fixedn) eq 0 THEN fixedn=0
   IF n_elements(edge) eq 0 THEN edge=0
   IF n_elements(transmission) eq 0 THEN transmission=0
   IF n_elements(holetransmission) eq 0 THEN holetransmission=1.0
   IF n_elements(areapercent) eq 0 THEN areapercent=0
   
   ;Need this geometrical correction factor since we want diam defined as maximum chord
   columncorrect=aspectratio/sqrt(1+float(aspectratio)^2)  

   diamfloat=float(diam)*columncorrect
   
   ;Find n based on diameter or user-specified
   IF fixedn eq 0 THEN n=round(diamfloat/res*1.5)>2 ELSE n=fixedn
   img=fltarr(n,n)+1
   center=(n-1)/2.0   ;Centerpoint of image
   
   d=fltarr(n,n)           ;Distance from center
   
   ;Use line in form of Ax+By=0 (goes through origin)
   a=cos(phimajor)
   b=sin(phimajor)
   
   den=sqrt(a^2+b^2)  ;Denominator in distance formula

   xpos=fltarr(n,n)
   ypos=fltarr(n,n)
   FOR i=0,n-1 DO BEGIN
      xpos[i,*]=i-n/2.0
      ypos[*,i]=i-n/2.0
   ENDFOR
   
   distmajor=abs(a*xpos*res +  b*ypos*res)/den
   distminor=abs(b*xpos*res -  a*ypos*res)/den
   
   good=where((distmajor le diamfloat/2.0) and (distminor le (diamfloat/aspectratio/2.0)),ngood)
   IF ngood gt 0 THEN img[good]=transmission
    
   IF edge ne 0 THEN BEGIN
      ;Make a transparent column with 'edge' thickness around border
      good=where((distmajor le diamfloat/2.0-edge) and (distminor le (diamfloat/aspectratio/2.0)-edge),ngood)
      IF ngood gt 0 THEN img[good]=holetransmission 
   ENDIF
   
   IF areapercent ne 0 THEN BEGIN
      ;Make a transparent column with 'areapercent' area removed from the center
      
      ;This method results in unequal border thickness
      ;factor=1.0/(areapercent/100.0)
      ;good=where((distmajor le diamfloat/2.0/sqrt(factor)) and (distminor le (diamfloat/aspectratio/2.0/sqrt(factor))),ngood)
      ;IF ngood gt 0 THEN img[good]=holetransmission
   
      ;This method results in equal border thickness, use quadratic equation to solve
      dmaj=diamfloat
      dmin=diamfloat/aspectratio
      dmaj_inner=0.5*(dmaj-dmin + sqrt((dmaj-dmin)^2+4*dmaj*dmin*areapercent/100.0))
      dmin_inner=dmin-dmaj+dmaj_inner
      good=where((distmajor le dmaj_inner/2.0) and (distminor le (dmin_inner/2.0)),ngood)
      IF ngood gt 0 THEN img[good]=holetransmission
   ENDIF
 
   return,img
END
   
   