FUNCTION make_dendrite, diam, aspectratio, phimajor, res=res, center=center, fixedn=fixedn, edge=edge, transmission=transmission
   ;FUNCTION to return an image of a dendrite with given orientations.
   ;Diam and res in microns
   ;psi is angle in/out of image
   ;phimajor is angle of major axis in radians
   ;AB 5/2012
   
   IF n_elements(res) eq 0 THEN res=1.0
   IF aspectratio lt 4 THEN stop,'Suggest higher aspect ratio for dendrite arms'
   IF n_elements(fixedn) eq 0 THEN fixedn=0
   IF n_elements(edge) eq 0 THEN edge=0
   IF n_elements(transmission) eq 0 THEN transmission=0
  
   ;Need this geometrical correction factor since we want diam defined as maximum chord
   columncorrect=aspectratio/sqrt(1+float(aspectratio)^2)  

   diamfloat=float(diam)*columncorrect
   
   ;Find n based on diameter or user-specified
   IF fixedn eq 0 THEN n=round(diamfloat/res*1.5)>2 ELSE n=fixedn
   img=fltarr(n,n)+1
   center=(n-1)/2.0   ;Centerpoint of image
   
   d=fltarr(n,n)           ;Distance from center
   
   ;Use line in form of Ax+By=0 (goes through origin)
   theta=!pi/3
   a1=cos(phimajor)
   b1=sin(phimajor)
   a2=cos(phimajor+theta)
   b2=sin(phimajor+theta)
   a3=cos(phimajor-theta)
   b3=sin(phimajor-theta)

   den1=sqrt(a1^2+b1^2)  ;Denominator in distance formula
   den2=sqrt(a2^2+b2^2)  ;Denominator in distance formula
   den3=sqrt(a3^2+b3^2)  ;Denominator in distance formula

   xpos=fltarr(n,n)
   ypos=fltarr(n,n)
   FOR i=0,n-1 DO BEGIN
      xpos[i,*]=i-n/2.0
      ypos[*,i]=i-n/2.0
   ENDFOR
   
   ;Find the intersection of three long columns to make the plate
   distmajor1=abs(a1*xpos*res +  b1*ypos*res)/den1
   distmajor2=abs(a2*xpos*res +  b2*ypos*res)/den2
   distmajor3=abs(a3*xpos*res +  b3*ypos*res)/den3
   distminor1=abs(b1*xpos*res -  a1*ypos*res)/den1
   distminor2=abs(b2*xpos*res -  a2*ypos*res)/den2
   distminor3=abs(b3*xpos*res -  a3*ypos*res)/den3
   
   
   
   ;The aspect ratio will be used to decrease width of axis-1
   ;good=where((distmajor1 le diamfloat*platecorrect/2.0*cos(psi)) or (distmajor2 le diamfloat/2.0*sin(theta)) or (distmajor3 le diamfloat/2.0*sin(theta)),ngood)
   
   good=where(((distmajor1 le diamfloat/2.0) and (distminor1 le (diamfloat/aspectratio/2.0))) or $
              ((distmajor2 le diamfloat/2.0) and (distminor2 le (diamfloat/aspectratio/2.0))) or $
              ((distmajor3 le diamfloat/2.0) and (distminor3 le (diamfloat/aspectratio/2.0))), ngood)
   IF ngood gt 0 THEN img[good]=transmission
   
   IF edge ne 0 THEN BEGIN
       print, 'Edge option not implemented in make_dendrite.pro', edge
   ;   ;Make a transparent plate with 'edge' thickness around border
   ;   good=where((distmajor1 le diamfloat*platecorrect/2.0*cos(psi) -edge) and (distmajor2 le diamfloat/2.0*sin(theta)-edge) and (distmajor3 le diamfloat/2.0*sin(theta)-edge),ngood)
   ;   IF ngood gt 0 THEN img[good]=1 
   ENDIF

   return,img
END
   
   