FUNCTION make_disc, diam, res=res, center=center, fixedn=fixedn, edge=edge, transmission=transmission
   ;FUNCTION to return an image of a rectangle with given aspect ratio.
   ;Diam and res in microns
   ;phimajor is angle of major axis in radians
   ;AB 5/2012
   
   IF n_elements(res) eq 0 THEN res=1.0
   IF n_elements(fixedn) eq 0 THEN fixedn=0
   IF n_elements(edge) eq 0 THEN edge=0
   IF n_elements(transmission) eq 0 THEN transmission=0
   
   diamfloat=float(diam)
   
   ;Find n based on diameter or user-specified
   IF fixedn eq 0 THEN n=round(diamfloat/res*1.5)>2 ELSE n=fixedn
   img=fltarr(n,n)+1
   center=(n-1)/2.0   ;Centerpoint of image
   
   d=fltarr(n,n)           ;Distance from center
   im=fltarr(n,n)+1
   remainder=n mod 2     ;If n is odd then need to add this to the length of x
   x=[reverse(findgen(n/2.0+remainder)),findgen(n/2.0)+1]*res
   FOR i=0,n-1 DO d[i,*]=sqrt(x[i]^2+x^2)
     
   good=where((d le diamfloat/2.0),ngood)
   IF ngood gt 0 THEN im[good]=transmission
   
   IF edge ne 0 THEN BEGIN
      ;Make a transparent image with 'edge' thickness around border
      good=where((d le diamfloat/2.0-edge),ngood)
      IF ngood gt 0 THEN im[good]=1 
   ENDIF
  
   return,im
END
   
   