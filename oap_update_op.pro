PRO oap_update_op, op
   ;Updates the 'op' structure if any new tags are added
   ;All units in meters, seconds
   
   ;-----------------Model settings-------------------------
   ;Number of particles to model
   IF total(where(tag_names(op) eq 'N')) eq -1 THEN BEGIN & op=create_struct(op,'n',10000) & print, 'Default N set.' & ENDIF
   
   ;TAS   
   IF total(where(tag_names(op) eq 'TAS')) eq -1 THEN BEGIN & op=create_struct(op,'tas',100.0) & print, 'Default TAS set.' & ENDIF
   
   ;Lambda (1/m), exponential
   IF total(where(tag_names(op) eq 'LAM')) eq -1 THEN BEGIN & op=create_struct(op,'lam',1.0e4) & print, 'Default lambda set.' & ENDIF
   
   ;N0 (1/m4), exponential
   IF total(where(tag_names(op) eq 'N0')) eq -1 THEN BEGIN & op=create_struct(op,'n0',1.0e8) & print, 'Default N0 set.' & ENDIF
   
   ;Mu, default 0
   IF total(where(tag_names(op) eq 'MU')) eq -1 THEN BEGIN & op=create_struct(op,'mu',0.0) & print, 'Default mu set.' & ENDIF
   
   ;All particles a single size
   IF total(where(tag_names(op) eq 'MONODISPERSE')) eq -1 THEN BEGIN & op=create_struct(op,'monodisperse',0.0) & print, 'Monodisperse off.' & ENDIF
   
   ;Add some 'wobble' to the monodisperse value
   IF total(where(tag_names(op) eq 'MONODISPERSERANGE')) eq -1 THEN BEGIN & op=create_struct(op,'monodisperserange',0.0) & print, '' & ENDIF
   
   ;Padding to make the particle box bigger than array by this amount on each side
   IF total(where(tag_names(op) eq 'PAD')) eq -1 THEN op=create_struct(op,'pad',0.01)
   
   ;Image output
   IF total(where(tag_names(op) eq 'WRITEIMAGES')) eq -1 THEN op=create_struct(op,'writeimages',0)
   
   ;Default particle shape
   IF total(where(tag_names(op) eq 'SHAPE')) eq -1 THEN BEGIN & op=create_struct(op,'shape','disc') & print, 'Default shape set to disc.' & ENDIF
   
   ;Default aspect ratio (ice particles only)
   IF total(where(tag_names(op) eq 'COLUMNASPECT')) eq -1 THEN BEGIN & op=create_struct(op,'columnaspect',2.0) & print, 'Default aspect ratio set.' & ENDIF
   
   ;Default transparency (grey level set on original particle)
   IF total(where(tag_names(op) eq 'TRANSPARENCY')) eq -1 THEN op=create_struct(op,'transparency',0)
   
   ;Spinning disk mode (x position fixed)
   IF total(where(tag_names(op) eq 'SPINNINGDISK')) eq -1 THEN op=create_struct(op,'spinningdisk',0)
   
   ;Shattering simulation on/off
   IF total(where(tag_names(op) eq 'SHATTER')) eq -1 THEN op=create_struct(op,'shatter',0)
   
   ;Smoothing window in meters
   IF total(where(tag_names(op) eq 'SMOOTHPARTICLE')) eq -1 THEN op=create_struct(op,'smoothparticle',0)
   
   ;Minimum size in meters, skip all particles smaller than this
   IF total(where(tag_names(op) eq 'MINSIZE')) eq -1 THEN op=create_struct(op,'minsize',0)
   
   ;File name for bitmap shapes
   IF total(where(tag_names(op) eq 'BITMAPFILE')) eq -1 THEN op=create_struct(op,'bitmapfile',strarr(n_elements(op.monodisperse)))
   
   ;Put all particles in the center of the DOF
   IF total(where(tag_names(op) eq 'CENTERZ')) eq -1 THEN op=create_struct(op,'centerz',0)
   
   ;Put all particles at a fixed position in DOF relative to the center
   IF total(where(tag_names(op) eq 'FIXEDZ')) eq -1 THEN op=create_struct(op,'fixedz',0)
   
   
   
   ;-----------------Probe settings, CAPS default------------------------
   IF total(where(tag_names(op) eq 'NDIODES')) eq -1 THEN BEGIN & op=create_struct(op,'ndiodes',64) & print, 'Default nDiodes set.' & ENDIF
   ;Resolution in meters
   IF total(where(tag_names(op) eq 'RES')) eq -1 THEN BEGIN & op=create_struct(op,'res',25.0e-6) & print, 'Default resolution set.' & ENDIF
   ;Arm width in meters
   IF total(where(tag_names(op) eq 'ARMWIDTH')) eq -1 THEN BEGIN & op=create_struct(op,'armwidth',0.1) & print, 'Default arm width set.' & ENDIF
   ;;Baumgardner (1997 JTECH) time response 
   IF total(where(tag_names(op) eq 'TAU')) eq -1 THEN BEGIN & op=create_struct(op,'tau',0.4e-6) & print, 'Default tau set.' & ENDIF
   ;Diode sensitivity
   IF total(where(tag_names(op) eq 'SENSITIVITY')) eq -1 THEN BEGIN & op=create_struct(op,'sensitivity',fltarr(512)+1) & print, 'Default sensitivity set.' & ENDIF
   ;TAS error factor for particle stretching
   IF total(where(tag_names(op) eq 'TASERROR')) eq -1 THEN BEGIN & op=create_struct(op,'taserror',1.0) & print, 'Default TAS error set.' & ENDIF
   ;Skip slice flag for legacy OAPs
   IF total(where(tag_names(op) eq 'SKIPSLICE')) eq -1 THEN BEGIN & op=create_struct(op,'skipslice',0) & print, 'Default slice skip flag set.' & ENDIF
   ;Wavelength in meters
   IF total(where(tag_names(op) eq 'WAVELENGTH')) eq -1 THEN BEGIN & op=create_struct(op,'wavelength',0.658e-6)  & print, 'Default wavelength set.' & ENDIF
   ;Greyscale shadow thresholds, lightest to darkest
   IF total(where(tag_names(op) eq 'GREYTHRESHOLDS')) eq -1 THEN BEGIN & op=create_struct(op,'greythresholds',[0.25, 0.5, 0.75]) & print, 'Default grey thresholds set.' & ENDIF
   ;Diode spacing ratio (x/y), fraction of diode space that is actually active
   IF total(where(tag_names(op) eq 'DIODESPACING')) eq -1 THEN BEGIN & op=create_struct(op,'diodespacing',0.8) & print, 'Default diode spacing set.' & ENDIF
   ;Deprecated ;IF total(where(tag_names(op) eq 'DOFTHRESHOLD')) eq -1 THEN BEGIN & op=create_struct(op,'dofthreshold',0) ;Shadow level (%) required to register a particle

   ;Check for some incompatibilites and give a warning
   IF (op.smoothparticle ne 0) and (op.shape eq 'disc') THEN print,'NOTE: Particle smoothing will not work for this setup, use a Fresnel shape.'


END
