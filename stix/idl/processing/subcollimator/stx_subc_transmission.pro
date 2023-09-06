;+
;
; NAME:
;
;   stx_subc_transmission
;
; PURPOSE:
;
;   Compute the trasmission of a STIX subcollimator corrected for internal shadowing
;
; CALLING SEQUENCE:
;
;   subc_transmission = stx_subc_transmission(flare_loc)
;
; INPUTS:
;
;   flare_loc: bidimensional array containing the X and Y coordinate of the flare location 
;             (arcsec, in the STIX coordinate frame)
;
; OUTPUTS:
;
;   A float number that represent the subcollimator transmission value
;
; HISTORY: August 2022, Massa P., first version (working only for detectors 3 to 10)
;
; CONTACT:
;   paolo.massa@wku.edu
;-


function stx_subc_transmission, flare_loc

restore,loc_file( 'grid_temp.sav', path = getenv('STX_GRID') )
fff=read_ascii(loc_file( 'grid_param_front.txt', path = getenv('STX_GRID') ),temp=grid_temp)
rrr=read_ascii(loc_file( 'grid_param_rear.txt', path = getenv('STX_GRID') ),temp=grid_temp)

grid_orient_front = 180.-fff.o ;; Orientation of the slits of the grid as seen from the detector side
grid_pitch_front  = fff.p
grid_slit_front   = fff.slit 
grid_thick_front  = fff.thick
bridge_p_front    = fff.bpitch
bridge_w_front    = fff.bwidth

grid_orient_rear = 180.-rrr.o ;; Orientation of the slits of the grid as seen from the detector side
grid_pitch_rear  = rrr.p
grid_slit_rear   = rrr.slit
grid_thick_rear  = rrr.thick
bridge_p_rear    = rrr.bpitch
bridge_w_rear    = rrr.bwidth
index_f = where(grid_thick_front ne 0.)
mean_thick_f = mean(grid_thick_front[index_f])
index_r = where(grid_thick_rear ne 0.)
mean_thick_r = mean(grid_thick_rear[index_r])

transm = fltarr(32)
sc = fff.sc

for i=0,n_elements(grid_orient_front)-1 do begin
    
 this_list = where(sc eq i+1)
  
  ;; Exclude detectors 1 and 2
  
  if (n_elements(this_list) eq 1) AND (this_list(0) ne -1) then begin
  
  ;if (sc[i] ne 11) and (sc[i] ne 12) and (sc[i] ne 13) and (sc[i] ne 17) and (sc[i] ne 18) and (sc[i] ne 19) then begin
   
  transm_front = stx_grid_transmission(flare_loc[0], flare_loc[1], grid_orient_front[this_list[0]], $
                                       grid_pitch_front[this_list[0]], grid_slit_front[this_list[0]], grid_thick_front[this_list[0]])
  ;transm_front = transm_front * (1.-f_div(bridge_w_front[this_list[0]],bridge_p_front[this_list[0]]))
  
  transm_rear  = stx_grid_transmission(flare_loc[0], flare_loc[1], grid_orient_rear[this_list[0]], $
                                       grid_pitch_rear[this_list[0]], grid_slit_rear[this_list[0]], grid_thick_rear[this_list[0]])
  ;transm_rear  = transm_rear * (1.-f_div(bridge_w_rear[this_list[0]],bridge_p_rear[this_list[0]]))
    
  transm[i] = transm_front * transm_rear
    
  endif
  
  if (n_elements(this_list) eq 2) then begin
 
    transm_front_double_layers = stx_grid_transmission_double_layers(flare_loc[0], flare_loc[1], grid_orient_front[this_list], $
                                                                    grid_pitch_front[this_list], grid_slit_front[this_list], mean_thick_f)
    ;transm_front_double_layers = transm_front_double_layers * (1.-f_div(bridge_w_front[this_list[0]],bridge_p_front[this_list[0]]))

    transm_rear_double_layers = stx_grid_transmission_double_layers(flare_loc[0], flare_loc[1], grid_orient_rear[this_list], $
                                                                    grid_pitch_rear[this_list], grid_slit_rear[this_list], mean_thick_r)
    ;transm_rear_double_layers  = transm_rear_double_layers * (1.-f_div(bridge_w_rear[this_list[0]],bridge_p_rear[this_list[0]]))

    transm[i] = transm_front_double_layers * transm_rear_double_layers

   endif
  
  if (n_elements(this_list) eq 3) then begin

    transm_front_triple_layers = stx_grid_transmission_triple_layers(flare_loc[0], flare_loc[1], grid_orient_front[this_list], $
                                                                    grid_pitch_front[this_list], grid_slit_front[this_list], mean_thick_f)
    ;transm_front_triple_layers = transm_front_triple_layers * (1.-f_div(bridge_w_front[this_list[0]],bridge_p_front[this_list[0]]))

    transm_rear_triple_layers = stx_grid_transmission_triple_layers(flare_loc[0], flare_loc[1], grid_orient_rear[this_list], $
                                                                    grid_pitch_rear[this_list], grid_slit_rear[this_list], mean_thick_r)
    ;transm_rear_triple_layers  = transm_rear_triple_layers * (1.-f_div(bridge_w_rear[this_list[0]],bridge_p_rear[this_list[0]]))

    transm[i] = transm_front_triple_layers * transm_rear_triple_layers
 
  endif  
  
endfor

transm[where(transm eq 0.)] = 1

return, transm

end