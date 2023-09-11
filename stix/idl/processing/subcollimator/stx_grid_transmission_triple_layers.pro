;+
;
; NAME:
;
;   stx_grid_transmission_triple_layers
;
; PURPOSE:
;
;   Compute the trasmission of a STIX two_layers grid corrected for internal shadowing
;
; CALLING SEQUENCE:
;
;   int_shadow = stx_grid_transmission_triple_layers(x_flare, y_flare, grid_orient, grid_pitch, grid_slit, grid_thick)
;
; INPUTS:
;
;   x_flare: X coordinate of the flare location (arcsec, in the STIX coordinate frame)
;
;   y_flare: Y coordinate of the flare location (arcsec, in the STIX coordinate frame)
;
;   grid_orient: orientation angle of the slits of each layer of the grid (looking from detector side, in degrees)
;
;   grid_pitch: dimension of pitches of each layer of the grid (mm)
;
;   grid_slit: dimension of slits of each layer of the grid (mm)
;
;   grid_thick: thickness of the grid (mm)
;
; OUTPUTS:
;
;   A float number that represent the grid transmission value
;
; HISTORY: October 2022, Volpara A., created
;
; CONTACT:
;   volpara@dima.unige.it
;-


function stx_grid_transmission_triple_layers, x_flare, y_flare, grid_orient, grid_pitch, grid_slit, mean_thick

  average_grid_orient = average(grid_orient)

  grid_slat = grid_pitch - grid_slit
  average_slat = average(grid_slat)
  average_slit = average(grid_slit)
  average_pitch = average(grid_pitch)

  this_slat = average_slat
  
  this_pitch = average_pitch/3.
  
  this_slit = this_pitch - this_slat

  ;; Distance of the flare on the axis perpendicular to the grid orientation
  flare_dist   = abs(x_flare * cos(average_grid_orient * !dtor) + y_flare * sin(average_grid_orient * !dtor))

  ;; Internal shadowing
  shadow_width = (mean_thick/3.)  * tan(flare_dist / 3600. * !dtor)

  return, (this_slit - shadow_width) / this_pitch


end
