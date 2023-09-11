;+
;
; NAME:
;
;   stx_triple_layers_grid_values
;
; PURPOSE:
;
;   Compute the grids 1 values (pitch, slit, slat and orientation angle) to treat them as a single layer grid
;
; CALLING SEQUENCE:
;
;   average_values_f = stx_triple_layers_grid_values(grid_pitch_first_layer, grid_pitch_second_layer, grid_pitch_third_layer, $
;                                                     grid_slit_first_layer, grid_slit_second_layer, grid_slit_third_layer, $
;                                                     grid_orient_first_layer, grid_orient_second_layer, grid_orient_third_layer, $
;                                                     grid_phase_first_layer, grid_phase_second_layer, grid_phase_third_layer)
;
;
; INPUTS:
;     grid_pitch_first_layer: dimension of the first layer pitch (mm)
;     grid_pitch_second_layer: dimension of the second layer pitch (mm)
;     grid_pitch_third_layer: dimension of the third layer pitch (mm)
;
;     grid_slit_first_layer: dimension of the first layer slit (mm)
;     grid_slit_second_layer: dimension of the second layer slit (mm)
;     grid_slit_third_layer: dimension of the third layer slit (mm)
;
;     grid_orient_first_layer: orientation angle of the slit of first layer (looking from detector side, in degrees)
;     grid_orient_second_layer: orientation angle of the slit of second layer (looking from detector side, in degrees)
;     grid_orient_third_layer: orientation angle of the slit of third layer (looking from detector side, in degrees)
;
;     grid_phase_first_layer: phase of the first layer (mm)
;     grid_phase_second_layer: phase of the second layer (mm)
;     grid_phase_third_layer: phase of the third layer (mm)
;
; OUTPUTS:
;
;     pitch, slit, slat and orientation angle of the double-layers grid treated as a single layer grid
;
; HISTORY: September 2022, Volpara A., created
;
; CONTACT:
;   volpara@dima.unige.it
;-


function stx_triple_layers_grid_values, grid_pitch_first_layer, grid_pitch_second_layer, grid_pitch_third_layer, $
                                                    grid_slit_first_layer, grid_slit_second_layer, grid_slit_third_layer, $
                                                    grid_orient_first_layer, grid_orient_second_layer, grid_orient_third_layer, $
                                                    grid_phase_first_layer, grid_phase_second_layer, grid_phase_third_layer

  this_orientation = average([grid_orient_first_layer, grid_orient_second_layer, grid_orient_third_layer])

  slat_first_layer = grid_pitch_first_layer - grid_slit_first_layer
  slat_second_layer = grid_pitch_second_layer - grid_slit_second_layer
  slat_third_layer = grid_pitch_third_layer - grid_slit_third_layer
  this_slat = average([slat_first_layer, slat_second_layer, slat_third_layer])

  average_pitch = average([grid_pitch_first_layer, grid_pitch_second_layer, grid_pitch_third_layer])
  average_slit = average([grid_slit_first_layer, grid_slit_second_layer, grid_slit_third_layer])

  this_slit = (average_slit - 2*this_slat)/3.
   
;  this_pitch = average_pitch/3.
;  this_slit = this_pitch - this_slat

  this_pitch = this_slit + this_slat
  
  average_phase = (grid_phase_first_layer + grid_phase_second_layer + grid_phase_third_layer + average_pitch/3. + 2.*average_pitch/3.)/3.

  this_phase = (average_phase + average_pitch/6.) mod (average_pitch/3.)

  grid3_values = { pitch: this_pitch, orientation: this_orientation, phase: this_phase, slit: this_slit};, phase: this_phase}

  return, grid3_values

end
