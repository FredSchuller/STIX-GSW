;+
;
; :description:
;
;    This procedure returns a simulated model count spectrum using the background
;    continuum and calibration lines generated by stx_bkg_sim_spectrum and the background
;    version of the stix drm. The calibration line strength and detector response
;    parameters (Electron trapping length, hole trapping length and detector resolution) are
;    degraded in a very simple manner using the routine stx_degraded_detector_parameters to represent
;    changes over the mission lifetime
;
; :categories:
;
;    background, calibration, simulation
;
; :params:
;
;    degrade_per : in, required,  type="float"
;                  Fraction of the mission lifetime that has passed
;                  i.e. 0.0 = start of mission 1.0 = 100% expected total mission proton fluence
;
;
; :returns:
;
;    bkg - the background count spectrum in the bins defined by stx_bkg_sim_spectrum
;
; :examples:
;
;    background_counts = stx_degraded_detector_sim(0.5)
;
; :history:
;
;    28-Nov-2017 - ECMD (Graz), initial release
;
;-
function stx_degraded_detector_sim, degrade_per

  stx_degraded_detector_parameters, degrade_per, line_factor = line_factor, continuum_factor = continuum_factor, $
    trap_length_h = trap_length_h, trap_length_e = trap_length_e, func_par = func_par, tail = tail

  configuration_struct = stx_bkg_ecal_spec_config()
  background_sim = stx_bkg_sim_spectrum( edg2,per_sec = 1 , x4 = 0 ,$
    line_factor = line_factor, continuum_factor = continuum_factor, hecht_par = hecht_par, spectrogram = spectrogram_4096, _extra = _extra)
  edge_products, edg2, width = width, edges_1=edg1

  use = where(edg1 gt 1.0)

  dstr = stx_build_drm(edg1[use], func_par = func_par,trap_length_h = trap_length_h, trap_length_e = trap_length_e,$
    /back, tail = tail )

  phm  = dstr.pls_ht_mat

  bkg = phm#background_sim[use[0:-2]]

  return, bkg
end
