;---------------------------------------------------------------------------
;+
; PROJECT:
;       STIX
;
; NAME:
;       stx_build_drm
;
; PURPOSE:
;    Calculates the STIX detector response matrix with respect to the input energy bins
;    STIX detector parameters are predefined in this program and cannot be set from outside
;
;
; CATEGORY:
;       helper methods
;
; CALLING SEQUENCE:
;    IDL> stx_drm = stx_calc_drm(ct_edges1)
;    IDL> help, stx_drm,/st
;** Structure <c7b178>, 17 tags, length=13184, data length=13178, refs=1:
;   DETECTOR        STRING    'cdte'
;   AREA            FLOAT           1.00000
;   FUNC            STRING    'stx_fwhm'
;   FUNC_PAR        FLOAT     Array[4]
;   D               FLOAT          0.100000
;   Z               INT             13
;   GMCM            FLOAT          0.270000
;   NFLUX           LONG                32
;   ELOSS_MAT       FLOAT     Array[32, 32]
;   PLS_HT_MAT      FLOAT     Array[32, 32]
;   E_2D            FLOAT     Array[2, 32]
;   SMATRIX         FLOAT     Array[32, 32]
;   EDGES_IN        FLOAT     Array[33]
;   EDGES_OUT       FLOAT     Array[33]
;   EMEAN           FLOAT     Array[32]
;   EWIDTH          FLOAT     Array[32]
;   INFO            STRING    Array[1, 3]
;
;
; HISTORY:
;
;       29-Apr-2013 - ines.kienreich@uni-graz.at (INK)
;       26-jul-2013 - richard.schwartz@nasa.gov added func_par, more to come
;       6-mar-2014  - richard.schwartz@nasa.gov, added keywords,added efficiency to output structure
;       9-mar-2015  - richard.schwartz@nasa.gov, changed func_par for stx_fwhm
;       now it is only the electronic noise, Fano factor in stx_fwhm_fano and Hall tailing to be included separately
;       Hall tailing computation should be prior to resolution broadening
;       22-apr-2015  - ECMD (Graz), hole tailing due to incomplete charge collection applied to energy loss matrix
;       06-oct-2015  - richard.schwartz@nasa.gov, using _extra for keyword inheritance to resp_calc
;       14-jun-2017  - ECMD (Graz), Replaced stx_calc_pulse_tailing with stx_tailing_matrix for calculation of hole tailing
;       10-oct-2017  - ECMD (Graz), Background keyword added, if set DRM with no attenuation and is calculated large depth 
;       is used as it is assumed all model background photons will generate a count. For use with model background spectra 
;       such as those generated by stx_bkg_continuum_mdl  
;       29-nov-2017  - richard.schwartz@nasa.gov, if background is set then a background field is set
;       ;in the output structure
;
;-

;+
; :description:
;    computes the energy loss matrix, pulse height matrix and finally the detector response matrix
;    for a given energy binning; As per definition the counts and photons have the same energy binning,
;    the returned matrices are (n_ebins,n_ebins) arrays
;
; :keywords:
;
; :params:
;   ct_energy_edges - 1dim floating point array of energy edges of e-bins (for STIX 33 elements array from 4 to 150 keV)
; ;keywords:
;  PH_ENERGY_EDGES
;  ATTENUATOR = 0 - no attenuator, 1 - include the attenuator
;  D_AL - thickness of Al attenuator, default, 600 microns
;  D_BE - effective thickness of Be window, default, 3.5 mm
;
;  FUNC_PAR - parameters used with fwhm function, experts only
;  FUNC - function for resolution broadening default is stx_fwhm
;  EFFICIENCY - photopeak efficiency, multiply by area to get counts and not counts/cm2
;  TAILING - Include calculation of tailing due to incomplete charge collection
;  VERBOSE

; :returns:
;    Returns the structure stx_drm containing the detector response matrix, pulse height matrix, energy loss matrix, STIX energy edges,
;    STIX mean energies of bins, energy bin widths and input parameters
;-

function stx_build_drm, ct_energy_edges, $
    ph_energy_edges = ph_energy_edges, $
    attenuator = attenuator, $
    d_al = d_al, $
    d_be = d_be, $
    func_par=func_par, $
    func_name = func, $
    efficiency = efficiency, $
    tailing = tailing, $
    background = background, $
    verbose=verbose, $
    _extra = _extra
    
  detector='cdte'     ; Cadmium-Telluride
  area = 1.0     ; detector geometric area in cm^2
  default, func, 'stx_fwhm'     ; returns FWHM in keV
  default, func_par, 1.0 ;only electronic component, Fano factor in stx_fwhm_fano and Hall tailing to be included separately
  default, attenuator, 0 ; no attenuator
  attenuator = fix( attenuator >0 < 1)
  ;func_par = [0.05, 120, fltarr(2)]  &  print, fwhm( [5.0, 120.], par)
  ;      1.31928      1.82805
  ;func_par = [0.08, 120, fltarr(2)]  & print, fwhm( [5.0, 120.], par)
  ;      2.11085      2.92488
  d=0.10       ; detector thickness in cm
  d0 = d
  
  default, background, 0
  
  if keyword_set(background) then begin
    d_be = 0
    d_al = 0
    d_pt = 0
    d = 100.
  endif
  
  
  default, d_be, 0.35     ; 3.5 mm Be
  default, d_al, 0.06     ; 600 microns (0.6 mm) of Al
  default, d_pt, 0.00015
  
  gmcm_be= 1.85 * d_be    ; 3.5 mm Be
  gmcm_al= 2.7  * d_al * attenuator
  gmcm_pt = 21.09* d_pt
  z = [13, 4, 78]  ; atomic number of elements in detector window
  default, gmcm, [gmcm_al, gmcm_be,gmcm_pt]
  ;thickness of each element in grams per cm2
  ; (gmcm = 0, if no window is assumed)
  
  
  
  edges_out = ct_energy_edges
  default, ph_energy_edges, ct_energy_edges
  default, edges_in, ph_energy_edges
  
  edge_products, ct_energy_edges,  mean = emean, width = win, edges_2= ect2
  if keyword_set( verbose ) then begin
    print, emean
    print, win
  endif
  edge_products, ph_energy_edges, edges_2 = eph2, mean = phmean, width = wout
  
  
  resp_calc, detector, area, func, func_par, d, z, gmcm, nflux, elo, ehi, $
    eloss_mat, pls_ht_mat, ein, smatrix, edges_in=edges_in, edges_out=edges_out, eff=efficiency, $
    _extra = _extra
    
  drm_info=[['DRM calculated with resp_calc_bg. smatrix is DRM (counts/keV/photons)'],$
    ['pls_ht_mat is DRM (counts/photons) NOT normalized, eloss_mat is Energy Loss Matrix, NO response broadening'],$
    ['ein...2D energy array, edges_in...1D energy array']]
    
  if keyword_set(background) then d = d0
  
  ;EIN    channel energy edges for photon input
  ;ELOSS_MAT  energy loss matrix
  ;PLS_HT_MAT pulse-height response matrix (ELOSS_MAT convolved with energy resolution broadening)
  ;SMATRIX  detector response matrix;
  ;                     = PLS_HT_MAT normalized to 1/keV
  
  
  
  if keyword_set( tailing ) then begin
    ;calculate the effect of hole taling as a matirx
    tailing_matrix =  stx_tailing_matrix(phmean, depth = d, include_damage = 1 )
    
    ;calculate the pulse height matrix and smatrix from the energly loss and tailing matrices
    eloss_mat = stx_tailing_products( eloss_mat, tailing_matrix, edges_in, edges_out, win, wout, func, func_par, area, $
      pls_ht_mat = pls_ht_mat_tail, smatrix = smatrix_tail )
      
    ;replace smatrix and pls_ht_mat with versions including tailing
    smatrix = smatrix_tail
    pls_ht_mat = pls_ht_mat_tail
    
  endif
  
  if keyword_set( verbose ) then help, smatrix
  
  stx_drm =      {detector:detector,$
    area:area,$
    func:func,$
    func_par:func_par,$
    d:d,$
    z:z,$
    gmcm:gmcm,$
    nflux:nflux,$
    eloss_mat:eloss_mat,$
    pls_ht_mat:pls_ht_mat, $
    e_2D:ein,$
    smatrix:smatrix,$ ;smatrix is DRM (counts/keV/photons)
    efficiency: efficiency, $ ;relative photopeak efficiency, multiply by area to get effective area for the peak
    edges_in:edges_in,$ ;1d
    edges_out:edges_out,$
    emean:emean,$ ;referenced to edges_out
    ewidth:win,$  ;referenced to edges_out
    background: background, $ ;if set then this response is used for background line calibration as there is no
    ;attenuation by windows and the detector is assumed to stop all photons and have
    ;no fluorescence
    info:drm_info}
  ;Compute the total efficiency for each input photon, almost, but not the same as diagonal eff
    
  return, stx_drm
  
end