## Mapping storing default color maps (see R. Nelson notes)
# xco2      => viridis
# xco2_bias => RdBu_r (not directly used here)
# SIF       => YlGn
# TCWV      => Blues
# Radiance  => Purple_r
# TimeOfDay => YlGnBu_r (not directly used here)

default_cmaps = dict()

# L1 variables
default_cmaps['SoundingMeasurements/rad_continuum_o2'] = 'Purples_r'
default_cmaps['SoundingMeasurements/rad_continuum_weak_co2'] = 'Purples_r'
default_cmaps['SoundingMeasurements/rad_continuum_strong_co2'] = 'Purples_r'
default_cmaps['SoundingMeasurements/radiance_o2'] = 'Purples_r'
default_cmaps['SoundingMeasurements/radiance_weak_co2'] = 'Purples_r'
default_cmaps['SoundingMeasurements/radiance_strong_co2'] = 'Purples_r'

# LtCO2 variables
default_cmaps['xco2'] = 'viridis'
default_cmaps['Retrieval/xco2_raw'] = 'viridis'
default_cmaps['Preprocessors/xco2_weak_idp'] = 'viridis'
default_cmaps['Preprocessors/xco2_strong_idp'] = 'viridis'

default_cmaps['Retrieval/fs'] = 'YlGn'

default_cmaps['Retrieval/tcwv'] = 'Blues'

default_cmaps['Sounding/snr_o2a'] = 'Purples_r'
default_cmaps['Sounding/snr_wco2'] = 'Purples_r'
default_cmaps['Sounding/snr_sco2'] = 'Purples_r'

# LtSIF variables
default_cmaps['SIF_757nm'] = 'YlGn'
default_cmaps['SIF_771nm'] = 'YlGn'
default_cmaps['uncorrected_SIF_757nm'] = 'YlGn'
default_cmaps['uncorrected_SIF_771nm'] = 'YlGn'

default_cmaps['continuum_radiance_757nm'] = 'Purples_r'
default_cmaps['continuum_radiance_771nm'] = 'Purples_r'

