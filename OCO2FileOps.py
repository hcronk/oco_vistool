import h5py
import sys
import numpy as np

class LiteFile:

    def __init__(self, lite_file):
        self.lite_file = lite_file

    def open_file(self):
        self.lf = h5py.File(self.lite_file, 'r')
    
    def get_lat(self):
        return self.lf['latitude'][:]
	
    def get_lon(self):
        return self.lf['longitude'][:]

    def get_xco2(self):
        return self.lf['xco2'][:]
	
    def get_warn(self):
        return self.lf['warn_level'][:]

    def get_qf(self):
        return self.lf['xco2_quality_flag'][:]
    
    def get_sid(self):
	return self.lf['sounding_id'][:]
	
    def get_orbit(self):
        return self.lf['/Sounding/orbit'][:]

    def get_sfcP(self):
        return self.lf['/Retrieval/psurf'][:]
	
    def get_sfcP_apriori(self):
        return self.lf['/Retrieval/psurf_apriori'][:]
    
    def get_delta_sfcP(self):
        return self.lf['/Preprocessors/dp_abp'][:]

    def get_co2_ratio(self):
        return self.lf['/Preprocessors/co2_ratio'][:]
    
    def close_file(self):
        self.lf.close()

class L2File:

    def __init__(self, l2_file):
        self.l2_file = l2_file

    def open_file(self):
        self.l2f = h5py.File(self.l2_file, 'r')
    
    def get_delta_sfcP(self):
        return self.l2f['/PreprocessingResults/surface_pressure_delta_abp'][:]
	
    def get_abp_cloud_flag(self):
        return self.l2f['/PreprocessingResults/cloud_flag_abp'][:]
    
    def get_co2_ratio_idp(self):
        return self.l2f['/PreprocessingResults/co2_ratio_idp'][:]
	
    def get_h2o_ratio_idp(self):
        return self.l2f['/PreprocessingResults/h2o_ratio_idp'][:]
    
    def get_l2_sid(self):
        return self.l2f['/RetrievalHeader/sounding_id'][:]

    def get_l1b_sid(self):
        return self.l2f['/L1bScSoundingReference/sounding_id_l1b'][:]
    
    def get_outcome_flag(self):
        return self.l2f['/RetrievalResults/outcome_flag'][:]
    
    def get_sounding_QF(self):
    	return self.l2f['/L1bScSoundingReference/sounding_qual_flag'][:]
    
    def get_warn_level(self):
        """
	This is a filled legacy field. For real warn level, use the lite file.
	"""
        return self.l2f['/RetrievalHeader/warn_level'][:]
	
    def get_packaging_QF(self):
        return self.l2f['/L1bScSoundingReference/packaging_qual_flag'][:]
	
    def get_aerosol_types(self):
        return self.l2f['/AerosolResults/aerosol_types'][:]
	
    def get_aerosol_1_aod(self):
        return self.l2f['/AerosolResults/aerosol_1_aod'][:]

    def get_aerosol_2_aod(self):
        return self.l2f['/AerosolResults/aerosol_2_aod'][:]

    def get_aerosol_3_aod(self):
        return self.l2f['/AerosolResults/aerosol_3_aod'][:]
	
    def get_aerosol_4_aod(self):
        return self.l2f['/AerosolResults/aerosol_4_aod'][:]
	
    def get_albedo_o2_fph(self):
        return self.l2f['/AlbedoResults/albedo_o2_fph'][:]
	
    def get_albedo_strong_co2_fph(self):
        return self.l2f['/AlbedoResults/albedo_strong_co2_fph'][:]
	
    def get_albedo_weak_co2_fph(self):
        return self.l2f['/AlbedoResults/albedo_weak_co2_fph'][:]
	
    def get_sfcP_fph(self):
        return self.l2f['/RetrievalResults/surface_pressure_fph'][:]
	
    def get_sfcP_apriori_fph(self):
        return self.l2f['/RetrievalResults/surface_pressure_apriori_fph'][:]
	
    def get_dof_co2_profile(self):
        return self.l2f['/RetrievalResults/dof_co2_profile'][:]

    def get_co2_vertical_gradient_delta(self):
        return self.l2f['/RetrievalResults/co2_vertical_gradient_delta'][:]
    
    def close_file(self):
        self.l2f.close()
	
class PreprocFile:

    def __init__(self, pp_file):
        self.pp_file = pp_file

    def open_file(self):
        self.ppf = h5py.File(self.pp_file, 'r')
    
    def get_delta_sfcP(self):
        return self.ppf['/ABandRetrieval/surface_pressure_delta_abp'][:]
	
    def get_abp_cloud_flag(self):
        return self.ppf['/ABandRetrieval/cloud_flag_abp'][:]
    
    def get_co2_ratio(self):
        return self.ppf['/DOASCloudScreen/co2_ratio_idp'][:]
    
    def get_sid(self):
        return self.ppf['/SoundingGeometry/sounding_id'][:]
	
    def get_albedo(self):
        return self.ppf['/ABandRetrieval/albedo_o2_abp'][:, :, 0]
    
    def close_file(self):
        self.ppf.close()
	   

class SubsetFile:

    def __init__(self, subset_file):
        self.subset_file = subset_file

    def open_file(self):
        self.sf = h5py.File(self.subset_file, 'r')
    
    def get_modis_lat(self):
        return self.sf['MODIS_Latitude'][:]
	
    def get_modis_lon(self):
        return self.sf['MODIS_Longitude'][:]

    def get_oco2_lat(self):
        return self.sf['OCO_Latitude'][:]
	
    def get_oco2_lon(self):
        return self.sf['OCO_Longitude'][:]
    
    def get_xmatch(self):
        return self.sf['matchup_Xindex'][:]
	
    def get_ymatch(self):
    	#try:
        ymatch = self.sf['matchup_Yindex'][:]
	return ymatch
	#except ValueError:
	    #print
	    
    def get_match_dist(self):
        return self.sf['matchup_distance_km'][:]
    
    def get_sid(self):
	return self.sf['OCO_sounding_id'][:]
	
    def get_1km_overall_cloud_flag(self):
        modis_cloud_flags = self.sf['Cloud_Mask_1km'][:, :, 0]
	
	#Get bits 1-2 for overall MODIS cloud flag
        modis_cloud_flag = np.right_shift(np.bitwise_and(modis_cloud_flags, 6), 1)
	#modis_cloud_flag = np.bitwise_and(modis_cloud_flags, 6)
	
	return modis_cloud_flag
	
    def get_1km_cirrus_refl(self):
        modis_cirrus_refl_obj = self.sf['Cirrus_Reflectance']
        modis_cirrus_refl = modis_cirrus_refl_obj[:]
        modis_cirrus_refl_scale = modis_cirrus_refl_obj.attrs.get('scale_factor')[0]
	modis_cirrus_refl_offset = modis_cirrus_refl_obj.attrs.get('add_offset')[0]
	
	#Scale the cirrus reflectance
        modis_cirrus_refl = modis_cirrus_refl.astype(float)
        modis_cirrus_refl = np.ma.masked_where(modis_cirrus_refl == -9999.0, modis_cirrus_refl)
        modis_cirrus_refl *= modis_cirrus_refl_scale
	modis_cirrus_refl += modis_cirrus_refl_offset
	
	return modis_cirrus_refl
	
    def get_calipso_lat(self):
        return self.sf['CALIPSO_Latitude'][:]
	
    def get_calipso_lon(self):
        return self.sf['CALIPSO_Longitude'][:]
	
    def get_modis_rgb(self):
        """
	Red = Band 1
	Green = Band 4
	Blue = Band 3
	"""
        refl_250_obj = self.sf['EV_250_Aggr500_RefSB']
	refl_250_scale = refl_250_obj.attrs.get('reflectance_scales')
	refl_250_offset = refl_250_obj.attrs.get('reflectance_offsets')
	refl_250_fill = refl_250_obj.attrs.get('_FillValue')
	refl_500_obj = self.sf['EV_500_RefSB']
	refl_500_scale = refl_500_obj.attrs.get('reflectance_scales')
	refl_500_offset = refl_500_obj.attrs.get('reflectance_offsets')
	refl_500_fill = refl_500_obj.attrs.get('_FillValue')
	
	b1 = refl_250_obj[:, :, 0]
	b1 = np.ma.masked_where(b1 == refl_250_fill, b1)
	b1 = b1.astype(float)
	b1 = b1 * refl_250_scale[0] + refl_250_offset[0]
	
	b4 = refl_500_obj[:, :, 1]
	b4 = np.ma.masked_where(b4 == refl_500_fill, b4)
	b4 = b4.astype(float)
	b4 = b4 * refl_500_scale[1] + refl_500_offset[1]
	
	b3 = refl_500_obj[:, :, 0]
	b3 = np.ma.masked_where(b3 == refl_500_fill, b3)
	b3 = b3.astype(float)
	b3 = b3 * refl_500_scale[0] + refl_500_offset[0]
	
	return b1, b4, b3
    
    def close_file(self):
        self.sf.close()
