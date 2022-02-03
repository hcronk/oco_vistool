import h5py
import netCDF4
import sys
import numpy as np
import datetime

class LiteSIFFile:

    def __init__(self, lite_file):
        self.lite_file = lite_file

    def open_file(self):
        # many things were renamed in OCO3 B10... as an attempt to
        # simplify reading, determine if this is a B10 format file,
        # and set variable names accordingly.
        self.lf = netCDF4.Dataset(self.lite_file, 'r')
        if 'Metadata' in self.lf.groups:
            self._time_name = 'Delta_Time'
            self._mode_name = 'Metadata/MeasurementMode',
            self._orbit_name = 'Metadata/OrbitId'
            self._footprint_name = 'Metadata/FootprintId'
            self._sid_name = 'Metadata/SoundingId'
        else:
            self._time_name = 'time'
            self._mode_name = 'measurement_mode'
            self._orbit_name = 'orbit_number'
            self._footprint_name = 'footprint'
            self._sid_name = 'sounding_id'

    def get_product_version(self):
        # this does not work with the method setup in the open_file,
        # because the old way used an attribute.
        if 'Metadata' in self.lf.groups:
            return self.lf['Metadata/CollectionLabel'].getValue()
        return self.product_version

    # LtSIF stores the target name/ID very differently than
    # the LtCO2 - in each case, derive an array that matches
    # how it is done in LtCO2.
    def get_target_id(self):
        num_sounding = self.lf.dimensions['sounding_dim'].size
        target_id = np.zeros(num_sounding, dtype=object)
        target_id[:] = 'none'
        if 'Sequences' in self.lf.groups:
            for s in range(self.lf.dimensions['sequences_dim'].size):
                msk = self.lf['Sequences/SequencesIndex'][:] == s
                target_id[msk] = self.lf['Sequences/SequencesId'][s]
        return target_id

    def get_target_name(self):
        num_sounding = self.lf.dimensions['sounding_dim'].size
        target_name = np.zeros(num_sounding, dtype=object)
        target_name[:] = 'none'
        if 'Sequences' in self.lf.groups:
            for s in range(self.lf.dimensions['sequences_dim'].size):
                msk = self.lf['Sequences/SequencesIndex'][:] == s
                target_name[msk] = self.lf['Sequences/SequencesName'][s]
        return target_name

    # note the Lite SIF has a different code -> string mapping for
    # the operation mode.
    def get_operation_mode(self):
        opmode_code = self.lf[self._mode_name][:]
        opmode = np.zeros(opmode_code.shape[0], dtype='S3')
        opmode_map = ['ND', 'GL', 'TG', 'SAM', 'XS']
        for code, opmode_str in enumerate(opmode_map):
            msk = opmode_code == code
            opmode[msk] = opmode_str
        return opmode

    def get_sid(self):
        return self.lf[self._sid_name][:]
        
    def get_orbit(self):
        return self.lf[self._orbit_name][:]
        
    def get_footprint(self):
        return self.lf[self._footprint_name][:]

    def get_time(self):

        # gets time and shifts it to match the 1970 reference time
        # (standard unix time stamp)
        # this matches LtCO2 except that there is a leap second
        # mismatch (the same sounding ID in LtCO2 will be off by 
        # the leap second change.)

        t = self.lf[self._time_name][:]

        # B10 SIF has a different time reference
        if 'Delta_Time' in nc.lf.variables:
            dt = datetime.datetime(1990,1,1) - datetime.datetime(1970,1,1)
        else:
            dt = datetime.datetime(1993,1,1) - datetime.datetime(1970,1,1)

        t = t + dt.total_seconds()

        return t

    def close_file(self):
        self.lf.close()


class LiteCO2File:

    def __init__(self, lite_file):
        self.lite_file = lite_file
        
    def open_file(self):
        self.lf = netCDF4.Dataset(self.lite_file, 'r')

    def get_product_version(self):
        # try several different variables... there is not consistency
        # among versions here.
        try:
            product_version = self.lf.getncattr('CollectionLabel')
        except AttributeError:
            try:
                product_version = self.lf.getncattr('BuildId')
            except:
                product_version = 'unknown'
        return product_version

    def get_target_id(self):
        if 'target_id' in self.lf['Sounding'].variables:
            target_id = self.lf['Sounding/target_id'][:]
        else:
            target_id = np.zeros(self.lf.dimensions['sounding_id'].size, object)
            target_id[:] = 'none'
        return target_id

    def get_target_name(self):
        if 'target_name' in self.lf['Sounding'].variables:
            target_name = self.lf['Sounding/target_name'][:]
        else:
            target_name = np.zeros(self.lf.dimensions['sounding_id'].size, object)
            target_name[:] = 'none'
        return target_name

    def get_operation_mode(self):
        opmode_code = self.lf['/Sounding/operation_mode'][:]
        opmode = np.zeros(opmode_code.shape[0], dtype=object)
        opmode_map = ['ND', 'GL', 'TG', 'XS', 'SAM']
        for code, opmode_str in enumerate(opmode_map):
            msk = opmode_code == code
            opmode[msk] = opmode_str
        return opmode

    def get_warn(self):
        return self.lf['warn_level'][:]

    def get_qf(self):
        return self.lf['xco2_quality_flag'][:]
    
    def get_sid(self):
        return self.lf['sounding_id'][:]
        
    def get_orbit(self):
        return self.lf['/Sounding/orbit'][:]
        
    def get_footprint(self):
        return self.lf['/Sounding/footprint'][:]

    def get_time(self):
        return self.lf['time'][:]
    
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
        
    def get_total_aod(self):
        return self.l2f['/AerosolResults/aerosol_total_aod'][:]
        
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
        
    def get_sfc_type(self):
        return self.l2f['/RetrievalResults/surface_type'][:]
    
    def get_solar_zenith(self):
        return self.l2f['/RetrievalGeometry/retrieval_solar_zenith'][:]
    
    def get_sensor_zenith(self):
        return self.l2f['/RetrievalGeometry/retrieval_zenith'][:]
        
    def get_lat(self):
        return self.l2f['/RetrievalGeometry/retrieval_latitude'][:]
    
    def get_lon(self):
        return self.l2f['/RetrievalGeometry/retrieval_longitude'][:]
        
    def get_sfc_rough(self):
        return self.l2f['/RetrievalGeometry/retrieval_surface_roughness'][:]
    
    def get_rel_resid_mean_sq_strongCO2(self):
        return self.l2f['/SpectralParameters/relative_residual_mean_square_strong_co2'][:]
    
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
        
    def get_modis_solar_zenith(self):
        solar_zenith_obj = self.sf['SolarZenith']
        solar_zenith = solar_zenith_obj[:]
        solar_zenith_scale = solar_zenith_obj.attrs.get('scale_factor')[0]
        try: 
            solar_zenith_offset = solar_zenith_obj.attrs.get('add_offset')[0]
        except: 
            solar_zenith_offset = 0
        solar_zenith = solar_zenith.astype(float)
        solar_zenith *= solar_zenith_scale
        solar_zenith += solar_zenith_offset
        
        return solar_zenith
        
    def get_oco2_file(self):
        return self.sf['OCO_Input_File'][:][0]
    
    def get_oco2_lat(self):
        return self.sf['OCO_Latitude'][:]
        
    def get_oco2_lon(self):
        return self.sf['OCO_Longitude'][:]
    
    def get_xmatch(self):
        self.xmatch = self.sf['matchup_Xindex'][:]
        return self.xmatch
        
    def get_ymatch(self):
        #self.ymatch = self.sf['matchup_Yindex'][:]
        return self.sf['matchup_Yindex'][:]
            
    def get_match_dist(self):
        return self.sf['matchup_distance_km'][:]
    
    def make_matched_array(self, var):
        xmatch = self.get_xmatch()
        ymatch = self.get_ymatch()
        dmatch = self.get_match_dist()
        match_mask = np.logical_or(np.logical_or(xmatch == -999, ymatch == -999), dmatch > 5)
        xmatch = np.ma.array(xmatch, mask = match_mask)
        ymatch = np.ma.array(ymatch, mask = match_mask)
        
        n_oco_soundings = xmatch.size
        
        xmatch_1D = np.reshape(xmatch, n_oco_soundings)
        ymatch_1D = np.reshape(ymatch, n_oco_soundings)
        
        matched_var = np.full_like(xmatch_1D, -999, dtype = var.dtype)
        
        for i, x_idx in enumerate(xmatch_1D):
            if x_idx:
                matched_var[i] = var[x_idx, ymatch_1D[i]]
        
        return matched_var.reshape(xmatch.shape)
    
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

class GeoScFile:

    def __init__(self, geo_file):
        self.geo_file = geo_file

    def open_file(self):
        self.gsf = h5py.File(self.geo_file, 'r')
    
    def get_lat(self):
        return self.gsf['/SoundingGeometry/sounding_latitude']
    
    def get_lon(self):
        return self.gsf['/SoundingGeometry/sounding_longitude']
    
    def get_start_time(self):
        self.id = self.gsf['/SoundingGeometry/sounding_id'][:]
        start = (self.id[0])[0]
        count=1
        while True:
            if start <= 0:
                start = (self.id[count])[0]
                count+=1
            else:
                break

        return(start)
        
    def get_end_time(self):
        self.id = self.gsf['/SoundingGeometry/sounding_id'][:]
        end = (self.id[-1])[-1]
        count=2
        while True:
            if end <= 0:
                end = (self.id[-count])[-1]
                count+=1
            else:
                break
        
        return(end)
    
    def get_collection(self):
        return(self.gsf['/Metadata/CollectionLabel'].value[0])
    
    def get_vertex_lat(self):
        return self.gsf['/FootprintGeometry/footprint_vertex_latitude'][:,:,0,:]
        
    def get_vertex_lon(self):
        return self.gsf['/FootprintGeometry/footprint_vertex_longitude'][:,:,0,:]
    
    def close_file(self):
        self.gsf.close()
