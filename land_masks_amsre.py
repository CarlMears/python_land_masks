
class AMSRE_Landmask:

    def __init__(self):

        import numpy as np
        from  scipy.interpolate import RectBivariateSpline

        percent_land_file = 'L:/sea_ice_ancillary_for_Ball/static_ancillary/land_masks/land_tables_amsra.dat'
    
        num_lats = 2160
        num_lons = 4320
    
        dy = 180.0/num_lats
        dx = 360.0/num_lons

        lm = np.fromfile(percent_land_file,dtype = 'uint8')

        assert(lm.shape[0] == 5*num_lats*num_lons)

        land_mask_temp = np.reshape(lm,(5,num_lats,num_lons))
        land_mask_temp = 0.4*land_mask_temp
        land_mask = np.zeros((5,num_lats+2,num_lons+2))
        for i in range(0,5):
            land_mask[i,1:num_lats+1,1:num_lons+1] = land_mask_temp[i,:,:]
            land_mask[i,0,:] = 100.0
            land_mask[i,num_lats+1,:] = 0.0
            land_mask[i,:,0] = land_mask[i,:,num_lons]
            land_mask[i,:,num_lons+1] = land_mask[i,:,1]
  
        lons = np.linspace(0.0 - dx/2.0,360.0 + dx/2.0,num = num_lons+2)
        lats = np.linspace(-90.0  - dy/2.0,90.0  + dy/2.0,num = num_lats+2)

        self.land_mask= land_mask
        self.lons = lons
        self.lats = lats
        self.dlon = dx
        self.dlat = dy

        # this sets up the interpolator to use simple linear interpolation
        # kx and ky are the degree of the spline -- 1 means linear
        # using land mask 2 means we are using the 19/23 GHz amsre footprint size

        self.interp_vlo_res = RectBivariateSpline(self.lats,self.lons,self.land_mask[0,:,:],kx=1,ky=1)
        self.interp_low_res = RectBivariateSpline(self.lats,self.lons,self.land_mask[1,:,:],kx=1,ky=1)
        self.interp_md_res =  RectBivariateSpline(self.lats,self.lons,self.land_mask[2,:,:],kx=1,ky=1)
        self.interp_hi_res =  RectBivariateSpline(self.lats,self.lons,self.land_mask[3,:,:],kx=1,ky=1)

    def percent_land_interp(self,lat=None,lon=None,res = 'md'):

        #vlo corresponds to AMSR2 7 GHz footprints
        #low corresponds to AMSR2 11 GHz footprints
        #md corresponds to AMSR2 19 GHz footprints
        #hi corresponds to AMSR2 37 GHz footprints
        

        if res == 'md':
            percent_land = self.interp_md_res.ev(lat,lon)
        elif ((res == 'low') or (res == 'lo')):
            percent_land = self.interp_low_res.ev(lat,lon)
        elif res == 'vlo':
            percent_land = self.interp_vlo_res.ev(lat,lon)
        elif res == 'hi':
            percent_land = self.interp_hi_res.ev(lat,lon)
        else:
            raise ValueError(f"Resolution: {res} not implemented")
        return percent_land


# test code
if __name__ == '__main__':

    import sys
    sys.path.append('//mears/C/job_CCMP/python/plotting_and_analysis/')
    import numpy as np
    from global_map import global_map
    import matplotlib.pyplot as plt

    land_percent = AMSRE_Landmask()
  
    lons = np.linspace(0.5, 359.5, num=360)
    lats = np.linspace(-89.5, 89.5, num=180)

    #much faster to call the interpolator all at once using np arrays
    lonsv,latsv = np.meshgrid(lons,lats)
    test_map = land_percent.percent_land_interp(lat=latsv, lon=lonsv,res='vlo')

    fig,ax = global_map(test_map, vmin=0.0, vmax=100.0)

    test_map = land_percent.percent_land_interp(lat=latsv, lon=lonsv,res='hi')

    fig,ax = global_map(test_map, vmin=0.0, vmax=100.0)
    plt.show()

'''
# global map, in case you can't access //mears/c/job_ccmp

def global_map(a, vmin=0.0, vmax=30.0, cmap=None, plt_colorbar=False,title='',extent=None,central_longitude = 0.0):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import cartopy.crs as ccrs

    img_extent = [-180.0, 180.0, -90.0, 90.0]
    fig = plt.figure(figsize=(10, 5))  # type: Figure
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude = central_longitude),title=title)
    for item in ([ax.title]):
        item.set_fontsize(16)
    sz = a.shape
    num_lons = sz[1]

    map = ax.imshow(np.flipud(np.roll(a, int(num_lons/2), axis=1)), cmap=cmap, origin='upper', transform=ccrs.PlateCarree(),
                    norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax), extent=img_extent)
    if plt_colorbar:
        cbar = fig.colorbar(map, shrink=0.7, orientation='horizontal')
        cbar.ax.tick_params(labelsize=14)
    ax.coastlines()
    ax.set_global()
    if extent is not None:
        ax.set_extent(extent,crs = ccrs.PlateCarree())
    return fig, ax
'''

'''
Here is the original Fortran Module:
Fortran does not interpolate, so results should not be exactly the same....

module percent_land
    
    use SeaIceCtypes,    only: c_int, c_float, c_double, c_int32_t, c_int8_t
    use SeaIceFilePaths, only: land_percent_file
    use logger,          only : set_filter,log_debug,log_info,log_error,log_crit,log_msg_string
    
    character(1)       :: percent_land_tab(4320,2160,5)
    integer(c_int32_t) :: table_loaded = 0

    contains
    
    subroutine init_percent_land(error)
        integer(c_int32_t),intent(out)  :: error
        
        logical(4)  :: file_exists
        
        inquire(file=land_percent_file,exist=file_exists)
        if (file_exists) then
            open(unit=3, file=land_percent_file, status='old', access='stream', form='unformatted')
            read(3) percent_land_tab
            close(3)
            table_loaded = 1
            error = 0
        else
            call log_error('Error Opening Land Percentage File')
            call log_error('Expected Location: ')
            call log_error(land_percent_file)
            table_loaded = 0
            error = -1
        endif
    end subroutine init_percent_land
    
    subroutine fd_percent_land(xlat,xlon,percent_land)
        
        implicit none
        real(c_float),intent(in)     :: xlat
        real(c_float),intent(in)     :: xlon
        real(c_float),intent(out)    :: percent_land(5)
        
        integer(c_int32_t)           :: ilat
        integer(c_int32_t)           :: ilon
        integer(c_int32_t)           :: error
        
        if(table_loaded .ne. 1) then
            call init_percent_land(error)
            if (error .ne. 0) then
                call log_error('Fatal Error')
                stop
            endif
        endif
        
        ilat=1+nint(12.*(xlat+89.95833))
        ilon=1+nint(12.*(xlon- 0.04167))
        if(ilat.lt.1) ilat=1
        if(ilon.lt.1) ilon=1
        if(ilat.gt.2160) ilat=2160
        if(ilon.gt.4320) ilon=4320
        
        percent_land=0.4*ichar(percent_land_tab(ilon,ilat,:))
        
        return
    end subroutine fd_percent_land
end module percent_land
'''