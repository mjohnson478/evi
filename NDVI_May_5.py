import numpy as np
from binit import fastbin
import fasthist as fh
from matplotlib.colors import Normalize
from matplotlib import cm
from plot_rads import make_dir
import matplotlib.pyplot as plt
import pyhdf.SD
from mpl_toolkits.basemap import Basemap
import read_mod13
import os
import scipy.io as si

# #List of hdf files 2002
# modis_file_lst = ['MODIS_data_2/May9_2002/MOD13A1.A2002129.h10v04.005.2007150195436.hdf', 
# 'MODIS_data/16-Day_EVI_NDVI_May9_2002/MOD13A1.A2002129.h10v03.005.2007151130943.hdf',
# 'MODIS_data_2/May9_2002/MOD13A1.A2002129.h11v03.005.2007151004041.hdf', 
# 'MODIS_data_2/May9_2002/MOD13A1.A2002129.h11v04.005.2007150172839.hdf', 
# 'MODIS_data_2/May9_2002/MOD13A1.A2002129.h12v03.005.2007151005048.hdf',
# 'MODIS_data_2/May9_2002/MOD13A1.A2002129.h13v03.005.2007151131431.hdf']

#List of hdf files 2009
modis_file_lst = ['MODIS_data_2/May9_2009/MOD13A1.A2009129.h10v03.005.2009148223942.hdf',
'MODIS_data_2/May9_2009/MOD13A1.A2009129.h10v04.005.2009150010445.hdf',
'MODIS_data_2/May9_2009/MOD13A1.A2009129.h11v03.005.2009149221953.hdf',
'MODIS_data_2/May9_2009/MOD13A1.A2009129.h11v04.005.2009149193728.hdf',
'MODIS_data_2/May9_2009/MOD13A1.A2009129.h12v03.005.2009149073255.hdf',
'MODIS_data_2/May9_2009/MOD13A1.A2009129.h13v03.005.2009149025655.hdf']


#check to see if lat/lon grids exist and if not, create them
if os.path.exists('MODIS_data_2/May9_2009/lat_lon_grids.mat'):
    print 'lat lon grid file exists'
    data = si.loadmat('MODIS_data_2/May9_2009/lat_lon_grids.mat')
    lat_grids = data['lat_grids'].astype(np.float32)
    lon_grids = data['lon_grids'].astype(np.float32)

else:
    print ' No lat lon grid file exists, creating one'
    #loop to get lat/lon grids for each file
    lat_grids = []
    lon_grids = []
    for modis_file_name in modis_file_lst:
        data = read_mod13.lat_lon_data(modis_file_name)
        lat_grid = data['lat_grid'].astype(np.float32)
        lon_grid = data['lon_grid'].astype(np.float32)
        lat_grids.append(lat_grid)
        lon_grids.append(lon_grid)

    savedict = {'lat_grids': lat_grids, 'lon_grids': lon_grids}
    si.savemat('MODIS_data_2/May9_2009/lat_lon_grids.mat',savedict)


#print out corners of grids
for i in range(0,len(lat_grids)):

    lat_grid = lat_grids[i]
    lon_grid = lon_grids[i]

    #print corners of grid
    print (lon_grid[0,0], lat_grid[0,0])
    print (lon_grid[2399,0], lat_grid[2399,0])
    print (lon_grid[0,2399], lat_grid[0,2399])
    print (lon_grid[2399,2399], lat_grid[2399,2399])
    print ' '


#Set up NDVI grids
NDVI_grids = []
lat_centers_lst = []
lon_centers_lst = []
figcount = 1
for j in range(0,len(modis_file_lst)):
#for j in range(0,1):
    print 'Loop Index', j

    lat_grid = lat_grids[j]
    lon_grid = lon_grids[j]
    modis_file_name = modis_file_lst[j]

    numlatbins=1000
    numlonbins=1000

    #set limits for fastbin to grid corners
    north = lat_grid[0,0]
    south = lat_grid[2399,0]
    east = lon_grid[2399,2399]
    west = lon_grid[0,0]

    #lat/lon bins
    bin_lats=fastbin(south,north,numlatbins,-999,-888)
    bin_lons=fastbin(west,east,numlonbins,-999,-888)

    #lat/lon centers
    lon_centers=bin_lons.get_centers()
    lat_centers=bin_lats.get_centers()
    lon_centers_lst.append(lon_centers)
    lat_centers_lst.append(lat_centers)

    new_hist=fh.pyhist(lat_grid,lon_grid,bin_lats,bin_lons)
    lat_lon_counts=new_hist.get_hist2d()

    dirname='plots'
    make_dir(dirname)
    granule_info='count'

    #Figure 1: Lat/lon bin count
    fig1 = plt.figure(figcount)
    fig1.clf()
    cmap=cm.RdBu_r
    cmap.set_over('y')
    cmap.set_under('k')
    vmin= 0
    vmax= 20
    the_norm=Normalize(vmin=vmin,vmax=vmax,clip=False)
    axis1=fig1.add_subplot(111)
    im=axis1.pcolormesh(lon_centers,lat_centers,lat_lon_counts,cmap=cmap,norm=the_norm)
    cb=plt.colorbar(im,extend='both')
    the_label=cb.ax.set_ylabel('counts',rotation=270)
    axis1.set_title('{}: 2-d histogram (pixel count in each lat/lon bin'.format(granule_info))
    fig1.canvas.draw()
    #fig1.savefig('{0:s}/{1:s}_hist2d.png'.format(dirname,granule_info))
    figcount += 1

    #get NDVI data 
    model13_file = modis_file_name 
    NDVI=pyhdf.SD.SD(model13_file)
    NDVI_data=NDVI.select('500m 16 days NDVI')
    NDVI=NDVI_data.get()

    #scale NDVI data
    scale_NDVI=10000
    offset_NDVI=0
    NDVI = (NDVI*1.e-4).astype(np.float32)

    #set up NDVI grid
    NDVI_grid=new_hist.get_mean(NDVI)
    NDVI_grids.append(NDVI_grid)

    #Figure 2:
    fig1 = plt.figure(figcount)
    fig1.clf()
    del cmap
    cmap=cm.RdBu_r
    cmap.set_over('y')
    cmap.set_under('k')
    vmin= -1.
    vmax= 1.
    the_norm=Normalize(vmin=vmin,vmax=vmax,clip=False)

    axis1=fig1.add_subplot(111)
    im=axis1.pcolormesh(lon_centers,lat_centers,NDVI_grid,cmap=cmap,\
                      norm=the_norm)
    cb=plt.colorbar(im,extend='both')
    the_label=cb.ax.set_ylabel('NDVI',rotation=270)
    axis1.set_title('{}: NDVI'.format(granule_info))
    fig1.canvas.draw()
    #fig1.savefig('{0:s}/{1:s}_NDVI.png'.format(dirname,granule_info))
    figcount += 1


#Figure 3: NDVI plotted over base map.
fig = plt.figure(figcount)
fig.clf()
del cmap
cmap=cm.jet
cmap.set_over('y')
cmap.set_under('k', alpha = 0)
vmin= 0.
vmax= 1.
the_norm=Normalize(vmin=vmin,vmax=vmax,clip=False)
axis1=fig.add_subplot(111)

#create base map
lcc_transform = Basemap(llcrnrlon=-125,llcrnrlat=40,urcrnrlon=-85,urcrnrlat=60,projection='lcc',
            resolution='c',lat_1=40.,lat_2=60,lon_0=-100.)
#lcc_transform.fillcontinents(color='coral',lake_color='aqua')



for k in range(0,len(NDVI_grids)):
#for k in range(0,1):
    print 'Index', k
    NDVI_grid = NDVI_grids[k]
    lon_centers = lon_centers_lst[k]
    lat_centers = lat_centers_lst[k]

    #convert lat/lon to x,y coordinates
    lon_array,lat_array=np.meshgrid(lon_centers, lat_centers)
    x,y=lcc_transform(lon_array,lat_array)

    #plot NDVI on base map
    im=axis1.pcolormesh(x,y,NDVI_grid,cmap=cmap,\
                  norm=the_norm)


# draw coastlines.
lcc_transform.drawcoastlines(linewidth = 1)
lcc_transform.drawparallels(np.arange(40.,60.,1.),labels=[1,0,0,0],labelstyle='+/-',fontsize=8)
lcc_transform.drawmeridians(np.arange(-130.,-85.,2.),labels=[0,0,0,1],labelstyle='+/-',fontsize=8)
lcc_transform.drawstates(linewidth = 1)
lcc_transform.drawcountries(linewidth = 1)
cb=lcc_transform.colorbar(im,location='bottom',pad='10%',fig=fig,extend='both')
cb.set_ticks([0,0.25,0.5,0.75,1], update_ticks=True)
cb.set_label('NDVI')
axis1.set_title('NDVI May 2009')

#fig.savefig('{0:s}/{1:s}_NDVI_Mapped.png'.format(dirname,granule_info))
# si.savemat('MODIS_data_2/May9_2002/lat_lon_grids.mat',savedict)

fig.savefig('MODIS_data_2/May9_2009/NDVI_May_Figure2.png')

#plt.show()
