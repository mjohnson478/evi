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

#List of hdf files 2002
modis_file_lst = ['MODIS_data_2/June10_2002/MOD13A1.A2002161.h09v04.005.2008241032814.hdf',
'MODIS_data_2/June10_2002/MOD13A1.A2002161.h10v04.005.2008241023056.hdf',
'MODIS_data_2/June10_2002/MOD13A1.A2002161.h11v03.005.2008241021235.hdf',
'MODIS_data_2/June10_2002/MOD13A1.A2002161.h11v04.005.2008241020511.hdf',
'MODIS_data_2/June10_2002/MOD13A1.A2002161.h12v03.005.2008241025942.hdf',
'MODIS_data_2/June10_2002/MOD13A1.A2002161.h13v03.005.2008241062259.hdf',
'MODIS_data_2/June10_2002/MOD13A1.A2002161.h10v03.005.2008241065032.hdf']

#List of hdf files 2009
modis_file_lst = ['MODIS_data_2/June10_2009/MOD13A1.A2009161.h10v03.005.2009180094703.hdf',
'MODIS_data_2/June10_2009/MOD13A1.A2009161.h10v04.005.2009181023650.hdf',
'MODIS_data_2/June10_2009/MOD13A1.A2009161.h11v03.005.2009181002829.hdf',
'MODIS_data_2/June10_2009/MOD13A1.A2009161.h11v04.005.2009180211016.hdf',
'MODIS_data_2/June10_2009/MOD13A1.A2009161.h12v03.005.2009182053254.hdf',
'MODIS_data_2/June10_2009/MOD13A1.A2009161.h13v03.005.2009181170826.hdf']


#check to see if lat/lon grids exist and if not, create them
if os.path.exists('MODIS_data_2/June10_2009/lat_lon_grids.mat'):
    print 'lat lon grid file exists'
    data = si.loadmat('MODIS_data_2/June10_2009/lat_lon_grids.mat')
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
    si.savemat('MODIS_data_2/June10_2009/lat_lon_grids.mat',savedict)


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


#Set up EVI grids
EVI_grids = []
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

    #get EVI data 
    model13_file = modis_file_name 
    EVI=pyhdf.SD.SD(model13_file)
    evi_data=EVI.select('500m 16 days EVI')
    evi=evi_data.get()

    #scale EVI data
    scale_evi=10000
    offset_evi=0
    evi = (evi*1.e-4).astype(np.float32)

    #set up EVI grid
    evi_grid=new_hist.get_mean(evi)
    EVI_grids.append(evi_grid)

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
    im=axis1.pcolormesh(lon_centers,lat_centers,evi_grid,cmap=cmap,\
                      norm=the_norm)
    cb=plt.colorbar(im,extend='both')
    the_label=cb.ax.set_ylabel('EVI',rotation=270)
    axis1.set_title('{}: EVI'.format(granule_info))
    fig1.canvas.draw()
    #fig1.savefig('{0:s}/{1:s}_EVI.png'.format(dirname,granule_info))
    figcount += 1


#Figure 3: EVI plotted over base map.
fig = plt.figure(figcount)
fig.clf()
del cmap
cmap=cm.jet
cmap.set_over('y')
cmap.set_under('k', alpha = 0)
vmin= 0.
vmax= 0.75
the_norm=Normalize(vmin=vmin,vmax=vmax,clip=False)
axis1=fig.add_subplot(111)

#create base map
lcc_transform = Basemap(llcrnrlon=-125,llcrnrlat=40,urcrnrlon=-85,urcrnrlat=60,projection='lcc',
            resolution='c',lat_1=40.,lat_2=60,lon_0=-100.)
#lcc_transform.fillcontinents(color='coral',lake_color='aqua')



for k in range(0,len(EVI_grids)):
#for k in range(0,1):
    print 'Index', k
    evi_grid = EVI_grids[k]
    lon_centers = lon_centers_lst[k]
    lat_centers = lat_centers_lst[k]

    #convert lat/lon to x,y coordinates
    lon_array,lat_array=np.meshgrid(lon_centers, lat_centers)
    x,y=lcc_transform(lon_array,lat_array)

    #plot EVI on base map
    im=axis1.pcolormesh(x,y,evi_grid,cmap=cmap,\
                  norm=the_norm)


# draw coastlines.
lcc_transform.drawcoastlines(linewidth = 1)
lcc_transform.drawparallels(np.arange(40.,60.,1.),labels=[1,0,0,0],labelstyle='+/-',fontsize=8)
lcc_transform.drawmeridians(np.arange(-130.,-85.,2.),labels=[0,0,0,1],labelstyle='+/-',fontsize=8)
lcc_transform.drawstates(linewidth = 1)
lcc_transform.drawcountries(linewidth = 1)
cb=lcc_transform.colorbar(im,location='bottom',pad='10%',fig=fig,extend='both')
cb.set_ticks([0,0.25,0.5,0.75,1], update_ticks=True)
cb.set_label('EVI')
axis1.set_title('EVI June 2009')

#fig.savefig('{0:s}/{1:s}_EVI_Mapped.png'.format(dirname,granule_info))
# si.savemat('MODIS_data_2/MB_May9_2002/lat_lon_grids.mat',savedict)

fig.savefig('MODIS_data_2/June10_2009/EVI_June_Figure2.png')

#plt.show()