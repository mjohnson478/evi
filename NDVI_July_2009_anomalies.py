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
import glob
from scipy import stats



# modis_file_lst1 = ['MODIS_data_3/2000/MOD13A1.A2000193.h10v04.005.2008274154726_1.hdf',
# 'MODIS_data_3/2001/MOD13A1.A2001193.h10v04.005.2007039120415_1.hdf',
# 'MODIS_data_3/2002/MOD13A1.A2002193.h10v04.005.2007177192540_1.hdf',
# 'MODIS_data_3/2003/MOD13A1.A2003193.h10v04.005.2008314110654_1.hdf',
# 'MODIS_data_3/2004/MOD13A1.A2004193.h10v04.005.2007298211504_1.hdf',
# 'MODIS_data_3/2005/MOD13A1.A2005193.h10v04.005.2008215134939_1.hdf',
# 'MODIS_data_3/2006/MOD13A1.A2006193.h10v04.005.2008137022300_1.hdf',
# 'MODIS_data_3/2007/MOD13A1.A2007193.h10v04.005.2007212170755_1.hdf',
# 'MODIS_data_3/2008/MOD13A1.A2008193.h10v04.005.2008210214211_1.hdf',
# 'MODIS_data_3/2009/MOD13A1.A2009193.h10v04.005.2009212091120_1.hdf']

modis_dic = {}
for x in range(0,6):
    y = x+1
    os.chdir('MODIS_data_3/' + str(y) + '/')
    modis_file_lst = glob.glob('*.hdf')

    for z in range(0,len(modis_file_lst)):
        modis_file_lst[z] = 'MODIS_data_3/' + str(y) + '/' + modis_file_lst[z] 

    modis_dic[str(y)] = modis_file_lst
    os.chdir('..')
    os.chdir('..')



#check to see if lat/lon grids exist and if not, create them
if os.path.exists('MODIS_data_3/lat_lon_grids.mat'):
    print 'lat lon grid file exists'
    data = si.loadmat('MODIS_data_3/lat_lon_grids.mat')
    lat_grids = data['lat_grids'].astype(np.float32)
    lon_grids = data['lon_grids'].astype(np.float32)

else:
    print ' No lat lon grid file exists, creating one'
    #loop to get lat/lon grids for each file
    lat_grids = []
    lon_grids = []

    for i in range(0, len(modis_dic)):
        j = i+1
        modis_file_name = modis_dic[str(j)][0]

        data = read_mod13.lat_lon_data(modis_file_name)
        lat_grid = data['lat_grid'].astype(np.float32)
        lon_grid = data['lon_grid'].astype(np.float32)
        lat_grids.append(lat_grid)
        lon_grids.append(lon_grid)

    savedict = {'lat_grids': lat_grids, 'lon_grids': lon_grids}
    si.savemat('MODIS_data_3/lat_lon_grids.mat',savedict)



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
NDVI_2009_anomalies_grids = []
lat_centers_lst = []
lon_centers_lst = []
figcount = 1
#loop over 6 images
for j in range(0,len(modis_dic)):
#for j in range(1,2):
#for j in range(0,1):
    print 'Loop Index', j

    lat_grid = lat_grids[j]
    lon_grid = lon_grids[j]
    #modis_file_name = modis_file_lst[j]

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

    #loop over years
    modis_file_names_lst = modis_dic[str(j+1)]
    NDVI_grids = []
    for a in range(0,len(modis_file_names_lst)):
        modis_file_name = modis_file_names_lst[a]

        #get NDVI data 
        model13_file = modis_file_name 
        NDVI=pyhdf.SD.SD(model13_file)
        NDVI_data=NDVI.select('500m 16 days NDVI')
        NDVI_int=NDVI_data.get()
        #print sum(NDVI_int.ravel() == -3000)
        NDVI = NDVI_int.astype('float32')
        

        # if j==2:
        #     plt.figure(figcount)
        #     figcount += 1
        #     plt.hist(NDVI, bins = 10000)
        #     s = 'MODIS_data_3/image2_' + str(a) + ".png"
        #     fig.savefig(s)

        hit = NDVI == -3000
        #print sum(hit.ravel())
        NDVI[hit] = np.nan
        

        #scale NDVI data
        scale_NDVI=10000
        offset_NDVI=0
        NDVI = (NDVI*1.e-4).astype(np.float32)
        #print NDVI.max(), NDVI.min()

        #set up NDVI grid
        NDVI_grid=new_hist.get_mean(NDVI)
        hit2 = NDVI_grid == -999
        NDVI_grid[hit2] = np.nan
        NDVI_grids.append(NDVI_grid)

        if j ==1:
            NDVI_grid_dct = {}
            NDVI_grid_dct[str(a)] = NDVI_grid
            s2 = 'MODIS_data_3/NDVI_grid_image2_' + str(a) + '.mat'
            si.savemat(s2,NDVI_grid_dct)

        #for grid in NDVI_grids:
            #print grid.shape

    NDVI_cube = np.empty([1000,1000,10], dtype='float32')
    for g in range(0,len(NDVI_grids)):
        NDVI_cube[:,:,g] = NDVI_grids[g]

    NDVI_avg_grid = stats.nanmean(NDVI_cube, axis=2)

    for c in range(0,len(modis_file_names_lst)):
        if modis_file_names_lst[c][24:28] == '2009':
            index = c

    NDVI_2009_anomalies = NDVI_grids[index] - NDVI_avg_grid

    # if j ==1:
    #         NDVI_anomalies_dct = {}
    #         NDVI_anomalies_dct[str(a)] = NDVI_2009_anomalies
    #         s3 = 'MODIS_data_3/NDVI_anomalies_image2_' + str(a) + '.mat'
    #         si.savemat(s3,NDVI_anomalies_dct)


        # NDVI_sum_grid = []
        # for b in range(0,len(NDVI_grids)):
        #     print NDVI_grids[b].max(), NDVI_grids[b].min()
        #     if b==0:
        #         NDVI_sum_grid = NDVI_grids[b]
        #     else:
        #         NDVI_sum_grid += NDVI_grids[b]

        # NDVI_avg_grid = NDVI_sum_grid/len(NDVI_grids)

    NDVI_2009_anomalies_grids.append(NDVI_2009_anomalies)

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
vmin= -0.3
vmax= 0.1
the_norm=Normalize(vmin=vmin,vmax=vmax,clip=False)
axis1=fig.add_subplot(111)

#create base map
lcc_transform = Basemap(llcrnrlon=-125,llcrnrlat=40,urcrnrlon=-85,urcrnrlat=60,projection='lcc',
            resolution='c',lat_1=40.,lat_2=60,lon_0=-100.)
#lcc_transform.fillcontinents(color='coral',lake_color='aqua')



for k in range(0,len(NDVI_2009_anomalies_grids)):
#for k in range(0,1):
    print 'Index', k
    NDVI_grid = NDVI_2009_anomalies_grids[k]
    #print NDVI_avg_grids[k].max(), NDVI_avg_grids[k].min()
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
cb.set_ticks([-0.3,-0.2,-0.1,0,0.1], update_ticks=True)
cb.set_label('NDVI')
axis1.set_title('NDVI Anomalies 2009')

#fig.savefig('{0:s}/{1:s}_NDVI_Mapped.png'.format(dirname,granule_info))
# si.savemat('MODIS_data_2/MB_May9_2002/lat_lon_grids.mat',savedict)

fig.savefig('MODIS_data_3/NDVI_July_2009_anomalies_Figure.png')

#plt.show()