# modified DP 25 June 2023
# Modified Ash - 17th Nov 2023 


# Imports
execfile("get_simobs_mods.py")
import os 
import numpy as np

# Define
min_baseline =  [14.6, 14.6, 14.6, 14.6, 14.6, 14.6, 64.0, 110.4, 367.6, 244.0]
mrs =           [28.5, 22.6, 16.2, 11.2, 6.7, 4.11, 2.58, 1.42, 0.814, 0.496]
beam_size =     [3.38, 2.3, 1.42, 0.918, 0.545, 0.306, 0.211, 0.096, 0.057, 0.042]
configuration = ['conf1', 'conf2', 'conf3', 'conf4', 'conf5', 'conf6', 'conf7', 'conf8', 'conf9', 'conf10']


# Run simulate and clean
min_baseline = np.array(min_baseline)
mrs = np.array(mrs)
beam_size = np.array(beam_size)

pix = beam_size/5.
imagesize=420 #size in pix of conf1=2.5'
image_size=np.zeros(10)
base_source_size = []

for i in range (0,10):

    base_source_size.append(1.5*mrs[0]*min_baseline[0]/min_baseline[i]/15)
    image_size[i] = int(imagesize*pix[0]/pix[i])

    if image_size[i] % 2 == 0:
        image_size[i] = image_size[i]+0 
    else:
        image_size[i] = image_size[i]-1

for i in range (0,10):
    for j in range (0,15):

        run_only = 'conf5'
        if configuration[i] != run_only: 
            print('[INFO] Only running for %s - skipping...' %run_only)
            continue 

        # if j!=6:
        #     continue

        source_size = base_source_size[i]*(j+1)
        
        myskymodel_dir = './skymodels/'
        myproject_dir = './'

        for totaltime in ['1', '6', '30', '60']:

            myproject = myproject_dir+configuration[i]+'_'+str('%.1f'%(source_size))+'mrs0_project_%stotaltime' %totaltime
            myskymodel = myskymodel_dir+configuration[i]+'_'+'ref_image'+str('%.1f'%(source_size))+'mrs0'

            os.system('rm -rf %s' %myproject)
            os.system('mkdir %s' %myproject)

            print('[INFO] myproject %s' %myproject)
            print('[INFO] myskymodel %s' %myskymodel)
            
            print('-----------------')
            print('Simulating ' + myproject)
            print('-----------------')
            print('Starting simobserve')

            simobserve(project=myproject, 
                        skymodel=myskymodel, 
                        obsmode="int", 
                        antennalist='alma.cycle5.'+str(i+1)+'.cfg', 
                        thermalnoise='', 
                        setpointings=True, 
                        mapsize='2.arcmin', 
                        maptype='hexagonal', 
                        pointingspacing='nyquist', 
                        integration='30s', 
                        totaltime=totaltime,
                        graphics='file',
                        overwrite=True)

            print('-----------------')
            print('First cleaning to obtain the dirty image and estimate the peak')

            dirty_image=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+".dirty.image"
            
            tclean(vis=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+".ms", 
                    imagename=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+".dirty", 
                    imsize=[int(image_size[i])],
                    specmode='mfs',  
                    phasecenter = "J2000 12h00m00.0s -23d00m00.0s",
                    cell=str(pix[i])+'arcsec',
                    gridder='mosaic', 
                    deconvolver='hogbom',
                    weighting='briggs', 
                    robust=0.5, 
                    niter=0,
                    pbcor=False)

            get_convo2target(myskymodel, dirty_image)
            drop_axis(myskymodel+"_conv")
            CASA2fits(myskymodel+"_conv_subimage")
            convert_JypB_JypP(myskymodel+"_conv")
            drop_axis(myskymodel+"_conv.Jyperpix")
            CASA2fits(myskymodel+"_conv.Jyperpix_subimage")

            stats = imstat(imagename=myskymodel+"_conv")
            peak = stats["max"]  
            my_threshold = 0.1*peak[0]

            immath(imagename=[myskymodel+"_conv"], 
                    outfile=myproject+'/mymask', 
                    expr='iif( IM0 >='+str(my_threshold)+', 1., 0.)')
            
            print('-----------------')
            print('Cleaning')

            tclean(vis=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+".ms", 
                    imagename=myproject+'/'+myproject+".alma.cycle5."+str(i+1), 
                    imsize=[int(image_size[i])], 
                    specmode='mfs',  
                    phasecenter = "J2000 12h00m00.0s -23d00m00.0s",
                    cell=str(pix[i])+'arcsec',
                    gridder='mosaic', 
                    deconvolver='hogbom',
                    weighting='briggs', 
                    robust=0.5, 
                    niter=100000000,
                    pbcor=False,
                    threshold=str(my_threshold)+'Jy',
                    mask=myproject+'/mymask')

            print('-----------------')
            print('PB correction')

            impbcor(imagename=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+'.image',
                        pbimage=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+'.pb',
                        outfile=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+'.image.pbcor',
                        cutoff=0.2)

            drop_axis(myproject+'/'+myproject+'.alma.cycle5.'+str(i+1)+'.image.pbcor')
            CASA2fits(myproject+'/'+myproject+'.alma.cycle5.'+str(i+1)+'.image.pbcor_subimage')
            convert_JypB_JypP(myproject+'/'+myproject+'.alma.cycle5.'+str(i+1)+'.image.pbcor')
            drop_axis(myproject+'/'+myproject+'.alma.cycle5.'+str(i+1)+'.image.pbcor.Jyperpix')
            CASA2fits(myproject+'/'+myproject+'.alma.cycle5.'+str(i+1)+'.image.pbcor.Jyperpix_subimage')

            psf=myproject+'/'+myproject+'.alma.cycle5.'+str(i+1)+'.psf'
            pb=myproject+'/'+myproject+'.alma.cycle5.'+str(i+1)+'.pb'
            CASA2fits(psf)
            CASA2fits(pb)

            print('-----------------')
            print('Simulations of '+ myproject+ ' completed')

os.system('mv *.logs ./logs/')
os.system('mv *.last ./logs/')
os.system('mv *.pre ./logs/')