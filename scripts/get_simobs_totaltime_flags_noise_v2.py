# modified DP 25 June 2023
# Modified Ash - 17th Nov 2023 


# Imports
execfile("get_simobs_mods.py")
exec(open("./get_simobs_mods.py").read())

import os 
import numpy as np
import time
import warnings
warnings.filterwarnings('ignore')
# 

# Define
min_baseline =  [14.6, 14.6, 14.6, 14.6, 14.6, 14.6, 64.0, 110.4, 367.6, 244.0]
mrs =           [28.5, 22.6, 16.2, 11.2, 6.7, 4.11, 2.58, 1.42, 0.814, 0.496]
beam_size =     [3.38, 2.3, 1.42, 0.918, 0.545, 0.306, 0.211, 0.096, 0.057, 0.042]
configuration = ['conf1', 'conf2', 'conf3', 'conf4', 'conf5', 'conf6', 'conf7', 'conf8', 'conf9', 'conf10']
uvrange =       ['20m', '20m', '20m', '20m', '20m', '20m', '70m', '115m', '390m', '260m']
rerun = False

# Run simulate and clean
min_baseline = np.array(min_baseline)
mrs = np.array(mrs)
beam_size = np.array(beam_size)

pix = beam_size/5.
imagesize=200 
image_size=np.zeros(10)
base_source_size = []

####
# run_only_conf = ['conf1', 'conf2', 'conf3', 'conf4', 'conf5', 'conf6', 'conf7']
# run_only_conf = ['conf1', 'conf2', 'conf3', 'conf4', 'conf5']
run_only_conf = ['conf8', 'conf9', 'conf10']
run_only_mrs = ['']
run_only_totaltime = ['6']
####


for i in range (0,10):

    base_source_size.append(1.5*mrs[0]*min_baseline[0]/min_baseline[i]/15)
    image_size[i] = int(imagesize*pix[0]/pix[i])

    if image_size[i] % 2 == 0:
        image_size[i] = image_size[i]+0 
    else:
        image_size[i] = image_size[i]-1

for i in range (0,10):
    for j in range (0,15):
        start_time = time.time()

        if configuration[i] not in run_only_conf: 
            print('[INFO] Skipping... %s' %(configuration[i]))
            continue 

        source_size = base_source_size[i]*(j+1)
        
        myskymodel_dir = './skymodels/'
        myproject_dir = './'

        for totaltime in run_only_totaltime:

            myproject = myproject_dir+configuration[i]+'_'+str('%.1f'%(source_size))+'mrs0_project_%stotaltime_flagged_noise' %totaltime
            myskymodel = myskymodel_dir+configuration[i]+'_'+'ref_image'+str('%.1f'%(source_size))+'mrs0'

            if run_only_mrs != ['']: 
                if str('%.1f'%(source_size)) not in run_only_mrs: 
                    print('[INFO] Skipping... %s' %(str('%.1f'%(source_size))+'mrs0'))
                    continue        
            
            if rerun: 
                os.system('rm -rf %s' %myproject)
                os.system('mkdir %s' %myproject)
            elif os.path.exists(myproject): 
                print('[INFO] rerun = False, skipping myproject %s' %myproject)
                continue
            else:
                os.system('mkdir %s' %myproject)
                
            print('-----------------')
            myskymodel_new = '%s.myskymodel' %(myproject+'/'+myproject+".alma.cycle5."+str(i+1))
            os.system('rm -rf %s' %myskymodel_new)
            if totaltime == '6':
                model_f = 15000
            elif totaltime == '60':
                model_f = 15000*np.sqrt(10) 
            print('[INFO] Model scaling factor: 1/%i' %model_f)
            immath(imagename=myskymodel, expr='IM0/%i' %model_f, outfile=myskymodel_new)
            myskymodel = myskymodel_new

            print('-----------------')
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
                        thermalnoise='tsys-atm', 
                        user_pwv=5.0,
                        seed=11111,
                        setpointings=True, 
                        mapsize='2.arcmin', 
                        maptype='hexagonal', 
                        pointingspacing='nyquist', 
                        integration='30s', 
                        totaltime=totaltime,
                        graphics='file',
                        overwrite=True)

            plotms(vis=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+".noisy.ms",
                    xaxis='uvdist', 
                    yaxis='amp',
                    plotfile=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+".plotms.jpg",
                    showgui=False,
                    overwrite=True,
                    highres=True)

            print('-----------------')
            print('Flagging data')

            if uvrange[i] != '': 
                flagdata(vis=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+".noisy.ms",
                        uvrange='<%s' %uvrange[i])

                plotms(vis=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+".noisy.ms",
                        xaxis='uvdist', 
                        yaxis='amp',
                        plotfile=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+".plotms.flagged.jpg",
                        showgui=False,
                        overwrite=True,
                        highres=True)

            print('-----------------')
            print('First cleaning to obtain the dirty image and estimate the peak')

            dirty_image=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+".dirty.image"
            
            tclean(vis=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+".noisy.ms", 
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
            convert_JypB_JypP(myskymodel+"_conv")

            stats = imstat(imagename=myskymodel+"_conv")
            peak = stats["max"]  
            my_threshold = 0.1*peak[0]
            my_threshold_clean = 0.1*peak[0]

            immath(imagename=[myskymodel+"_conv"], 
                    outfile=myproject+'/mymask', 
                    expr='iif( IM0 >='+str(my_threshold)+', 1., 0.)')
            
            print('-----------------')
            print('Cleaning')

            tclean(vis=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+".noisy.ms", 
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
                    threshold=str(my_threshold_clean)+'Jy',
                    mask=myproject+'/mymask')

            print('-----------------')
            print('PB correction')

            impbcor(imagename=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+'.image',
                        pbimage=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+'.pb',
                        outfile=myproject+'/'+myproject+".alma.cycle5."+str(i+1)+'.image.pbcor',
                        cutoff=0.2)

            convert_JypB_JypP(myproject+'/'+myproject+'.alma.cycle5.'+str(i+1)+'.image.pbcor') 

            exportfits(myproject+'/'+myproject+'.alma.cycle5.'+str(i+1)+'.image.pbcor.Jyperpix', 
                        myproject+'/'+myproject+'.alma.cycle5.'+str(i+1)+'.image.pbcor.Jyperpix.fits',
                        overwrite=True)

            print('-----------------')
            print('Simulations of '+ myproject+ ' completed')
            print("Time taken %0.1f seconds" % (time.time() - start_time))
            print('-----------------')

os.system('mv *.log ./logs/')
os.system('mv *.last ./logs/')
os.system('mv *.pre ./logs/')