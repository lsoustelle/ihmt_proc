# ///////////////////////////////////////////////////////////////////////////////////////////////
# // L. SOUSTELLE, PhD, Aix Marseille Univ, CNRS, CRMBM, Marseille, France
# // 2020/12/27
# // Contact: lucas.soustelle@univ-amu.fr
# ///////////////////////////////////////////////////////////////////////////////////////////////

import sys
import os
import numpy
import nibabel
import scipy.optimize
import scipy.integrate
import scipy.linalg
import time
import multiprocessing
import argparse; from argparse import RawTextHelpFormatter
import subprocess
import collections

def main():
    ## parse arguments
    text_description = "Compute ihMT and MT saturation maps [1] from an ihMT-prepared experiment [2]. Outputs are in percentage unit.\
                        \nNotes: \
                        \n\t1) The estimation assumes a centric-out readout, in steady-state.\
                        \n\t2) Providing a B1 map is strongly recommended. \
                        \n\t3) The number of FLASH repetition in TFL should encompass acceleration settings (e.g., GRAPPA). \
                        \n\t4) Volume indices (see --R/--S/--D options) start at 1. \
                        \n\t5) To Siemens users using the C2P ihMT-RAGE sequence (Aix Marseille Univ, CNRS, CRMBM): \
                        \n\t - The last BTR duration can be calculated as: \
                        \n\t\t BTR_last = [ihMT Prep. time] - [# of bursts - 1]*[Burst TR]. \
                        \n\t - The TFL \"FLASH TR\" is equal to the [Echo Spacing] variable in Sequence/Part 1. \
                        \n\t - The TFL \"FLASH number of repetition\" is related to the [Turbo Factor] (TF) in Sequence/Part 1: \
                        \n\t\t - Linear encoding with integrated/separated GRAPPA: TF \
                        \n\t\t - Linear Rot. encoding with integrated GRAPPA: floor(TF/2)-1 + #ACS/2 \
                        \n\t\t - Linear Rot. encoding with separated GRAPPA: floor(TF/2)-1 \
                        \nReferences:\
                        \n\t [1] G. Helms et al., High-resolution maps of magnetization transfer with inherent correction for RF inhomogeneity and T1 relaxation obtained from 3D FLASH MRI, MRM 2008;60:1396-1407 \
                        \n\t [2] F. Munsch et al., Characterization of the cortical myeloarchitecture with inhomogeneous magnetization transfer imaging (ihMT), NeuroImage 2021;225:117442 \
                        " 
    parser = argparse.ArgumentParser(description=text_description,formatter_class=RawTextHelpFormatter)
    parser.add_argument('ihMT',         help="Input 4D NIfTI path")
    parser.add_argument('T1',           help="Input T1 (in sec) NIfTI path")
    parser.add_argument('ihMTsat',      help="Output ihMTsat NIfTI path")
    parser.add_argument('ihMTparx',   nargs="?",help="ihMT preparation parameters (comma-separated), in this order:   \n"
                                                                "\t *** ihMT-RAGE: \n"
                                                                "\t     1) Interpulse delay (\"delta_t\"; ms) \n"
                                                                "\t     2) Number of pulses per burst (int) \n"
                                                                "\t     3) Total number of bursts (int) \n"
                                                                "\t     4) Burst TR (BTR; ms) \n"
                                                                "\t     5) Last Burst TR (ms) \n"
                                                                "\t e.g. 0.8,10,5,100.0,10.0 \n"
                                                                "\n"
                                                                "\t *** ihMT-GRE: \n"
                                                                "\t     1) Interpulse delay (\"delta_t\"; ms) \n"
                                                                "\t     2) Number of pulses per preparation (int) \n"
                                                                "\t     3) Post-ihMT preparation delay (ms) \n"
                                                                "\t e.g. 0.8,10,2.0")
    parser.add_argument('TFLparx',      help="TurboFLASH readout and sequence parameters (comma-separated), in this order:   \n"
                                                                "\t     1) FLASH TR (ms) \n"
                                                                "\t     2) Readout flip angle (deg) \n"
                                                                "\t     3) FLASH number of repetition (int) \n"
                                                                "\t     4) Sequence Time to Repetition (TR; ms) \n"                                                               
                                                                "\t e.g. 8.3,7,128,4000.0")
    parser.add_argument('--MTdsat',     nargs="?",help="Output MTsat from dual-offset images NIfTI path")
    parser.add_argument('--MTssat',     nargs="?",help="Output MTsat from single-offset images NIfTI path")
    parser.add_argument('--ihMTsatB1sq',nargs="?",help="Output ihMTsat image normalized by squared B1 NIfTI path")
    parser.add_argument('--MTdsatB1sq', nargs="?",help="Output MTsat from single-offset images normalized by squared B1 NIfTI path")
    parser.add_argument('--MTssatB1sq', nargs="?",help="Output MTsat from dual-offset images normalized by squared B1 NIfTI path")
    parser.add_argument('--B1',         nargs="?",help="Input B1 map (in absolute unit) NIfTI path")
    parser.add_argument('--mask',       nargs="?",help="Input binary mask NIfTI path")
    parser.add_argument('--R',          nargs="?",help="Reference indices in 4D input (comma-separated integers; default: 1)")
    parser.add_argument('--S',          nargs="?",help="Single-offset indices in 4D input (comma-separated integers; default: 2,4,...,N-1)")
    parser.add_argument('--D',          nargs="?",help="Dual-offset indices in 4D input (comma-separated integers; default: 3,5,...,N)")
    parser.add_argument('--xtol',       nargs="?",type=float, default=1e-6, help="x tolerance for root finding (default: 1e-6)")
    parser.add_argument('--nworkers',   nargs="?",type=int, default=1, help="Use this for multi-threading computation (default: 1)")
    
    args                = parser.parse_args()
    ihMT_in_niipath     = args.ihMT
    ihMTsat_out_niipath = args.ihMTsat 
    T1map_in_niipath    = args.T1
    B1_in_niipath       = args.B1
    mask_in_niipath     = args.mask
    xtolVal             = args.xtol
    NWORKERS            = args.nworkers if args.nworkers <= get_physCPU_number() else get_physCPU_number()

    #### Sequence parx 
    print('')
    print('--------------------------------------------------')
    print('----- Checking entries for ihMTsat processing ----')
    print('--------------------------------------------------')
    print('')        
    
    args.ihMTparx = args.ihMTparx.split(',')        
    if len(args.ihMTparx) == 5: 
        EXP_RAGE_GRE = 'RAGE'; ## RAGE experiment
        ihMTparx_NT             = collections.namedtuple('ihMTparx_NT','dt Np NBTR BTR BTRlast')
        ihMTparx = ihMTparx_NT(dt       = float(args.ihMTparx[0])*1e-3, 
                               Np       = int(args.ihMTparx[1]),
                               NBTR     = int(args.ihMTparx[2]),
                               BTR      = float(args.ihMTparx[3])*1e-3,
                               BTRlast  = float(args.ihMTparx[4])*1e-3)
        print('Summary of input ihMT-RAGE preparation parameters:')
        print('\t Delta_t: {:.1f} ms'.format(ihMTparx.dt*1e3))
        print('\t Number of pulses per burst: {}'.format(ihMTparx.Np))
        print('\t Number of bursts: {}'.format(ihMTparx.NBTR))
        print('\t Burst TR: {:.1f} ms'.format(ihMTparx.BTR*1e3))
        print('\t Last Burst TR: {:.1f} ms'.format(ihMTparx.BTRlast*1e3))
        print('')
    elif len(args.ihMTparx) == 3: 
        EXP_RAGE_GRE = 'GRE'
        ihMTparx_NT             = collections.namedtuple('ihMTparx_NT','dt Np pihMTdelay')
        ihMTparx = ihMTparx_NT(dt           = float(args.ihMTparx[0])*1e-3, 
                               Np           = int(args.ihMTparx[1]),
                               pihMTdelay   = float(args.ihMTparx[2])*1e-3)
        print('Summary of input ihMT-GRE preparation parameters:')
        print('\t Delta_t: {:.1f} ms'.format(ihMTparx.dt*1e3))
        print('\t Number of pulses per preparation: {}'.format(ihMTparx.Np))
        print('\t Post-ihMT preparation delay: {}'.format(ihMTparx.pihMTdelay*1e3))
        print('')    
    else: 
        parser.error('Wrong amount of ihMT preparation parameters (ihMTparx \
                             --- expected 3 (ihMT-GRE) or 5 (ihMT-RAGE), found {})' .format(len(args.ihMTparx)))
    
    args.TFLparx = args.TFLparx.split(',')
    if len(args.TFLparx) != 4: 
        parser.error('Wrong amount of sequence/readout parameters (TFLparx \
                         --- expected 4, found {})' .format(len(args.TFLparx)))
    TFLparx_NT              = collections.namedtuple('TFLparx_NT', 'subTR FA Npart TR')
    TFLparx = TFLparx_NT(  subTR    = float(args.TFLparx[0])*1e-3, 
                           FA       = float(args.TFLparx[1]),
                           Npart    = int(args.TFLparx[2]),
                           TR       = float(args.TFLparx[3])*1e-3)    
    print('Summary of TFL parameters:')
    print('\t subTR: {:.1f} ms'.format(TFLparx.subTR*1e3))
    print('\t Readout FA: {:.1f} deg'.format(TFLparx.FA))
    print('\t Number of TFL rep.: {}'.format(TFLparx.Npart))
    print('\t Sequence TR: {:.1f} ms'.format(TFLparx.TR*1e3))
    print('')
    
    # last check 
    for field in ihMTparx._fields:
        if(getattr(ihMTparx, field) < 0):
            parser.error('All ihMTparx values should be positive')
    for field in TFLparx._fields:
        if(getattr(TFLparx, field) < 0):
            parser.error('All TFLparx values should be positive')
    
    
    #### check input data
    # check ihMT data
    if not os.path.isfile(ihMT_in_niipath) or len(nibabel.load(ihMT_in_niipath).shape) != 4:
        parser.error('Volume {} does not exist or is not 4D' .format(ihMT_in_niipath))
    print('ihMT provided volume exist')
        
    # check T1 map
    if not os.path.isfile(T1map_in_niipath):
        parser.error('T1 map volume {} does not exist' .format(T1map_in_niipath))
    print('T1 map provided volume exist')
    
    # check B1 map
    if args.B1 is None:
        print('No B1 map provided (this is highly not recommended)')
    elif args.B1 is not None and not os.path.isfile(B1_in_niipath):
        parser.error('B1 map volume {} does not exist' .format(B1_in_niipath))
    else:
        print('B1 map provided volume exist')
        
        
    # check mask
    if args.mask is None:
        print('No mask provided')
    elif args.mask is not None and not os.path.isfile(mask_in_niipath):
        parser.error('Mask map volume {} does not exist' .format(mask_in_niipath))
    else:
        print('Mask provided volume exist')
        

    #### load data
    # get ihMT data
    ihMT_data = nibabel.load(ihMT_in_niipath).get_fdata()

    # get B1 data
    if args.B1 is not None:
        B1_map = nibabel.load(B1_in_niipath).get_fdata()
    else:
        B1_map = numpy.ones(ihMT_data.shape[0:3])
        
    # get T1 data
    T1_map = nibabel.load(T1map_in_niipath).get_fdata()

    # get indices to process from mask
    if args.mask is not None:
        mask_data   = nibabel.load(mask_in_niipath).get_fdata()
    else:
        mask_data   = numpy.ones(ihMT_data.shape[0:3])
    mask_idx = numpy.asarray(numpy.where(mask_data == 1))


    #### Sorting and/or exclude with indices -R/-S/-D 
    # Default R=1; S=2,4,...N-1; D=3,5,...,N
    R_idx = args.R
    S_idx = args.S
    D_idx = args.D
    if args.R is not None:
        R_idx = numpy.array(args.R.split(','), dtype=numpy.int32)-1
    else:
        R_idx = [0]
    if args.S is not None:
        S_idx = numpy.array(args.S.split(','), dtype=numpy.int32)-1
    else:
        S_idx = numpy.arange(1,ihMT_data.shape[3],2)
    if args.D is not None:
        D_idx = numpy.array(args.D.split(','), dtype=numpy.int32)-1
    else:
        D_idx = numpy.arange(2,ihMT_data.shape[3],2)
    
    
    #### build common xData elements
    T1_data     = T1_map[mask_idx[0],mask_idx[1],mask_idx[2]][numpy.newaxis,:].T
    B1_data     = B1_map[mask_idx[0],mask_idx[1],mask_idx[2]][numpy.newaxis,:].T
    
    cosFA_RO    = numpy.cos(TFLparx.FA * B1_data *numpy.pi/180)
    Np          = numpy.full(T1_data.shape[0],ihMTparx.Np)[numpy.newaxis,:].T
    Npart       = numpy.full(T1_data.shape[0],TFLparx.Npart)[numpy.newaxis,:].T 
    E1_subTR    = numpy.exp(numpy.divide(-TFLparx.subTR,T1_data, \
                                         out=numpy.ones(T1_data.shape, dtype=float), where=T1_data!=0))
        
    ################ MT0 estimation
    #### build xData_MT0
    if EXP_RAGE_GRE == 'RAGE':
        E1_ihMT     = numpy.exp(numpy.divide(-(ihMTparx.BTR*(ihMTparx.NBTR-1)+ihMTparx.BTRlast),T1_data, \
                                             out=numpy.ones(T1_data.shape, dtype=float), where=T1_data!=0)) 
        E1_RD       = numpy.exp(numpy.divide(-(TFLparx.TR - (ihMTparx.NBTR-1)*ihMTparx.BTR - ihMTparx.BTRlast - TFLparx.Npart*TFLparx.subTR),T1_data, \
                                             out=numpy.ones(T1_data.shape, dtype=float), where=T1_data!=0))     
    else:
        E1_ihMT     = numpy.exp(numpy.divide(-(ihMTparx.Np*ihMTparx.dt+ihMTparx.pihMTdelay),T1_data, \
                                             out=numpy.ones(T1_data.shape, dtype=float), where=T1_data!=0))        
        E1_RD       = numpy.exp(numpy.divide(-(TFLparx.TR - ihMTparx.Np*ihMTparx.dt - ihMTparx.pihMTdelay - TFLparx.Npart*TFLparx.subTR),T1_data, \
                                             out=numpy.ones(T1_data.shape, dtype=float), where=T1_data!=0))  
            
    xData_MT0  = [*zip(numpy.hstack((E1_ihMT,E1_subTR,E1_RD,cosFA_RO,Npart)))]
    
    #### run pre-compute MT0 once & for all
    print('')
    print('--------------------------------------------------')
    print('-------- Proceeding to MT0 pre-computation -------')
    print('--------------------------------------------------')
    print('')
    
    start_time = time.time()
    with multiprocessing.Pool(NWORKERS) as pool:
        RES     = pool.starmap(func_MT0,xData_MT0)    
    delay = time.time()
    print("---- Done in {} seconds ----" .format(delay - start_time))
    Mz_MT0 = numpy.array([a_tup[0] for a_tup in RES],dtype=float)[numpy.newaxis,:].T
    A_RAGE = numpy.array([a_tup[1] for a_tup in RES],dtype=float)[numpy.newaxis,:].T
    B_RAGE = numpy.array([a_tup[2] for a_tup in RES],dtype=float)[numpy.newaxis,:].T
    

    ################ MTsat estimation
    #### build xData_MTw
    E1_dt       = numpy.exp(numpy.divide(-ihMTparx.dt,T1_data, \
                                         out=numpy.ones(T1_data.shape, dtype=float), where=T1_data!=0))
    if EXP_RAGE_GRE == 'RAGE':
        NBTR        = numpy.full(T1_data.shape[0],ihMTparx.NBTR)[numpy.newaxis,:].T 
        E1_BTR      = numpy.exp(numpy.divide(-(ihMTparx.BTR - ihMTparx.Np*ihMTparx.dt),T1_data, \
                                             out=numpy.ones(T1_data.shape, dtype=float), where=T1_data!=0))
        E1_BTRlast  = numpy.exp(numpy.divide(-(ihMTparx.BTRlast - ihMTparx.Np*ihMTparx.dt),T1_data, \
                                             out=numpy.ones(T1_data.shape, dtype=float), where=T1_data!=0))
            
        xData       = numpy.hstack((Mz_MT0,A_RAGE,B_RAGE,E1_dt,E1_BTR,E1_BTRlast,Np,NBTR))
    else:
        E1_pihMTdelay = numpy.exp(numpy.divide(-ihMTparx.pihMTdelay,T1_data, \
                                             out=numpy.ones(T1_data.shape, dtype=float), where=T1_data!=0))
            
        xData       = numpy.hstack((Mz_MT0,A_RAGE,B_RAGE,E1_dt,E1_pihMTdelay,Np))
    
    #### build yData
    MT0_data    = numpy.mean(ihMT_data[:,:,:,R_idx],axis=3)
    MTs_data    = numpy.mean(ihMT_data[:,:,:,S_idx],axis=3)
    MTd_data    = numpy.mean(ihMT_data[:,:,:,D_idx],axis=3)
    
    MT0_ydata   = MT0_data[mask_idx[0],mask_idx[1],mask_idx[2]][numpy.newaxis,:].T
    MTs_yData   = MTs_data[mask_idx[0],mask_idx[1],mask_idx[2]][numpy.newaxis,:].T / MT0_ydata
    MTd_yData   = MTd_data[mask_idx[0],mask_idx[1],mask_idx[2]][numpy.newaxis,:].T / MT0_ydata
    
    #### iterable lists
    EXP_RAGE_GRE = [EXP_RAGE_GRE]   * xData.shape[0]
    xtolVal      = [xtolVal]        * xData.shape[0]
    MTs_iterable = [*zip(xData,MTs_yData,EXP_RAGE_GRE,xtolVal)]
    MTd_iterable = [*zip(xData,MTd_yData,EXP_RAGE_GRE,xtolVal)]
        
    #### run
    print('')
    print('--------------------------------------------------')
    print('--------- Proceeding to MTssat estimation --------')
    print('--------------------------------------------------')
    print('')
    
    start_time = time.time()
    with multiprocessing.Pool(NWORKERS) as pool:
        MTssat = pool.starmap(fit_MTsat_brentq,MTs_iterable)
    delay = time.time()
    print("---- Done in {} seconds ----" .format(delay - start_time))
    MTssat = numpy.array(MTssat,dtype=float)
    
    print('')
    print('--------------------------------------------------')
    print('--------- Proceeding to MTdsat estimation --------')
    print('--------------------------------------------------')
    print('')
    
    start_time = time.time()
    with multiprocessing.Pool(NWORKERS) as pool:
        MTdsat = pool.starmap(fit_MTsat_brentq,MTd_iterable)
    delay = time.time()
    print("---- Done in {} seconds ----" .format(delay - start_time))
    MTdsat = numpy.array(MTdsat,dtype=float)
    
    
    ################ store & save NIfTI(s)
    ref_nii = nibabel.load(ihMT_in_niipath)
    ihMTsat_map = numpy.full(ref_nii.shape[0:3],0,dtype=float)
    ihMTsat_map[mask_idx[0],mask_idx[1],mask_idx[2]] = 2*(MTdsat - MTssat)*100
    ihMTsat_map[(ihMTsat_map < 0) | (ihMTsat_map > 100)] = 0
    new_img = nibabel.Nifti1Image(ihMTsat_map, ref_nii.affine, ref_nii.header)
    nibabel.save(new_img, ihMTsat_out_niipath)
    
    if args.MTdsat is not None:
        MTdsat_map = numpy.full(ref_nii.shape[0:3],0,dtype=float)
        MTdsat_map[mask_idx[0],mask_idx[1],mask_idx[2]] = MTdsat*100
        new_img = nibabel.Nifti1Image(MTdsat_map, ref_nii.affine, ref_nii.header)
        nibabel.save(new_img, args.MTdsat)
    if args.MTssat is not None:
        MTssat_map = numpy.full(ref_nii.shape[0:3],0,dtype=float)
        MTssat_map[mask_idx[0],mask_idx[1],mask_idx[2]] = MTdsat*100
        new_img = nibabel.Nifti1Image(MTssat_map, ref_nii.affine, ref_nii.header)
        nibabel.save(new_img, args.MTssat)
    
    if args.ihMTsatB1sq is not None and args.B1 is not None:
        ihMTsatB1sq_map = numpy.full(ref_nii.shape[0:3],0,dtype=float)
        ihMTsatB1sq_map = numpy.divide(ihMTsat_map, B1_map**2, out=numpy.zeros(ihMTsat_map.shape, dtype=float), where=B1_map!=0)
        new_img = nibabel.Nifti1Image(ihMTsatB1sq_map, ref_nii.affine, ref_nii.header)
        nibabel.save(new_img, args.ihMTsatB1sq)
    
    if args.MTdsatB1sq is not None and args.B1 is not None:
        MTdsatB1sq_map = numpy.full(ref_nii.shape[0:3],0,dtype=float)
        MTdsatB1sq_map[mask_idx[0],mask_idx[1],mask_idx[2]] = MTdsat*100
        MTdsatB1sq_map = numpy.divide(MTdsat_map, B1_map**2, out=numpy.zeros(MTdsat_map.shape, dtype=float), where=B1_map!=0)
        new_img = nibabel.Nifti1Image(MTdsatB1sq_map, ref_nii.affine, ref_nii.header)
        nibabel.save(new_img, args.MTdsatB1sq)
        
    if args.MTssatB1sq is not None and args.B1 is not None:
        MTssatB1sq_map = numpy.full(ref_nii.shape[0:3],0,dtype=float)
        MTssatB1sq_map[mask_idx[0],mask_idx[1],mask_idx[2]] = MTssat*100
        MTssatB1sq_map = numpy.divide(MTssat_map, B1_map**2, out=numpy.zeros(MTssat_map.shape, dtype=float), where=B1_map!=0)
        new_img = nibabel.Nifti1Image(MTssatB1sq_map, ref_nii.affine, ref_nii.header)
        nibabel.save(new_img, args.MTssatB1sq)
        
   
###################################################################
############## Fitting-related functions
###################################################################
def func_compute_AB(TILT,E1_SHORT,E1_LONG,N):
    ## compute closed forms of Mz+1=A*Mz+B for a N*(TILT+E1_SHORT)+E1_LONG sequence pattern 
    A = E1_LONG * (E1_SHORT*TILT)**N
    B = 1 + E1_LONG * ( E1_SHORT*TILT - (E1_SHORT*TILT)**N ) / ( 1 - E1_SHORT*TILT ) - \
        E1_LONG/TILT * ( E1_SHORT*TILT * (1 - (E1_SHORT*TILT)**N ) / ( 1 - E1_SHORT*TILT ) )
    return A,B

def func_MT0(xData):
    E1_ihMT,E1_subTR,E1_RD,cosFA_RO,Npart = xData
    A_RAGE,B_RAGE = func_compute_AB(cosFA_RO,E1_subTR,E1_RD,Npart)
    A_ihMT,B_ihMT = E1_ihMT,1-E1_ihMT
    
    Mz_MT0 = (B_ihMT + A_ihMT*B_RAGE) / (1 - A_ihMT*A_RAGE)
    
    return Mz_MT0,A_RAGE,B_RAGE

def func_MTsat_RAGE_root(delta,xData,yData):
    Mz_MT0,A_RAGE,B_RAGE,E1_dt,E1_BTR,E1_BTRlast,Np,NBTR = xData
    A_BTR,B_BTR   = func_compute_AB(1-delta,E1_dt,E1_BTR,Np)
    A_BTRL,B_BTRL = func_compute_AB(1-delta,E1_dt,E1_BTRlast,Np)
    
    Mz_MTw = (B_BTRL + B_BTR*A_BTRL*(A_BTR - A_BTR**NBTR)/(A_BTR - A_BTR**2) + B_RAGE*A_BTRL*A_BTR**(NBTR-1))/ \
             (1 - A_BTRL*A_RAGE*A_BTR**(NBTR-1))
    
    return Mz_MTw/Mz_MT0 - yData
        
def func_MTsat_GRE_root(delta,xData,yData):
    Mz_MT0,A_RAGE,B_RAGE,E1_dt,E1_pihMTdelay,Np = xData
    A_BTR,B_BTR   = func_compute_AB(1-delta,E1_dt,E1_pihMTdelay,Np)
    
    Mz_MTw = (B_BTR + A_BTR*B_RAGE) / (1 - A_BTR*A_RAGE)
    
    return Mz_MTw/Mz_MT0 - yData
    
def fit_MTsat_brentq(xData,yData,EXP_RAGE_GRE,xtolVal):
    # root-finding with Brent's method
    # see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html#scipy.optimize.brentq
    
    try:
        if EXP_RAGE_GRE == 'RAGE':
            x0 = scipy.optimize.brentq(func_MTsat_RAGE_root, 0, 0.3,
                                        (xData,yData),xtol=xtolVal)
        else:
            x0 = scipy.optimize.brentq(func_MTsat_GRE_root, 0, 0.3,
                                        (xData,yData),xtol=xtolVal)
        return x0
    except:
        return 0
    
    
###################################################################
############## Get CPU info
################################################################### 
def get_physCPU_number():
    # from joblib source code (commit d5c8274)
    # https://github.com/joblib/joblib/blob/master/joblib/externals/loky/backend/context.py#L220-L246
    if sys.platform == "linux":
        cpu_info = subprocess.run(
            "lscpu --parse=core".split(" "), capture_output=True)
        cpu_info = cpu_info.stdout.decode("utf-8").splitlines()
        cpu_info = {line for line in cpu_info if not line.startswith("#")}
        cpu_count_physical = len(cpu_info)
    elif sys.platform == "win32":
        cpu_info = subprocess.run(
            "wmic CPU Get NumberOfCores /Format:csv".split(" "),
            capture_output=True)
        cpu_info = cpu_info.stdout.decode('utf-8').splitlines()
        cpu_info = [l.split(",")[1] for l in cpu_info
                    if (l and l != "Node,NumberOfCores")]
        cpu_count_physical = sum(map(int, cpu_info))
    elif sys.platform == "darwin":
        cpu_info = subprocess.run(
            "sysctl -n hw.physicalcpu".split(" "), capture_output=True)
        cpu_info = cpu_info.stdout.decode('utf-8')
        cpu_count_physical = int(cpu_info)
    else:
        raise NotImplementedError(
            "unsupported platform: {}".format(sys.platform))
    if cpu_count_physical < 1:
            raise ValueError(
                "found {} physical cores < 1".format(cpu_count_physical))
    return cpu_count_physical
    
###################################################################
############## main
################################################################### 
if __name__ == "__main__":
    sys.exit(main()) 

   
    
    
    
    
    
    
    
    
    
    
    
    
    