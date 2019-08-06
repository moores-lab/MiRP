#!/usr/bin/env python

''' 
For use with the microtubule relion based pipeline (MiRP)

Protocols for smoothening the Rot angle and X/Y shifts
from widely distributed values.

Other handy functions such as per microtubule Euler
and X/Y shift plots, removing microtubules below a certain
particle number, and plotting the confidence of protofilament
and seam position classification, with the option to remove
microtubules below a certain confidence.
'''

__author__ = "Alex Cook" 
__email__ = "alexcook2@outlook.com"
__version__ = "1.0"


from collections import OrderedDict,Counter
from itertools import groupby,combinations
from operator import itemgetter
import ast
import re
import numpy as np
from scipy import stats
import argparse
from matplotlib import pyplot as plt


parser=argparse.ArgumentParser()

parser.add_argument('-ang', 
help = '_data.star on which to smoothen Euler angles')

parser.add_argument('-id', 
help = 'label of angle to smoothen (rlnAngleRot or rlnAnglePsi)')

parser.add_argument('-xy',
help = '_data.star file on which to smoothen X/Y shifts')

parser.add_argument('-p',
help = '_data.star to plot Euler angles and X/Y shifts for' 
       'each microtubule')

parser.add_argument('-c',
help = '_data.star from supervised classification for confidence plot')

parser.add_argument('-s',
help = '_data.star from which to remove microtubules with a certain' 
       'number of particles')

parser.add_argument('-lw', type=float,
help = 'remove microtubules with values below this cutoff (for confidence'
       'and size)')

args=parser.parse_args()


########################################################################
### General Functions ###

###
def parse_star(starfile):

    lines = [i for i in open(starfile,'r')
             if len(i.strip()) != 0 and i[0] != '#']

    starPar=OrderedDict()

    for line in lines:
        
        if line[:5] == 'data_': 
            pass        

        elif line[:5] == 'loop_':    
            pass                 

        elif line[0] == '_': 
            starPar[line.split()[0][1:]] = []

        elif len( line.split() ) == len(starPar): 
            data = line.split()

            for idx,field in enumerate(starPar):
                
                try:
                    val = ast.literal_eval(data[idx])
                except:
                    val = str(data[idx])

                starPar[field].append(val)
                
        else: print('This line could not be processed:\n%s' % line)
        
    return starPar
    
###
def write_star(name,lbls,array):

    outf = open('%s.star' % name,'w')
    
    outf.write('\ndata_\n\nloop_\n')
    
    for i,l in enumerate(lbls):
        outf.write('_%s #%i\n' % (l,i+1) )
        
    for vals in array:
 
        for val in vals:
            outf.write('%s\t' % str(val) )
        outf.write('\n')
    
    outf.close()

###    
def group(array,idx):
    
    grp = groupby( array, key = itemgetter(idx) )
    return [list(i[1]) for i in grp]        

###
def transpose(array):
    return [ list(i) for i in list( zip(*array) ) ]

###
def get_particles(starfile):

    starPar = parse_star(starfile) 
    MetaDataLabels = list( starPar.keys() )
    ptcls = transpose( starPar.values() )

    mgphIDX = MetaDataLabels.index('rlnMicrographName')
    tubeIDX = MetaDataLabels.index('rlnHelicalTubeID')
    trakIDX = MetaDataLabels.index('rlnHelicalTrackLength')

    ptcls.sort( key = itemgetter(mgphIDX, tubeIDX, trakIDX) )

    return ptcls,MetaDataLabels

###
def linear_fit(xax, vals, MTlen):
    
    slp,incpt=stats.linregress(xax, vals)[0:2]

    return [ incpt + slp * x for x in range(1, MTlen + 1) ]

###
def low_cut(MTs, vals, cutoff):

    rm = []

    for idx, p in enumerate(vals):

        if p < cutoff:
            rm.append(idx)
        else:
            pass

    for r in reversed(rm):
        del MTs[r]

    return MTs


########################################################################
''' fuctions to cluster lines with different y-axis intercepts with 
the same - but steep - slopes '''

###    
def cluster_numpy_bins(data,binned): 

    binIDs = np.digitize(data,binned[1])

    clst = {}

    for idx,a in enumerate(binIDs):
        if a in clst.keys():
            for b in clst.keys():
                if a == b:
                    clst[b].append(idx)
                else:
                    pass
        else:
            clst[a] = [idx]

    return clst

###
def flatness(array):

    diffs = []

    for i in range(len(array)-1):
    
        df = array[i]-array[i+1]

        if df < 0: 
            df =- df

        diffs.append(df)

    return sum(diffs)
    
###    
def flatten_and_cluster_shifts(shifts):

    flattening_factors = np.arange(-8, 8, 0.25)
    flatness_score = []
    b = []

    for factor in flattening_factors:
    
        flattened = [ ptcl - idx * factor 
                      for idx,ptcl in enumerate(shifts, start = 1) ] 
        
        try:
            binned = np.histogram(flattened, bins = 'auto')
        except MemoryError:
            binned = np.histogram(flattened, bins = 5)

        bins = cluster_numpy_bins(flattened, binned)

        flatness_score.append( sum([ flatness([flattened[i] for i in bin]) 
                                     for bin in bins.values()
                                   ]
                                  ))
        b.append(bins)

    topbin = max( b[np.argmin(flatness_score)].values(), key = len )

    fit = linear_fit(topbin, [shifts[i] for i in topbin], len(shifts) )

    return fit
        
        
#######################################################################
''' function to cluster lines with shallow slopes, but with many y-axis
intercepts ''' 
        
###
def cluster_shallow_slopes(angles, cutoff):

    linkMtrx = [ (i, i2) for ( (i, j), (i2, j2) ) 
                 in combinations( enumerate(angles), 2 ) 
                 if -cutoff <= float(j) - float(j2) <= cutoff
               ]
    
    cluster = []
        
    while linkMtrx:

        node = linkMtrx[-1]
        
        for idx,pair in reversed( list( enumerate(linkMtrx[:-1]) ) ):
                
            if any( i == j for i, j in combinations(node+pair, 2) ):

                node += pair
                node = tuple( set(node) )

                del linkMtrx[idx]

        del linkMtrx[-1]

        cluster.append( sorted(node) )
    
    if not cluster: return None,None
    
    topclst = cluster.pop( max( enumerate( [len(i) for i in cluster] ), 
                                key=itemgetter(1) )[0] 
                          )

    lwclst = [j for i in cluster for j in i]
        
    return topclst, lwclst
    

########################################################################

if args.ang:
    
    ptcls, MetaDataLabels = get_particles(args.ang)

    mgphIDX = MetaDataLabels.index('rlnMicrographName')
    tubeIDX = MetaDataLabels.index('rlnHelicalTubeID')
    angIDX = MetaDataLabels.index(args.id)

    MTs = [ MT for mgph in group(ptcls, mgphIDX) 
            for MT in group(mgph, tubeIDX) ]

    cutoff = 8

    bad_mts = []

    for mtIDX, MT in enumerate(MTs):
    
        angles = [i[angIDX] for i in MT]

        MTlen = len(MT)
  
        top_clstr, outliers = cluster_shallow_slopes(angles, cutoff)
    
        if not top_clstr:        

            print('MT %i %s cannot be fit, and is discarded' % 
                  ( mtIDX, args.id ) )

            bad_mts.append(mtIDX)

        else:
            top_clstr_vals = [ MT[i][angIDX] for i in top_clstr ]
            fitted = linear_fit( top_clstr, top_clstr_vals, MTlen ) 

            for ptclIDX in range(MTlen):
                MTs[mtIDX][ptclIDX][angIDX] = fitted[ptclIDX]

    for mtIDX in reversed(bad_mts):
        del MTs[mtIDX]

    ptcls = [ptcl for MT in MTs for ptcl in MT]                

    write_star( 'smoothened_%s_data' % args.id, MetaDataLabels, ptcls )


elif args.xy:

    ptcls, MetaDataLabels = get_particles(args.xy)

    mgphIDX = MetaDataLabels.index('rlnMicrographName')
    tubeIDX = MetaDataLabels.index('rlnHelicalTubeID')
    xshIDX = MetaDataLabels.index('rlnOriginX') 
    yshIDX = MetaDataLabels.index('rlnOriginY') 

    MTs = [ MT for mgph in group(ptcls, mgphIDX) 
            for MT in group(mgph, tubeIDX) ]

    for mtID,MT in enumerate(MTs):

        Xsh = [ i[xshIDX] for i in MT]
        Ysh = [ i[yshIDX] for i in MT]

        xax=[ i for i in range(1, len(Xsh) + 1 ) ]

        uniX = flatten_and_cluster_shifts(Xsh)
        uniY = flatten_and_cluster_shifts(Ysh)
                
        for i in range(len(MT)):
            MTs[mtID][i][xshIDX] = uniX[i]
            MTs[mtID][i][yshIDX] = uniY[i]    
        
    ptcls=[ptcl for MT in MTs for ptcl in MT]    
        
    write_star('smoothenedXY_data', MetaDataLabels, ptcls)


elif args.p:

    ptcls, MetaDataLabels = get_particles(args.p)

    mgphIDX = MetaDataLabels.index('rlnMicrographName')
    tubeIDX = MetaDataLabels.index('rlnHelicalTubeID')
    psiIDX=MetaDataLabels.index('rlnAnglePsi')
    theIDX=MetaDataLabels.index('rlnAngleTilt') 
    phiIDX=MetaDataLabels.index('rlnAngleRot') 
    xshIDX=MetaDataLabels.index('rlnOriginX') 
    yshIDX=MetaDataLabels.index('rlnOriginY') 

    MTs=[ MT for mgph in group(ptcls, mgphIDX) 
          for MT in group(mgph, tubeIDX) ]

    for MT in MTs:

        data = tuple( zip( *[ (i[phiIDX],i[theIDX],i[psiIDX],i[xshIDX],i[yshIDX]) 
                              for i in MT
                            ] ) )

        f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize = (12, 6) )    

        plt.setp( [ax1, ax2, ax3, ax4], 
                  xticks = [i for i in range(1, len(data[0]) + 1, 2) ],
                  xlabel = 'Particle Number'
                )

        ax1.plot(data[0], 'o')
        ax1.set_ylim([-181, 181])
        ax1.set_yticks( [i for i in range(-180, 180 + 1, 40) ] )
        ax1.set_ylabel('Phi Angle')
        ax1.set_title('Phi')

        ax2.plot(data[1], 'o')
        ax2.set_yticks( [i for i in range(-180, 180 + 1, 40) ] )
        ax2.set_ylim([-181, 181])
        ax2.set_ylabel('Theta Angle')
        ax2.set_title('Theta')

        ax3.plot(data[2], 'o')
        ax3.set_yticks( [i for i in range(-180, 180 + 1, 40) ] )
        ax3.set_ylim([-181, 181])
        ax3.set_ylabel('Psi Angle')
        ax3.set_title('Psi')

        ax4.plot(data[3], 'o', label = 'Xshift' )
        ax4.plot(data[4], 'o', label = 'Yshift' )
        ax4.set_title('X/Y Shifts')
        ax4.set_ylabel('Shift (pixels)')
        ax4.legend(loc='upper right')

        plt.show()


if args.c:

    ptcls, MetaDataLabels = get_particles(args.c)

    mgphIDX = MetaDataLabels.index('rlnMicrographName')
    tubeIDX = MetaDataLabels.index('rlnHelicalTubeID')
    clsIDX = MetaDataLabels.index('rlnClassNumber')

    MTs = [ MT for mgph in group(ptcls, mgphIDX)
            for MT in group(mgph, tubeIDX) ]

    cer = []

    for MT in MTs:
        CLSdist = Counter( [ptcl[clsIDX] for ptcl in MT] )
        CLStop = CLSdist.most_common(1)[0][1]
        total = len(MT)
        cer.append( ( float(CLStop) / float(total) ) * 100 )

    plt.hist(cer, bins = 10 )
    plt.xlabel('Percent Confidence')
    plt.ylabel('Frequency')
    plt.show()

    if args.lw:
        goodTubes = low_cut(MTs, cer, args.lw)
        ptcls = [ ptcl for tube in goodTubes for ptcl in tube ]
        write_star('microtubule_confidence_data', MetaDataLabels, ptcls)


elif args.s:

    ptcls, MetaDataLabels = get_particles(args.s)

    mgphIDX = MetaDataLabels.index('rlnMicrographName')
    tubeIDX = MetaDataLabels.index('rlnHelicalTubeID')

    MTs = [ MT for mgph in group(ptcls, mgphIDX)
            for MT in group(mgph, tubeIDX) ]

    MT_lens = [len(MT) for MT in MTs]

    ptcls=[ ptcl for tube in low_cut(MTs, MT_lens, args.lw) 
            for ptcl in tube ]

    write_star('microtubule_length', MetaDataLabels, ptcls)
