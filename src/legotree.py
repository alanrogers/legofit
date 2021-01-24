#! /usr/bin/env python3
###
#@file legotree.py
#@page legotree
#@brief Draw phylogenetic trees from lgo input files
#
#legotree.py: Draw phylogenetic trees from lgo input files
#==========================================================
#
#Legotree constructs a tree from legofit input files. Legotree is
#implemented in python and requires python packages numpy and matplotlib. 
#   
#Legotree is run with the command   
#   
#   legotree --lgo input.lgo
#
#\image html figure_1.png
#
#The basic command will produce an appropriately scaled tree. Admixture 
#events are shown in colored lines. By default, the tree is scaled to 
#generations. The timescale of the tree can be changed with the command 
#line argument "gentime," converting generations to desired generation
#length in years. 
#
#   legotree input.lgo --gentime 25
#
#\image html figure_2.png
#
#The text on the tree can be adjusted with "textsize"
#
#   legotree input.lgo --gentime 25 --textsize 18
#
##\image html figure_3.png
#
#If modelling populations with deep ancestry, you may want to scale the 
#tree so recent history events are more apparent visually. You can do 
#this with the "shrink" command. When this option is used, the labels 
#are moved from the y-axis onto the tree. This option is off by default.
#
#   legotree input.lgo --gentime 25 --textsize 18 --shrink True
#
#Additionally you can use the --method flag to specify the type of shrink.
#By default this is set to "log" and will take the natural log of the 
#time in generations or years. Setting --method to "long" will shrink long 
#branches near the root of the tree only.
#
#
##\image html figure_4.png
#
#By default legotree will include any admixture events you include in
#the lgo file. Sometimes you may want to include admixture events 
#in your lgo file that have an admixture value of 0. If you would like
#to omit those admixture events from your tree you can use the option
#"allmix," which is on by default. 
#
#   legotree input.lgo --gentime 25 --textsize 18 --shrink True --allmix False
#
##\image html figure_5.png
#
#Finally, legotree displays the tree rather than saving it as a file.
#This allows you to manipulate the size of the figure before saving. 
#If you wish to save directly, you can turn the "show" option to a
#file name. 
#
#   legotree input.lgo --gentime 25 --textsize 18 --shrink True --allmix False
#       --show outputtree.png
#
#Any questions or bug reports can be sent to
#   nathan.harris@anthro.utah.edu

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
from shapely.geometry.polygon import Polygon 
from descartes import PolygonPatch
import re
from random import randint

def get_slope(a,b):
    return (b.y - a.y) / (b.x - a.x)

def sub_input(input,string):
    inds = [i for i,j in enumerate(input[:-1]) if j.split()[0] == string]
    for i in inds:
        counter = 0
        while input[i+counter].split()[-1] == "-" or input[i+counter].split()[-1] == "+":
            if counter != 0:
                input[i] += input[i+counter]
            counter += 1
        if counter != 0:
            input[i] += input[i+counter]

    temp = [i for i in input if len(i.split()) > 0 and i.split()[0] == string]
    temp = [re.sub("="," = ",i) for i in temp]
    temp = [re.sub(" +"," ",i) for i in temp] #replace multiple spaces with one space
    return temp

def node_down(x):
    while len(x.children) == 1:
        x = x.children[0]
    return x

def node_up(x):
    while True:
        x = x.parents[0]
        if len(x.children) != 1:
            break
    return x

class piece:
    def __init__(self,name):
        self.name = name
        self.parents = []
        self.children = []
        self.x = None
        self.y = None
        self.time_label = None

class lego_tree:
    def __init__(self,name,input):
        self.name = name
        self.derive = sub_input(input,"derive")
        self.segs = np.array([i.split()[1] for i in input if len(i.split()) > 0 and i.split()[0] == 'segment'])
        self.mix = sub_input(input,"mix")
        self.pieces = np.array([piece(i) for i in self.segs]) #This should be something iterable containing objects of class piece
        
        self.N = sub_input(input,"twoN")
        self.times = sub_input(input,"time")
        self.full_seg = sub_input(input,"segment")

        parent = {}
        child = {}
        for i in self.derive:
            parent[i.split()[1]] = [i.split()[-1]]
            try:
                kid = child[i.split()[-1]] #i.split()[-1] is the parent
                child[i.split()[-1]].append(i.split()[1]) #check if parent has a child entry, and append child if so
            except:
                child[i.split()[-1]] = [i.split()[1]] #if not create parent and attach child to it
        if len(self.mix) > 0: 
            for i in self.mix:
                parent[i.split()[1] ] = [i.split()[3]]
                try:
                    kid = child[i.split()[3]]
                    child[i.split()[3]].append(i.split()[1])
                except:
                    child[i.split()[3]] = [i.split()[1]]

        for i in self.pieces:
            try:
                i.parents = parent[i.name]
            except:
                pass
            try:
                i.children = child[i.name]
            except:
                pass
            i.parents = [self.pieces[np.where(self.segs == j)[0][0]] for j in i.parents]
            i.children = [self.pieces[np.where(self.segs == j)[0][0]] for j in i.children]

        self.tips = [i for i in self.pieces if len(i.children) == 0]
        self.nodes = [i for i in self.pieces if len(i.children) ==2]
        
        #assign each tip an x value. These are arbitrary,  but affect the tree visually
        for i,j in enumerate(self.tips):
            j.x = i

        PCA = {}
        self.pcs = sub_input(input,"param")
        for i in self.pcs:
            PCA[i.split()[-3]] = float(i.split()[-1])

        time = {}
        for i in self.times:
            if i.split()[1] != 'constrained':
                time[i.split()[2]] = float(i.split()[4])
        for i in self.times:
            if i.split()[1] == 'constrained':
                if "pc" in i:
                    time[i.split()[2]] = int(eval(''.join(i.split('#')[0].split()[4:]),PCA))
                else:
                    time[i.split()[2]] = int(eval(''.join(i.split('#')[0].split()[4:]),time)) #feed it a dictionary for evalulating

        seg_time = {}
        for i in self.full_seg:
            seg_time[i.split()[1]] = i.split()[4]
        for i in self.pieces:
            i.y = time[seg_time[i.name]]
        for node in self.nodes:
            node.x = np.mean([node_down(i).x for i in node.children])
        #A parent is also between its children

        for i in self.pieces:
            i.y *= gentime
            i.time_label = i.y

        if len(self.mix) > 0:
            self.mids = np.array([i for i in self.pieces if i not in self.nodes and i not in self.tips])
            for i in self.mids:
                z = np.polyfit([node_up(i).x,node_down(i).x],[node_up(i).y,node_down(i).y],1)
                i.x = (i.y - z[1]) / z[0]

            self.mixes={}
            for m in self.mix:
                self.mixes[m.split()[1]] = [ m.split()[-1], m.split()[-3] ]

            oldest_mix = np.max([time[seg_time[i]] for i in self.mixes.keys()])
            node_times = np.sort([i.y for i in self.nodes])
            self.adj_time = node_times[np.searchsorted(node_times,oldest_mix,side='right')]

            self.frac = {}
            self.fracs = sub_input(input,"mixFrac")
            for i in self.fracs:
                self.frac[i.split()[2]] = float(i.split()[4])

            self.mixfrac = {}
            for i in self.mix:
                self.mixfrac[i.split()[1]] = i.split()[5]

def tellopt(lbl, msg):
    print("  %-20s %s" % (lbl, msg), file=sys.stderr)

def usage():
    print("Usage: legotree.py [options] <input_file>\n", file=sys.stderr)
    print("where options may include:", file=sys.stderr)
    tellopt("--log", "graph time on log scale")
    tellopt("--nomixtimes", "don't print admixture times")
    tellopt("--allmix", "draw arrows even if admixture is zero")
    tellopt("--textsize <x>",
            "set text size to integer <x> points (def: 12)")
    tellopt("--gentime <x>",
            "print time in years assuming <x> y/generation")
    tellopt("--width <x>", "set width of tree branches (def: 27)")
    tellopt("--nolegend", "suppress legend")
    tellopt("--outfile <xxx.png>",
             "specify name of output file (def: screen)")
    tellopt("-h or --help", "print this message")
    sys.exit(1)
    

if __name__ == "__main__":

    logscale = False
    mixtimes = True
    allmix = False
    legend = True
    arrow = True
    textsize = 12
    gentime = 1
    width = 27
    infile = None
    outfile = None

    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == "--log":
            logscale = True
        elif sys.argv[i] == "--nomixtimes":
            mixtimes = False
        elif sys.argv[i] == "--allmix":
            allmix = True
        elif sys.argv[i] == "--textsize":
            i += 1
            if i == len(sys.argv):
                print("Missing arg to --textsize", file=sys.stderr)
                usage()
            textsize = int(sys.argv[i])
        elif sys.argv[i] == "--gentime":
            i += 1
            if i == len(sys.argv):
                print("Missing arg to --gentime", file=sys.stderr)
                usage()
            gentime = float(sys.argv[i])
        elif sys.argv[i] == "--width":
            i += 1
            if i == len(sys.argv):
                print("Missing arg to --width", file=sys.stderr)
                usage()
            width = round(float(sys.argv[i]))
        elif sys.argv[i] == "--nolegend":
            legend = False
        elif sys.argv[i] == "--outfile":
            i += 1
            if i == len(sys.argv):
                print("Missing arg to --outfile", file=sys.stderr)
                usage()
            outfile = sys.argv[i]
        elif sys.argv[i] == "-h" or sys.argv[i] == "--help":
            usage()
        else:
            if infile != None:
                print("Only 1 input file is allowed", file=sys.stderr)
                usage()
            infile = sys.argv[i]
        i += 1
        
    if infile == None:
        print("Missing input file", file=sys.stderr)
        usage()
        
        
    colors = ['#e6194b', '#3cb44b', '#4363d8', '#f58231', '#911eb4', '#f032e6', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff']
    #badcolors = ['#fabebe', '#bcf60c', '#46f0f0']
    tree = lego_tree("tree",open(infile).readlines())
    get_obj = {}
    for i in tree.pieces:
        get_obj[i.name] = i

    ## Check if parents are older than children ##
    for piece in tree.pieces:
        for parent in piece.parents:
            if parent.y <= piece.y:
                print("Age of parent " + parent.name +
                       " <= that of child " +
                       piece.name + ".", file=sys.stderr )
                sys.exit()
            
    ####### And Graph #########
    fig, ax = plt.subplots(figsize=(16,9))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(axis="y",labelsize=textsize)

    ####shrink#####
    if logscale == True:
        for i in tree.pieces:
            if i.y != 0:
                i.y = np.log10(i.y)

    root = [i for i in tree.pieces if len(i.parents) == 0][0]
    plt.plot([root.x,root.x],[root.y,root.y+.1*root.y],
             linewidth = width,
             color='lightgray',zorder=0)  

    for node in tree.nodes:
        if not logscale:
            plt.annotate(node.name,xy=[node.x-0.3,node.y],
                         zorder = len(tree.segs)+3,
                         fontsize = textsize)
        for dest in node.children:
            dest = node_down(dest)
            if dest in tree.tips:
                m,b = np.polyfit([node.x,dest.x],[node.y,dest.y],1)
                final = -root.y
                plt.plot([node.x,dest.x],[node.y,dest.y],
                         solid_capstyle="round",linewidth = width,
                         color = 'lightgray',zorder = 0)#,label = dest.name)
                #ax.set_ylim([0,root.y])
            else:
                plt.plot([node.x,dest.x],[node.y,dest.y],
                         solid_capstyle="round",linewidth = width,
                         color = 'lightgray',zorder = 0)#,label = dest.name)
            
    if len(tree.mix) > 0:
        for col_i,key in enumerate(tree.mixes.keys()):
            m = tree.mixes[key]
            if m[1].startswith("m"):
                print(key, m, file=sys.stderr)
                sys.exit(0)
                pass
            else:
                print("unable to parse mix statement",
                      file=sys.stderr)
                print(key, m, file=sys.stderr)
                sys.exit(1)

            #admixture events occur mid branch, we use polyfit to get
            #the equation of the line so we can get the point where
            #the admixture lines should be.  
            p = tree.pieces[np.where([i.name == key for i in \
                                      tree.pieces])[0][0]].parents[0]
            q = tree.pieces[np.where([i.name == tree.mixes[key][0] for \
                                      i in tree.pieces])[0][0]]

            if allmix == False and \
            tree.frac[tree.mixfrac[p.children[0].name]] == 0:
                continue
            
            z = np.polyfit([node_up(p).x,node_down(p).x],
                           [node_up(p).y,node_down(p).y],1)
            new_px = (p.y - z[1]) / z[0]
            p.x = new_px

            w = np.polyfit([node_up(q).x,node_down(q).x],
                           [node_up(q).y,node_down(q).y],1)
            new_qx = (q.y - w[1]) / w[0]
            q.x = new_qx
            
            plt.plot([p.x,q.x],[p.y,q.y],color = colors[col_i],alpha =
                     .5, linewidth=3,zorder = len(tree.segs)+1,
                     ls=(0, (6, randint(2,7))), label = str(m[1]))
            
            if legend == True:
                plt.legend(handlelength=1,fontsize=textsize)
            
            ## Arrow ##
            sca = 0.017*np.sum(np.abs(plt.ylim()))

            if arrow == True:
                if q.x > p.x: #Check the direction of admixture
                    adj = 1
                else:
                    adj = -1
                plt.arrow(p.x+0.25*adj, p.y, -0.25*adj, p.y-q.y,
                          length_includes_head=True, head_width=sca,
                          head_length=.1, lw=0, color = colors[col_i])

            if mixtimes == True:
                plt.annotate(p.time_label,xy=[np.mean([p.x,q.x]),p.y],
                             zorder = len(tree.segs)+3,
                             fontsize = textsize,color = colors[col_i])


        for i in tree.mids:
            m,b = np.polyfit([node_up(i).x,node_down(i).x],
                             [node_up(i).y,node_down(i).y],1)
            i.x = (i.y-b) / m
            perp = -1/m
            new_b = i.y - perp * (i.x)

            #and add the segment label
            plt.annotate(i.name,xy=[i.x-0.1,i.y],
                         zorder = len(tree.segs)+3,
                         fontsize = textsize)
   
    plt.ylim([-0.1*root.y,root.y + 0.1 * root.y])
    plt.xlim([-0.5,len(tree.tips)-0.5])
    
    plt.xticks([])

    #Just an adjuster
    h = root.y*0.2

    limitlen = plt.ylim()[1] - plt.ylim()[0]
    for i in tree.tips:
        plt.annotate(i.name,xy=[i.x-.01,i.y-limitlen*.045],
                     fontsize = textsize)

    if gentime == 1:
        if logscale:
            plt.ylabel("$log_{10}$ generations ago",fontsize = textsize)
        else:
            plt.ylabel('Generations ago',fontsize = textsize)
    else:
        if logscale:
            plt.ylabel("$log_{10}$ years ago", fontsize = textsize)
        else:
            plt.ylabel('Thousand years ago',fontsize = textsize)

    plt.tight_layout()
    if outfile==None:
        plt.show()
    else:
        plt.savefig(outfile, dpi=200)

    #This is the white blocks near tips
'''for i,tip in enumerate(tree.tips):
    if i == len(tree.tips) -1:
        plt.fill_between([tip.x -0.33,tip.x+1],[tip.y,tip.y],[tip.y-h,tip.y-h],color = 'white',zorder=2)
    else:
        plt.fill_between([tip.x -0.33,tip.x+0.33],[tip.y,tip.y],[tip.y-h,tip.y-h],color = 'white',zorder=2)
    
    m,b = np.polyfit([tip.x,node_up(tip).x],[tip.y,node_up(tip).y],1)
    final = -root.y

    ring_mixed = Polygon([(tip.x-.3,tip.y), (tip.x+.3,tip.y), ((final-b)/m +.3, final), ((final-b)/m -.3, final)])
    ax = fig.add_subplot(111)
    ring_patch = PolygonPatch(ring_mixed,color='white',zorder=2)
    ax.add_patch(ring_patch)

    line_x = [(tip.y-h-b)/m,(final-b)/m]
    line_y = [tip.y-h,final]

    plt.plot(line_x,line_y,solid_capstyle="round",linewidth = 30,color = 'blue',label = dest.name,zorder =2)
    #plt.scatter(tip.x,tip.y,s=900,marker='v',color='red',zorder=2)'''

'''class data_linewidth_plot():
    def __init__(self, x, y, **kwargs):
        self.ax = kwargs.pop("ax", plt.gca())
        self.fig = self.ax.get_figure()
        self.lw_data = kwargs.pop("linewidth", 1)
        self.lw = 1
        self.fig.canvas.draw()

        self.ppd = 72./self.fig.dpi
        self.trans = self.ax.transData.transform
        self.linehandle, = self.ax.plot([],[],**kwargs)
        if "label" in kwargs: kwargs.pop("label")
        self.line, = self.ax.plot(x, y, **kwargs)
        self.line.set_color(self.linehandle.get_color())
        self._resize()
        self.cid = self.fig.canvas.mpl_connect('draw_event', self._resize)

    def _resize(self, event=None):
        lw =  ((self.trans((1, self.lw_data))-self.trans((0, 0)))*self.ppd)[1] #changed this last index
        #If it is a 1 lw will be in terms of y axis size, 0 will be in terms of x axis size
        if lw != self.lw:
            self.line.set_linewidth(lw)
            self.lw = lw
            self._redraw_later()

    def _redraw_later(self):
        self.timer = self.fig.canvas.new_timer(interval=10)
        self.timer.single_shot = True
        self.timer.add_callback(lambda : self.fig.canvas.draw_idle())
        self.timer.start()'''
