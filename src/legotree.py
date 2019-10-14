#!/usr/bin/python

###
#@file lego_tree.py
#@page diverg
#@brief Draw phylogenetic trees from lgo input files
#
#diverg.py: Draw phylogenetic trees from lgo input files
#=================================================================
#
#Legotree constructs a tree from legofit input files. Legotree is
#implemented in python and requires python packages numpy and matplotlib. 
#   
#Legotree is run with the command   
#   
#   legotree input.lgo
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
#The text on the tree can be adjusted with "text_size"
#
#   legotree input.lgo --gentime 25 --text_size 18
#
##\image html figure_3.png
#
#If modelling populations with deep ancestry, you may want to scale the 
#tree so recent history events are more apparent visually. You can do 
#this with the "shrink" command. When this option is used, the labels 
#are moved from the y-axis onto the tree. This option is off by default.
#
#   legotree input.lgo --gentime 25 --text_size 18 --shrink True
#
##\image html figure_4.png
#
#By default legotree will include any admixture events you include in
#the lgo file. Sometimes you may want to include admixture events 
#in your lgo file that have an admixture value of 0. If you would like
#to omit those admixture events from your tree you can use the option
#"allmix," which is on by default. 
#
#   legotree input.lgo --gentime 25 --text_size 18 --shrink True --allmix False
#
##\image html figure_5.png
#
#Finally, legotree displays the tree rather than saving it as a file.
#This allows you to manipulate the size of the figure before saving. 
#If you wish to save directly, you can turn the "show" option to a
#file name. 
#
#   legotree input.lgo --gentime 25 --text_size 18 --shrink True --allmix False
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

class data_linewidth_plot():
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
        self.timer.start()

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
            parent[i.split(' ')[1]] = [i.split()[-1]]
            try:
                kid = child[i.split()[-1]]
                child[i.split()[-1]].append(i.split()[1])
            except:
                child[i.split()[-1]] = [i.split()[1]]
        if len(self.mix) > 0: 
            for i in self.mix:
                parent[i.split()[1] ] = [i.split()[3]]
                try:
                    kid = child[i.split()[-1]]
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
                self.mixes[m.split()[1]] = m.split() [-1]

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


if __name__ == "__main__":
    arg_vals = {}
    target =  sys.argv[1]
    print(target)
    try:
        args = sys.argv[2:]
    except:
        pass
    
    inds = [i for i,j in enumerate(args) if "--" in j]
    for i in inds:
        arg_vals[args[i].split("--")[-1]] = args[i+1]
        #arg_vals[i.split('=')[0]] = i.split("=")[1]


    #set defaults
    try:
        shrink = eval(arg_vals["shrink"])
    except:
        shrink=False
    try:
        method = arg_vals["method"]
    except:
        method = "log"
    try:
        allmix=eval(arg_vals['allmix'])
    except:
        allmix = True
    try:
        text_size = float(arg_vals["textsize"])
    except:
        text_size = 12
    try:
        gentime=float(arg_vals["gentime"])
    except:
        gentime = 1.0
    try:
        tlabels = arg_vals["tlabels"].split(',')
    except:
        tlabels = []
    try:
        show = eval(arg_vals["show"])
    except:
        try:
            show = arg_vals["show"]
        except:
            show = True    

    for i in arg_vals.keys():
        if i not in ['shrink','allmix','textsize','gentime','tlabels','show','method']:
            print("Argument " + str(i) + " not understood")

    colors = ['#e6194b', '#3cb44b', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff']
    tree = lego_tree("tree",open(target).readlines())
    get_obj = {}
    for i in tree.pieces:
        get_obj[i.name] = i


    ####### And Graph #########

    #
    fig, ax = plt.subplots(figsize=(16,9))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)


    ####shrink#####
    if shrink == True:
        for i in tree.pieces:
            if method == "log":
                if i.y != 0:
                    i.y = np.log(i.y)
            if i.y > tree.adj_time:
                i.time_label = i.y
                i.y = tree.adj_time + 0.1*i.y
            if i not in tree.mids and i not in tree.tips:
                plt.annotate(i.name +'  ' + str(int(i.time_label)),xy=[i.x+0.1,i.y],fontsize = text_size)
            elif i in tree.tips and i.y !=0:
                plt.annotate(str(int(i.time_label)),xy=[i.x+0.1,i.y],fontsize = text_size)
            #elif i in tree.mids and i.name in tree.mixes.keys() and i.name not in tree.mixfrac.values():
            #    plt.annotate(str(int(i.time_label)),xy=[i.x+0.1,i.y],fontsize = text_size)


        ax.spines['left'].set_visible(False)
        plt.yticks([])
        #plt.yticks([i.time_label for i in pieces])

    root = [i for i in tree.pieces if len(i.parents) == 0][0]
    plt.plot([root.x,root.x],[root.y,root.y+.1*root.y],linewidth = 27,color='lightgray',zorder=0)  
    #plt.plot([root.x,root.x],[root.y,root.y+.1*root.y], label='some 1 data unit wide line', linewidth=7, alpha=1,color='lightgray',zorder=0)
    for node in tree.nodes:
        if shrink != True:
            plt.annotate(node.name,xy=[node.x-0.3,node.y],zorder = len(tree.segs)+3,fontsize = text_size)
        for dest in node.children:
            dest = node_down(dest)
            if dest in tree.tips:
                m,b = np.polyfit([node.x,dest.x],[node.y,dest.y],1)
                final = -root.y
                #plt.plot([node.x,(final-b)/m],[node.y,final],solid_capstyle="round",linewidth = 27,color = 'lightgray',label = dest.name,zorder = 0)
                plt.plot([node.x,dest.x],[node.y,dest.y],solid_capstyle="round",linewidth = 27,color = 'lightgray',label = dest.name,zorder = 0)
                #ax.set_ylim([0,root.y])
            else:
                plt.plot([node.x,dest.x],[node.y,dest.y],solid_capstyle="round",linewidth = 27,color = 'lightgray',label = dest.name,zorder = 0)
            
    if len(tree.mix) > 0:
        for col_i,m in enumerate(tree.mixes.keys()):
            p = tree.pieces[np.where([i.name == m for i in tree.pieces])[0][0]].parents[0]
            q = tree.pieces[np.where([i.name == tree.mixes[m] for i in tree.pieces])[0][0]]

            if allmix == False and tree.frac[tree.mixfrac[p.children[0].name]] == 0:
                continue

            z = np.polyfit([node_up(p).x,node_down(p).x],[node_up(p).y,node_down(p).y],1)
            new_px = (p.y - z[1]) / z[0]
            p.x = new_px

            w = np.polyfit([node_up(q).x,node_down(q).x],[node_up(q).y,node_down(q).y],1)
            new_qx = (q.y - w[1]) / w[0]
            q.x = new_qx

            #if time[seg_time[m]] == 0:
            #    continue
            plt.plot([p.x,q.x],[p.y,q.y],color = colors[col_i],alpha = .5,linewidth=2,zorder = len(tree.segs)+1,ls='--')
            plt.annotate(p.time_label,xy=[np.mean([p.x,q.x]),p.y],zorder = len(tree.segs)+3,fontsize = text_size,color = colors[col_i])

            ### make lines going up and down from the admixture event
            plt.plot([p.x,node_down(p).x],[p.y,node_down(p).y],color = colors[col_i],alpha = .5,linewidth=2,solid_capstyle='round',zorder = len(tree.segs)+1,ls='--')
            
            #plt.plot([q.x,(q.y + .5*(node_up(q).y -q.y)-w[1])/w[0]],[q.y,q.y + .5*(node_up(q).y -q.y)],color = colors[col_i],alpha = .5,linewidth=2,solid_capstyle='round',zorder = len(tree.segs)+1,ls='--')
            #if q.name =='xy2':
            #    break

        for i in tree.mids:
            m,b = np.polyfit([node_up(i).x,node_down(i).x], [node_up(i).y,node_down(i).y],1)
            i.x = (i.y-b) / m
            #if i.x is None:

            #plt.plot([i.x-.2,i.x+.2],[i.y,i.y],linewidth = 3,color = 'black',alpha = 0.4)
            perp = -1/m
        
            new_b = i.y - perp * (i.x)

            ax.plot([i.x-.1,i.x+.1],[perp*(i.x-.1)+new_b,perp*(i.x+.1)+new_b],linewidth = 2,color = 'black',alpha = 0.5)

            #and add the segment label
            plt.annotate(i.name,xy=[i.x-0.1,i.y],zorder = len(tree.segs)+3,fontsize = text_size)
   
    plt.ylim([-0.1*root.y,root.y + 0.1 * root.y])
    plt.xlim([-0.5,len(tree.tips)-0.5])
    
    #moved xlabels onto the tree
    #plt.xticks([i.x for i in tree.tips],[i.name for i in tree.tips],fontsize = text_size)
    plt.xticks([])

    #Just an adjuster
    h = root.y*0.2

    for i in tree.tips:
        plt.annotate(i.name,xy=[i.x-.01,i.y-plt.ylim()[1]*.06],fontsize = text_size)

    if gentime == 1:
        plt.ylabel('Generations ago',fontsize = text_size)
    else:
        plt.ylabel('Years ago',fontsize = text_size)
    #ax.plot(ax.get_xlim(),ax.get_ylim(),ls='--')

    for lab in tlabels:
        i = get_obj[lab]
        plt.annotate(str(int(i.time_label)),xy=[i.x+0.1,i.y],fontsize = text_size)

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

    if show==True:
        plt.show()
    else:
        plt.savefig(show)
