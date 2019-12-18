import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


#PLOTMATRIX Scatter plot matrix.
#   PLOTMATRIX(X,Y) scatter plots the columns of X against the columns
#   of Y.  If X is P-by-M and Y is P-by-N, PLOTMATRIX will produce a
#   N-by-M matrix of axes. PLOTMATRIX(Y) is the same as PLOTMATRIX(Y,Y)
#   except that the diagonal will be replaced by HIST(Y(:,i)). 
#
#   PLOTMATRIX(...,'LineSpec') uses the given line specification in the
#   string 'LineSpec' '.' is the default (see PLOT for possibilities).  
#
#   PLOTMATRIX(AX,...) uses AX as the BigAx instead of GCA.
#
#   [H,AX,BigAx,P,PAx] = PLOTMATRIX(...) returns a matrix of handles
#   to the objects created in H, a matrix of handles to the individual
#   subaxes in AX, a handle to big (invisible) axes that frame the
#   subaxes in BigAx, a matrix of handles for the histogram plots in
#   P, and a matrix of handles for invisible axes that control the
#   histogram axes scales in PAx.  BigAx is left as the CurrentAxes so
#   that a subsequent TITLE, XLABEL, or YLABEL will be centered with
#   respect to the matrix of axes.
#
#   Example:
#       x = randn(50,3) y = x*[-1 2 12 0 11 -2 3]'
#       plotmatrix(y)
#
#   Clay M. Thompson 10-3-94
#   Copyright 1984-2005 The MathWorks, Inc.
#   $Revision: 1.1 $  $Date: 2008-11-29 12:04:58 $
def plotmatrix(varargin):
    
    # Parse possible Axes input
    cax,args,nargs = axescheck(varargin[:])
    error(nargchk(1, 4, nargs, 'struct'))
    nin = nargs
    
    sym = '.' # Default scatter plot symbol.
    dohist = 0
    
    xx = 1
    linespec = 0
    putlabels = 0
    if nin==1: # plotmatrix(x)
        xx = 1
    elif: nin==2:
        if ischar(args[2]): #plotmatrix(x,line)
            xx = 1
            linespec = 1
            sym = args[2]
            l, c, m, msg = colstyle(sym) ##ok
        if ~isempty(msg):
           error(msg) ##ok
        elif iscellstr(args[2]): # plotmatrix(x,text)
            xx = 1
            putlabels = 1
            varnames = args[2]
        else: # plotmatrix(x,y)
            xx = 0
    elif nin==3:
        linespec = 1
        if (iscellstr(args[3])) & (ischar(args[2])): #plotmatrix(x,line,text)
            xx = 1
            putlabels = 1
            varnames = args[3]
            linespec = 1
            sym = args[2]
            l, c, m, msg = colstyle(sym) ##ok
            if ~isempty(msg)
               error(msg) ##ok
        elif (isfloat(args[2])) & (iscellstr(args[3])): # plotmatrix(x,y,text)
            xx = 0
            putlabels = 1
            varnames = args[3]
        
        else:  # plotmatrix(x,y,line)
            xx = 0
            putlabels = 0
            sym = args[3]
            l, c, m, msg = colstyle(sym) ##ok
            if ~isempty(msg)
               error(msg) ##ok
    else: # plotmatrix(x,y,text,line)
        linespec = 1
        xx = 0
        putlabels = 1
        varnames = args[3]
        sym = args[4]
        l, c, m, msg = colstyle(sym) ##ok
        if ~isempty(msg)
           error(msg) ##ok
    
    #if ischar(linespec):
    #  sym = args[nin]
    #  l, c, m, msg = colstyle(sym) ##ok
    #  if ~isempty(msg)
    #     error(msg) ##ok
    #  nin -=1
    
    if xx: # plotmatrix(y)
      rows = np.shape(args[1], 2)
      cols = rows
      x = args[1]
      y = args[1]
      nvar = np.shape([x, 2])
      dohist = 1
    elif ~xx: # plotmatrix(x,y)
      rows = np.shape(args[2], 2)
      cols = np.shape(args[1], 2)
      x = args[1]
      y = args[2]
      nvar = np.shape([x, 2]) + np.shape([y, 2])
    else:
      print 'MATLAB:plotmatrix:InvalidLineSpec: Invalid marker specification. Type ''help plot''.'
    if (putlabels) & (len(varnames)~=nvar):
        print 'texto sem valores suficientes'
    
    # Don't plot anything if either x or y is empty
    patches = []
    pax = []
    if (isempty(rows)) | (isempty(cols)):
        if nargout>0:
            h = []
            ax = []
            BigAx = []
        return
    
    if (ndims(x)>2) | (ndims(y)>2):
        print 'InvalidXYMatrices: X and Y must be 2-D.'
    if (np.shape([x, 1])~=np.shape([y, 1])) | (np.shape([x, 3])~=np.shape([y, 3])):
        print 'XYSizeMismatch: X and Y must have the same number of rows and pages.'
    
    # Create/find BigAx and make it invisible
    #BigAx = newplot(cax)
    #fig = ancestor(BigAx,'figure')
    #hold_state = ishold(BigAx)
    #set(BigAx,'Visible','off','color','none')
    
    #if any(sym=='.'):
    #    units = get(BigAx,'units')
    #    set(BigAx,'units','pixels')
    #    pos = get(BigAx,'Position')
    #    set(BigAx,'units',units)
    #    markersize = max(1,min(15,round(15*min(pos(3:4))/max(1,size(x,1))/max(rows,cols))))
    #else:
    #    markersize = get(0,'defaultlinemarkersize')
    
    # Create and plot into axes
    ax = np.zeros(rows, cols)
    #pos = get(BigAx,'Position')
    #width = pos(3)/cols
    #height = pos(4)/rows
    #space = .02 # 2 percent space between axes
    #pos(1:2) = pos(1:2) + space*[width height]
    m = np.shape([y, 1])
    k = np.shape([y, 3])
    xlim = np.zeros([rows, cols, 2])
    ylim = np.zeros([rows, cols, 2])
    #BigAxHV = get(BigAx,'HandleVisibility')
    #BigAxParent = get(BigAx,'Parent')
    for i in range(rows: -1: 1):
        for j in range(cols: -1: 1):
            axPos = [pos[1]+[j-1]*width, pos[2]+[rows-i]*height, width*[1-space], height*[1-space]]
            findax = findobj(fig, 'Type', 'axes', 'Position', axPos)
            if isempty(findax):
                ax[i,j] = axes('Position', axPos, 'HandleVisibility', BigAxHV, 'parent', BigAxParent)
                set(ax[i,j], 'visible', 'on')
            else:
                ax[i,j] = findax[1]
            hh[i,j,:] = plot(np.reshape(x[:,j,:],[m,k]), reshape(y[:,i,:],[m,k]), sym, 'parent', ax[i,j])
            set(hh[i,j,:], 'markersize', markersize)
            set(ax[i,j], 'xlimmode', 'auto', 'ylimmode', 'auto', 'xgrid', 'off', 'ygrid', 'off')
            xlim(i,j,:) = get(ax[i,j],'xlim')
            ylim(i,j,:) = get(ax[i,j],'ylim')
    
    if (putlabels) & (~xx):
        count = 1
        # xlabel
        for i in range(np.shape(ax)[2]):
            set(get(ax[end, i], 'xlabel'), 'String', varnames[count])
            count += 1
        # ylabel
        for i in range(np.shape(ax)[1]):
            set(get(ax[i, 1], 'ylabel'), 'String', varnames[count])
            count += 1
    
    #xlimmin = np.min(xlim[:,:,1],[],1)
    #xlimmax = np.max(xlim[:,:,2],[],1)
    #ylimmin = np.min(ylim[:,:,1],[],2)
    #ylimmax = np.max(ylim[:,:,2],[],2)
    xlimmin = np.min(x,[],1)  
    xlimmax = np.max(eps, np.max(x,[],1))
    ylimmin = np.min(y,[],1)
    ylimmax = np.max(eps, np.max(y,[],1))
    
    # Try to be smart about axes limits and labels.  Set all the limits of a
    # row or column to be the same and inset the tick marks by 10 percent.
    inset = 0.15
    for i in range(rows):
        set(ax[i,1], 'ylim', [ylimmin[i,1], ylimmax[i,1]])
        dy = np.diff(get(ax[i,1], 'ylim'))*inset
        set(ax[i,:], 'ylim', [ylimmin[i,1)-dy, ylimmax[i,1]+dy])
    dx = np.zeros([1,cols])
    for j in range(cols):
        set(ax[1,j], 'xlim', [xlimmin[1,j], xlimmax[1,j]])
        dx(j) = np.diff(get(ax[1,j], 'xlim'))*inset 
        set(ax[:,j],'xlim',[xlimmin[1,j]-dx[j] xlimmax[1,j]+dx[j]])
    
    set(ax[1:rows-1, :], 'xticklabel', '')
    set(ax[:, 2:cols],   'yticklabel', '')
    set(BigAx,'XTick', get(ax[rows,1], 'xtick'), 'YTick', get(ax[rows,1], 'ytick'), 'userdata', ax, 'tag', 'PlotMatrixBigAx')
    
    if dohist: # Put a histogram on the diagonal for plotmatrix(y) case
        for i in range(rows:-1:1):
            histax = axes('Position', get(ax[i,i],'Position'), 'HandleVisibility', BigAxHV, 'parent', BigAxParent)
            nn,xx = np.hist(np.reshape(y[:,i,:],[m,k]))
            patches[i,:] = plt.bar(histax, xx, nn, 'hist')
            if putlabels:
                xt = 0.5*(np.max(y[:,i,:])+np.min(y[:,i,:]))
                yt = 0.9*np.max(nn)
                txt = varnames[i]
                plt.text(xt, yt, txt)
            set(histax, 'xtick', [], 'ytick', [], 'xgrid', 'off', 'ygrid', 'off')
            set(histax, 'xlim', [xlimmin[1,i]-dx[i], xlimmax[1,i]+dx[i]])
            pax[i] = histax # ax handles for histograms
      patches = patches
    
    # A bug seems to occur when plotmatrix is ran to produce a plot inside a GUI 
    # whereby the default fig menu items and icons appear. Commenting out the code below fixed the issue.
    # Tim Peterson -  April 2016
    # Make BigAx the CurrentAxes
    # set(fig,'CurrentAx',BigAx)
    # if ~hold_state,
    #    set(fig,'NextPlot','replace')
    # end
    
    # Also set Title and X/YLabel visibility to on and strings to empty
    set([get(BigAx,'Title') get(BigAx,'XLabel') get(BigAx,'YLabel')], 'String','','Visible','on')
     
     if nargout~=0:
       h = hh
     
    def id(strg):
        return 'MATLAB:plotmatrix:',strg
    
    return h, ax, BigAx, patches, pax