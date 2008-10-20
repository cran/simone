## __________________________________________________________
##
## Gplot
##
## INPUT
##	A   : graph precision matrix
##	cl        : node classification vector 
##	cols      : colours to use (per class if cl is specified) 
##	coord	  : matrix of node coordinates for gplot
##	labels  : vector of node labels (to override dimnames(A))
##	degree.threshold : threshold under which degrees are considered null ;
##		null-degree nodes are considered as dutbin nodes
##	display.isolate : boolean : plot dustbin nodes?
##
##	max.edges : max number of edges below which the graph is plotted
##
##	main	  : graphical parameter main
##	sub       : graphical parameter sub
##
## OUTPUT :
##	gplot matrix of node coordinates
##	gplot graphic
##
## Wrapper for gplot in our case :
## Our graphs have real-valued undirected edges, have 
## coloured vertices, and don't have self-loops.
## __________________________________________________________
##
Gplot <- function (A, cl=NULL, ...) {
 
  ## defaults for hidden optional arguments ##
  cols   <- sub.param("cols" , NULL , ...) # colours for each class
  coord  <- sub.param("coord", NULL , ...) # coordinates for nodes
  labels <- sub.param("labels", NULL, ...) # labels for nodes
  degree.threshold <- sub.param("degree.threshold", .Machine$double.xmin, ...)
  display.isolate  <- sub.param("display.isolate", TRUE, ...)
  max.edges <- sub.param("max.edges", 10000, ...) # max. number of plottable edges 
  main <- sub.param("main", NULL, ...) # graph title
  sub  <- sub.param("sub", NULL, ...) # graph subtitle

  if (sum(abs(A))==0) {
    cat("Gplot warning : no edges to plot \n")
    return(NULL)
  }
  
  display.labels <- !is.null(labels)
  
  # dimension check
  p <- dim(A)[1]
  if (p != dim(A)[2]) {
    cat("Gplot warning : A matrix must be square \n")
    return(coord)
  }
  old.p <- p
  
  if (!isSymmetric(A)) {
    cat("Gplot warning : A matrix must be symmetric \n")
    cat("\t automatically symmetrizing using weak rule. \n")
    A <- Symmetrize(A)
  }

  # identify dustbin nodes
  b.dust <- rep(FALSE, p)
  if (!is.null(cl)) {
    b.dust[cl=="D"] <- TRUE
  } else {
    b.dust  <- ( CalcDegrees(abs(A)) <= degree.threshold )
  }

  # gplot handles dustbin nodes poorly
  # -> if display.isolate == FALSE, we'll just drop those nodes!
  if (!display.isolate) {
    A <- A[!b.dust,!b.dust]
    p <- dim(A)[1]
    if (!is.null(coord)) {
      coord     <- coord [!b.dust,]
    }
    if (!is.null(cl)) {
      cl        <- cl [!b.dust]
    }
    if (!is.null(labels)) {
      labels <- labels [!b.dust]
    }
  } 


  # gplot threshold and line weight handling
  e.col <- A
  e.col[A<0]  <- "red"		# red : negative partial correlation (inhibition)
  e.col[A>0]  <- "blue"		# blue : positive partial correlation (activation)
  e.col[A==0] <- "black"	

  lwd <- abs(A)
  lwd <- lwd * 5 / max( lwd )	        # max line weight is 5
	
  e.lty <- lwd
  e.lty[abs(lwd)>1]  <- "solid"
  e.lty[abs(lwd)<=1] <- "dotted"
  e.lty[abs(lwd)==0] <- "blank"


  # get node labels
  if (is.null(labels)) {
   labels <- dimnames(A)[[1]]
  }


  # get node colours
  if (is.null(cl)) {
    if (is.null(cols)) {
      vertex.col <- rep(1, p)
    } else {
      vertex.col <- cols
    }
  } else {
    cl <- as.factor(cl)
    if (is.null(cols)) {
      cols <- rainbow(length(levels(cl)))
    }
    vertex.col <- cols[as.numeric(cl)]
  }
  if (display.isolate) {
    vertex.col[b.dust] <- NA
  }


  # strictly undirected graph, matrix considered symmetric
  A <- abs(A)
  A[lower.tri(A)] <- 0

  # only plot if not too many edges
  num.edges <- sum(A>.Machine$double.xmin)
  if (num.edges>max.edges) {

    cat("Gplot warning : too many edges to plot (",num.edges,">",max.edges,") \n")
    g  <- coord
    g2 <- g

  } else {
    gr <- Gplot.graphics( A,
                          thresh		= 0,	    # no edge blocking
                          margin		= 0.1,
                          main		= main,
                          sub		= sub)
    gr <- Gplot.network( gr,
                         drawloops       = FALSE,
                         displayisolates = display.isolate,
                         coord           = coord)
    gr <- Gplot.edge( gr,
                      col       = e.col,
                      lty       = e.lty,
                      lwd       = lwd,
                      usearrows = FALSE)
    gr <- Gplot.vertex( gr, col= vertex.col )

    if ( display.labels )
      gr <- Gplot.vertex.label( gr,
                                label=labels,
                                cex=0.7,
                                useboxes=FALSE
                               )
    if (display.isolate) {
      g2 <- gr$network$coord
    } else {
      g2 <- matrix(0, old.p, 2)
      g2[!b.dust,] <- gr$network$coord
    }
  }

  invisible(g2)

}

## __________________________________________________________
##
## Mplot:
##        Display an image of a matrix A.
##
## Authors : Alexander Smith and Gilles Grasseau
##
## __________________________________________________________

# Defining "grey" and "light" palette
grey.palette <- function( n, ascending ) {
  
  if( n > 1 ) {
    if ( ascending >= 0 ) 
      ret <- grey( 0:(n-1) / (n-1) )
    else
      ret <- grey( (n-1):0 / (n-1) )
  } else {
    ret <- grey(0)
  }
  return(ret)
}

light.palette <- function( n ){

  R1 <- vector()
  G1 <- vector()
  B1 <- vector()

  R1 <- c( 0    , 0.08 , 0.15 ,
          0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    ,
          0    , 0.376, 0.6  , 0.8  , 0.9  , 0.95 , 1    , 1    , 1    ,
          1    , 1    , 1    , 1
        )
  G1 <- c( 0    , 0    , 0    ,
          0    , 0.294, 0.498, 0.725, 0.933, 1    , 1    , 1    ,
          1    , 1    , 1    , 1    , 1    , 1    , 1    , 0.9  ,0.796,
          0.7  , 0.569, 0.365,
          0.0
        )

  B1 <- c( 0    , 0.5  , 0.750,
          1    , 1    , 1    , 1    , 1    , 0.863, 0.659, 0.355,
          0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    ,
          0    , 0    , 0    , 0
        )
  if( n > 1 ) {
    
    X1 <- 1:length( R1 )
    X <- (0:(n-1)) * (length(R1)-1) / (n-1) + 1

    R <- approx(X1, R1, X )
    G <- approx(X1, G1, X )
    B <- approx(X1, B1, X )
  } else {
    return( rgb(0,0,0) )
  }
  
  return( rgb(R$y, G$y, B$y) )
}

                                        #
#  Environment "Mplot.env"
#       used to store some global data
#  --------------------------------------
Mplot.env <- new.env( parent=baseenv() )
assign("Mplot.CLegend", list( min=0, max=1, col=light.palette(32) ),
           envir=Mplot.env )
                              
## __________________________________________________________
##
## Mplot
##
##
##  Display a image of a matrix A. Each value A(i,j) corresponds
##  to a color contained in the given color table called here palette.
##  Offers the possibility to permute rows and colomns according to
##  the classification argument "cl".
##
##  Usage:
##
##     Mplot( A,  cl = NULL, A.diag = TRUE, cl.order = TRUE)
##
##  Arguments:
##
##           A: Matrix 'A' to display. 'NA' values are colored with the
##              background color 'bg'.
##
##      A.diag: Logical value. If 'TRUE', diagonal terms are set to zero. 
##
##          cl: Vector of classes identifier (integers). 'A' must be a square
##              matrix and 'cl' must be a vector of the same dimension.
##
##    cl.order: Logical value. It true, the matrix is reordered to display
##              classes in class size ascending order. 
##
##         ...: see "details" section. 
##
##  Details:
##     Graphical parameters can be specified in the call to 'Mplot' to
##     change the color palette, the axis labels, and general parameters. 
##
##
##     * Palette parameters
##
##        'colors' (default "light.64") Defines the palette used when
##             drawing the matrix. 'rainbow', 'heat.colors',
##             'terrain.colors', 'topo.colors', 'cm.colors' are available
##             'R' (the numer of colors is set with the mandatory
##             argument of these functions - e.g. 'rainbow(10)' generates
##             10 colors from the  'rainbow' palette). A grey palette
##             'grey' and a light spectrum-like palette 'light' are
##             provided in this function. The number of color levels are
##             chosen by adding ".xxx" where 'xxx' is the level number
##             (e.g. '"light.256"' ). 
##
##        'normalize'(default 'TRUE') 'A' values are normalized  so
##             that generated colors belong to interval
##             '1:length(colors)'. If 'normalize = FALSE' then the
##             normalization of the previous 'Mplot' is used. In this
##             case, values out of bounds are set to colors[1] (if the
##             color is less than one) or 'colors[length(colors)]' (if
##             the color is greater than 'length(colors)'). 
##
##        'color.legend' (default 'FALSE') Display the color map.
##
##
##     * Label parameters
##
##        Axis labels are set by using 'rownames()' and 'colnames()'
##        functions on the initial matrix. Default values are the row
##         and column matrix indexes. Axis label parameters:
##
##
##        'lab.col.cex' (default 1.0) Expansion value of the column
##             labels (y axis).
##
##        'lab.row.cex' (default 1.0) Expansion value of the row
##             labels (x axis).
##
##        'lab.col.pos' (default "h") Label orientation of the
##             columns "v" for vertical or '"h"' for horizontal).
##
##        'lab.row.pos' (default "h") Label orientation of the rows
##             ("v" for vertical or "h" for horizontal). 
##
##     * Other parameters: 
##
##        'margin' (default 0.0) Margin between the matrix
##             representation and the axis.
##
##        'new' (default 'FALSE') Draw in a new window.
##
##        'main' (default "") Main title.
##
##        'sub' (default "") Subtitle.
##
##        ... Other arguments pass to functions 'image' and 'axis'. See
##             these functions for more.
##
##
## Value:
##     No value is returned.
##
## Author(s):
##
##     A. Smith and G. Grasseau
##
##
## __________________________________________________________
##
Mplot <- function (A=NULL, cl=NULL, ...) {

  
  # Optionnal arguments ...
  # -----------------------
  colors       <- sub.param("colors"      , "light.64", ...)
  normalize    <- sub.param("normalize"   , TRUE, ...)
  lab.display  <- sub.param("lab.display" , FALSE, ...)
  lab.col.cex  <- sub.param("lab.col.cex" , 1.0, ...)     
  lab.row.cex  <- sub.param("lab.row.cex" , 1.0, ...)     
  lab.col.pos  <- sub.param("lab.col.pos" , "h", ...)    
  lab.row.pos  <- sub.param("lab.row.pos" , "h", ...)    
  color.legend <- sub.param("color.legend", FALSE, ...)    
  margin       <- sub.param("margin"      , 0.0, ...)    
  main         <- sub.param("main"        , "", ...)    
  sub          <- sub.param("sub"         , "", ...)    
  A.diag       <- sub.param("A.diag"      , FALSE, ...)    
  cl.order     <- sub.param("cl.order"    , TRUE, ...)   

  # Test arguments and configure the function 
  # -----------------------------------------

  null.matrix   <- is.null(A)
  square.matrix <- FALSE
  if( !null.matrix )
    square.matrix <- (dim(A)[1] == dim(A)[2])
  
  # Verifying arguments
  if( ! is.null(cl) && is.null(A) ) {
    if( ! square.matrix || (length(cl) != dim(A)[1]) ) {
      cat(">> Warning: 'cl' and 'A' dimensions differ.\n")
      cat(">>          Classes vector 'cl' is ignored.\n")
      cl <- NULL
    }
  }
  
  lab.row.angle <- 0
  if ( lab.row.pos == "v")
    lab.row.angle <- 1

  lab.col.angle <- 1
  if ( lab.col.pos == "h")
    lab.col.angle <- 0
    
  layout.mat <- matrix( c(1), 1, 1 )
  layout.w   <- c(1)
  
  # Prepare two boxes for the color table
  if ( color.legend && ! (null.matrix) ) {
    layout.mat <- matrix( c(1,2), 1, 2 )
    layout.w   <- c(3,1)
    layout( layout.mat, widths  = layout.w, heights = 1)
  }

  # Colors Table
  if ( is.null(colors) )  {
    colors <- light.palette(32)
    if( ! is.null( cl ))
      colors <- light.palette(length( cl) )
    else
      colors <- light.palette( 32 )
  }
  else {
    if( length(colors) != length( grep("^#", colors ) ) ){
      # Not a Color table
      string <- strsplit(colors, '[.]')
        
      # Level number
      if ( length( grep("[^0-9]", string[[1]][2] ) == 1))
        level.nbr <- 32
      else 
        level.nbr <- as.numeric( string[[1]][2] )
    
      # Palette name
      if( string[[1]][1] == "inv-grey" )
        colors <- grey.palette( level.nbr, -1 )
      else if( string[[1]][1] == "grey" )
        colors <- grey.palette( level.nbr, 1 )
      else 
        colors <- light.palette( level.nbr )
    }
  }
  number.colors <- length( colors )
  
  # Matrix Diagonal terms
  if (A.diag==FALSE && square.matrix) {
    diag(A) <- 0 # drop diagonal
  }

  if( ! null.matrix ) {
     
    # Ticks position         
    x.axis <- 1:dim(A)[2]
    y.axis <- 1:dim(A)[1]
    n.x <- dim(A)[2]
    n.y <- dim(A)[1]
  
    # Limits
    # Margin + shift to enterely draw the sided boxes 
    xlim<-c( (-margin)*dim(A)[2] +0.5, (1+margin)*dim(A)[2] +0.5)
    ylim<-c( (-margin)*dim(A)[1] +0.5, (1+margin)*dim(A)[1] +0.5)

    # Classes  treatement
    # -------------------
    o.x <- 1:dim(A)[2]
    o.y <- 1:dim(A)[1]
    if ( !is.null(cl) ) {
      # order matrix A versus classes
      cl  <- as.factor(cl)

      if (cl.order) {
        # Indexes (l.o) which sort cl by the class frequency 
        l.o <- order(summary(cl))
        cl  <- factor(cl, levels=levels(cl)[l.o])
      }

      o.x <- order(cl)
      o.y <- o.x
    }

    # To set a fixed size of the box(es)
    w <- par()$pin[1] 
    h <- par()$pin[2]

    
    # Normalize colors
    if( normalize ) {

      A.min <- min( A, na.rm = TRUE )
      A.max <- max( A, na.rm = TRUE )
      
    } else {
      # Restore Environment
      Y <- get("Mplot.CLegend", envir=Mplot.env)
      A.min  <- Y$min
      A.max  <- Y$max
      A[ which( A > A.max ) ] = A.max
      A[ which( A < A.min ) ] = A.min
    }

    image( x.axis, y.axis, t(A[o.y,o.x]),
           col=colors, main=main,sub=sub,
           xlim=xlim, ylim=ylim, zlim=c(A.min, A.max),
           xlab="", ylab="",
           axes=FALSE)

    box()
  
    #      Draw class limits
    #      ------------------
  
    if ( ! is.null( cl ) ) {
    
      # Draw classe limits
      width.case <- 1
      retenue    <- -(1/2) * width.case
      
      for (i in 1:(nlevels(cl)-1)) {
        q <- levels(cl)[i]
        prop <- sum ( 1*(cl==q) ) 
        pq <- retenue + prop * width.case
        if (prop>0) {
          abline(v=pq, lty="dashed")
          abline(h=pq, lty="dashed")
        }
        retenue <- pq
      }   
    }

    #      Calculate Tick and Label position
    #      --------------------------------

    if( lab.display ) {
      
      if(is.null( colnames(A) )) {
        colnames(A) <- 1:dim(A)[2]
      }
      if(is.null( rownames(A) )) {
        rownames(A) <- 1:dim(A)[1]
      }
      
      # I indexes
      if( length( grep("[^0-9]", rownames(A) ) > 0) ) {
        # Character case
        char.case <- 1
      } else {
        # Numeric case
        char.case <- 0
      }
      ordered.names <-  rownames(A)[o.y]
      i.list <- get.label.index.list( y.axis, ordered.names,
                                     char.case, lab.row.cex,
                                     as.numeric(!lab.row.angle), 1 ) 

      axis( side=2, at=y.axis[ i.list ], labels= ordered.names[i.list],
            cex.axis=lab.row.cex, las=lab.row.angle*2+1 )
      
      # J indexes
      if( length( grep("[^0-9]", colnames(A) ) > 0) ) {
        # Character case
        char.case <- 1
      } else {
        # Numeric case
        char.case <- 0
      }

      ordered.names <-  colnames(A)[o.x]
      i.list <- get.label.index.list( x.axis, ordered.names,
                                   char.case, lab.col.cex, lab.col.angle, 0) 

      axis( side=1, at=x.axis[ i.list ], labels=ordered.names[i.list],
            cex.axis=lab.col.cex, las=lab.col.angle*2+1 )
    }

    #     Save context for the color table legend
    #     ---------------------------------------
    if( ! null.matrix ) {

      assign("Mplot.CLegend", list( min=A.min, max=A.max, col=colors ),
             envir=Mplot.env )
    }
  }
    
  #    Legend of Color table
  #    ---------------------
  if ( color.legend ) {

    legend.labels <- TRUE
    w  <- par()$pin[1] 
    h  <- par()$pin[2]
    wn <- w
    # Convert inches to cm ( * 2.54)
    if( w*2.54 > 1.0)
      wn = 1.0 / 2.54
    
    par( pin = c( wn, h) )
    
    # Restore Environment
    Y <- get("Mplot.CLegend", envir=Mplot.env)
    A.min  <- Y$min
    A.max  <- Y$max
    colors <- Y$col
 
    legend <- vector()
  
    legend[2:(number.colors+1)] <-  A.min +
      (A.max - A.min) * ((0:(number.colors-1)) + 1.0 ) / number.colors 
    legend[1] <- NA

    image( c(1), 0:(number.colors), matrix(legend, 1, number.colors+1),
           col=colors, main="Color legend",
           xlim=c(0,1), ylim=c( -0.5, number.colors+1-0.5),
           xlab="", ylab="",axes=FALSE)
    box()

    if ( legend.labels ) {
      pos <- 0:number.colors + 0.5
      labels <- vector()
      for ( i in 0:number.colors ) {
        labels[i+1] <- sprintf("%6.3g",
                       A.min + (A.max - A.min) * i /  number.colors )
      }
      char.case <- 1
      i.list <- get.label.index.list( pos, labels, char.case, lab.row.cex, 1, 1 )
      axis( side=2, at=pos[ i.list ], labels=labels[i.list],
            las=1, cex.axis=lab.col.cex)
    
      axis( side=4, at=0.0, labels="NA",tick = FALSE, las=1,
            cex.axis=lab.col.cex*0.8, line=-0.8 )
    }
    
    # Restore graphics context
    par( pin = c( w, h) )
  }
}    

get.label.index.list <- function( x.axis, labels, char.case, cex, rot, axis ){ 

  # rot = 0 - the labels is colinear with the axis direction
  # rot = 1 - label perpendicular to the axis 
  # axis = 0 - X axis.
  # axis = 1 - Y axis.
  n.x  <- length( x.axis )
  if( n.x !=  length( labels ) )
     cat("ERROR - get.label.index.list: bad arguments \n")
  
  if( char.case == 1) {
    # Character case
    x.start   <- 0
    x.end   <- n.x-1
  } else {
    # Numeric case
    x.start <- x.axis[1]
    x.end   <- x.axis[n.x]
  }

  par(cex.lab = cex)

  # Label margin
  l.margin.w <- 0.5
  l.margin.h <- 0.2
    
  # Label Max size (adding a label margin )
  if( rot == 1 ) {
    # labels and axis are perpendicular
    l.max.size <- ( max( strheight( labels , cex = cex) )
                   + 2*l.margin.h*par()$cxy[2] )
  } else {
    # labels and axis have the same direction
    l.max.size <- ( max( strwidth( labels , cex = cex) )
                   + 2*l.margin.w*par()$cxy[1] )
  }
  if( (axis == 1) && (rot == 0)) {
    # Y axis
    #
    # Tranform X length to Y length
    usr.coord <- par()$usr
    l.max.size <- l.max.size *
      (usr.coord[4] - usr.coord[3]) / (usr.coord[2] - usr.coord[1])
  }
  if( (axis == 0) && (rot == 1) ) {
    # X axis
    #
    # Tranform Y length to X length
    usr.coord <- par()$usr
    l.max.size <- l.max.size *
      (usr.coord[2] - usr.coord[1]) / (usr.coord[4] - usr.coord[3])
  }

  # Number max of labels on the axis
  t.max.labels <- (x.end - x.start)/ l.max.size
  
  t.space <- ticks.space( x.end - x.start, t.max.labels ) 
  t.start <- ceiling( x.start / t.space ) * t.space
  
  # x.axis are supposed to be consecutive (i,i+1, .., n-1,n)
  i.list  <- seq(t.start+char.case, n.x, t.space)
  return( i.list )
}

ticks.space <- function( interval, nb.of.label ) 
{
  val <- interval / nb.of.label 
  exp <- trunc( log10( val ) )
 
  # Mantisse
  ratio <- val / (10^exp)
  fact <- switch( ceiling(ratio) , 1, 2, 4, 4, 5, 10, 10, 10, 10, 10 )
  ret <- fact * 10^(exp)
 
  return( ret )
}


##########################################################################
#  
#  Routines to draw graphs.
#
#  These routines have been extracted from the SNA package (mainly gplot)
#  and reorganized to exploit more conveniently its possibilies.
#
#  Author : G.Grasseau and A.Smith (Statistics and Genome laboratory)
#
#  Entry points:
#  ------------
#
#    Gplot.graphics:     set default graphic parameters.
#                        Main steps:
#                        + store initial matrix as boolean one,
#                        + build the structure (list) which handle all
#                        graphics parameters.
#
#    Gplot.network:      specify the network representation (up to now only
#                        undirect graphs are taken into account).
#                        Main steps:
#                        + compute vertex position,
#                        + graph box limits and scaling factor,
#                        + define which vertices to display.
#
#    Gplot.vertex:       draw vertices.
#
#    Gplot.vertex.label: draw vertex labels. 
#
#    Gplot.edge:         draw edges
#                        Main steps:
#                        + identify edges to be drawn.
#                        + remove loops
#                        + draw edges/arrows.
#
#  Utilities:
#  ---------
#   + draw.edge:        (used by Gplot.edge) draw edges/arrows
#   + isolate.vertices: (used by Gplot.network) tag isolate vertices (TRUE)
#
#  To improve:
#  ----------
#   + sides of some boxed labels are not displayed
#   + 2 variables (displayisolates and use.isolates) should be introduced
#     to avoid taking into account of isolate vertices in the vertex coordinate
#     computation (4 cases have to be studied). 
#
#  To implement:
#  ------------
#   + drawing the loops
#
##########################################################################

Gplot.graphics<-function( mat, thresh=0, xlim=NULL, ylim=NULL, scale=0.01,
                          margin=0.2, main="", sub=""
                        ) {
# -------- Arguments -----------------------------------------------------
#
# mat        (matrix): initial matrix
# thresh     (scalar): matrix values are set to zero if they are < threshold
#                      (Pb: should be "abs( mat ) < 0")
# xlim, ylim (vector): x minimun and x maximum (same for y)
# scale  (scalar)    : scaling factor for the whole graphics.
# margin (scalar)    : margin value to add to all the graphics box sides.
# main  (char): title  
# sub   (char): subtitle
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters (default)
#
# ------------------------------------------------------------------------
  
   n<-dim(mat)[1]

   # Replace NAs with 0s
   mat[is.na(mat)]<-0
   
   # Save a copy of mat
   mat.raw<-mat
   
   # Binary matrix
   mat<-matrix(as.numeric(mat>thresh),n,n)

   l.graphics <- list( xlim=xlim, ylim=ylim, scale=scale, margin=margin,
                       main=main, sub=sub ) 

   l.network <- list( mode=NULL, drawloops=FALSE, vertex.pos.mode=NULL,
                      coord= NULL, displayisolates=TRUE, use=NULL )
   
   l.label   <- list( label=c(1:dim(mat)[1]), cex=1, col=1, pos=0,
                      useboxes=TRUE, box.margin=0.5, box.col=1,
                      box.bg="white", box.lty=NULL, box.lwd=par("lwd") )

   l.vertex <- list( cex=1, sides=8, col=2, border=1, lty=1 ,
                     label = l.label, baserad=0 )

   l.edge  <- list( col=1, lty=1, lwd=0,
                    usearrows=TRUE, arrow.cex=1, loop.cex=1 )

   graph   <- list( mat=mat, thresh=thresh, mat.raw = mat.raw,
                    graphics=l.graphics, network=l.network,
                    vertex=l.vertex, edge=l.edge )

   return( graph )
}


Gplot.network <-function( graph,
                          mode=NULL, drawloops=FALSE, vertex.pos.mode=NULL,
                          coord= NULL,
                          displayisolates=TRUE, ...)
{
# -------- Arguments -----------------------------------------------------
#
# graph (list): structure handling all graphics parameters
# mode  (char): kind of network to draw. 
# drawloops       (bool): draw network loops (not implemented).
# vertex.pos.mode (char): vertex positionning mode (not used).
# coord         (vector): vertex position.
#                         If NULL these coordinates are computed
# displayisolates (bool): display isolate vertices.
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters
#
# ------------------------------------------------------------------------

  # Arguments mode and vertex.pos.mode are not implemented
  if( ! missing(drawloops) )
    graph$network$drawloops = drawloops
  if( ! missing(coord) )
    graph$network$coord = coord
  if( ! missing(displayisolates) )
    graph$network$displayisolates = displayisolates

  if(is.null(graph$network$coord)){      
    
    # Provide default settings
    n           <- dim(graph$mat)[1]
    niter       <- 500
    max.delta   <- n
    area        <- n^2
    cool.exp    <- 3
    repulse.rad <- area*n
    
    # Set initial positions randomly on the circle
    tempa<-sample((0:(n-1))/n)
    x<-n/(2*pi)*sin(2*pi*tempa)
    y<-n/(2*pi)*cos(2*pi*tempa)

    layout<-.C( "vertex_coord_C",
                as.integer(graph$mat),
                as.double(n), as.integer(niter),
                as.double(max.delta), as.double(area),
                as.double(cool.exp), as.double(repulse.rad),
                x=as.double(x), y=as.double(y),
                PACKAGE="simone"
               )

    graph$network$coord <- cbind(layout$x,layout$y)
  }
  
  # Remove isolated vertex (if displayisolates FALSE)
  use <- displayisolates | ( ! isolate.vertices( graph$mat ) ) 

  graph$network$use = use
  
  x <- graph$network$coord[,1]
  y <- graph$network$coord[,2]

                         
  # Set limits for plotting region
  xlim = graph$graphics$xlim                
  ylim = graph$graphics$ylim
  margin = graph$graphics$margin
  if(is.null(xlim))
    xlim<-c(min(x[use])-margin,max(x[use])+margin)  # Save x, y limits
  if(is.null(ylim))
    ylim<-c(min(y[use])-margin,max(y[use])+margin)
  xrng<-diff(xlim)          
  yrng<-diff(ylim)
  xctr<-(xlim[2]+xlim[1])/2                 # Get center of plotting region
  yctr<-(ylim[2]+ylim[1])/2
  
  # Force scale to be symmetric
  if(xrng<yrng)
    xlim<-c(xctr-yrng/2,xctr+yrng/2)
  else
    ylim<-c(yctr-xrng/2,yctr+xrng/2)

  graph$graphics$xlim = xlim                
  graph$graphics$ylim = ylim               

  # Extract "base radius"
  graph$vertex$baserad <- min(diff(xlim),diff(ylim))* graph$graphics$scale
  
  # Configure the graphic box
  plot( 0,0,
        xlim=graph$graphics$xlim,
        ylim=graph$graphics$ylim, type="n", xlab="",ylab="", asp=1, axes=FALSE,
        main= graph$graphics$main, sub=graph$graphics$sub,
        ...
       )
  return( graph )
}


Gplot.vertex <-function( graph,
                         cex, sides, col, border, lty, ...)
{
# -------- Arguments -----------------------------------------------------
#
# graph           (list): structure handling all graphics parameters
# cex    (scalar/vector): scaling factors.
# sides  (scalar/vector): number of node sides.
# col    (scalar/vector): vertex colors.
# border (scalar/vector): color of node (vertex) borders.
# lty    (scalar/vector): type of nodes borders.
#                         (Pb: would be better with the line width
#                         border parameter "lwd".
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters 
#
# ------------------------------------------------------------------------

  if( ! missing(cex) )
   graph$vertex$cex    = cex
  if( ! missing(sides) )
   graph$vertex$sides  = sides
  if( ! missing(col) )
   graph$vertex$col    = col
  if( ! missing(border) )
   graph$vertex$border = border
  if( ! missing(lty) )
   graph$vertex$lty    = lty

   n <- dim(graph$mat)[1]
  
   # Build vectors describing vertex
   v.cex    <- rep( graph$vertex$cex                       , length=n)
   v.radius <- graph$vertex$baserad * v.cex
   v.sides  <- rep( graph$vertex$sides                     , length=n)
   v.col    <- rep( graph$vertex$col                       , length=n)
   v.border <- rep( graph$vertex$border                    , length=n)
   v.lty    <- rep( graph$vertex$lty                       , length=n)

   # remove unused
   use = graph$network$use
   v.radius <-  v.radius[use]
   v.sides  <-  v.sides[use]
   v.col    <-  v.col[use]
   v.border <-  v.border[use]
   v.lty    <-  v.lty[use]

   x <- graph$network$coord[use,1]
   y <- graph$network$coord[use,2]   
   n <- length(x)

 # Compute the coordinates
  coord<-vector()

  for(i in 1:n){
    ang <- (1:v.sides[i])/v.sides[i]*2*pi
    dx <- v.radius[i]*cos(ang)
    dy <- v.radius[i]*sin(ang)
    XY = rbind( cbind( x[i]+dx, y[i]+dy ), c(NA,NA) )
    coord<-rbind(coord, XY)
  }
  # Plot the vertices
  polygon(coord, col=v.col, border=v.border, lty=v.lty, ...)

  return( graph )
}


Gplot.vertex.label <-function( graph,
                               label, cex, col, pos,
                               useboxes, box.margin, box.col, box.bg,
                               box.lty, box.lwd,
                               ... )
{
# -------- Arguments -----------------------------------------------------
#
# graph           (list): structure handling all graphics parameters
# label         (vector): label titles (if NULL a default is provided)
# cex    (scalar/vector): scaling factors.
# col    (scalar/vector): label colors.
# pos           (scalar): label positionning mode
#                         + 0 labels are placed away from the graph
#                         + 1 labels are placed below the vertices      
#                         + 2 labels are placed on the vertex left.      
#                         + 3 labels are placed above the vertices      
#                         + 4 labels are placed on the vertex right.
# useboxes        (scalar): frame (box) the labels
# box.margin      (scalar): margin between the label titles and their boxes
#                           (in character size unit).
# box.col  (scalar/vector): box colors.
# box.bg   (scalar/vector): box backgroung color.
# box.lty  (scalar/vector): boxe line type.
# box.lwd  (scalar/vector): boxe line width.
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters
#
# ------------------------------------------------------------------------
  if( ! missing( label ) )
    if( ! is.null( label ) )
      graph$vertex$label$label = label
  if( ! missing( cex ) )
   graph$vertex$label$cex   = cex
  if( ! missing( col ) )
   graph$vertex$label$col   = col
  if( ! missing( pos ) )
   graph$vertex$label$pos   = pos
  if( ! missing( useboxes ) )
   graph$vertex$label$useboxes   = useboxes
  if( ! missing( box.margin ) )
   graph$vertex$label$box.margin = box.margin
  if( ! missing( box.col ) )
   graph$vertex$label$box.col    = box.col
  if( ! missing( box.bg ) )
   graph$vertex$label$box.bg     = box.bg
  if( ! missing( box.lty ) )
   graph$vertex$label$box.lty    = box.lty
  if( ! missing( box.lwd ) )
   graph$vertex$label$box.lwd    = box.lwd

  # Plot vertex labels
  use  <- graph$network$use
  x <- graph$network$coord[use,1]
  y <- graph$network$coord[use,2]
   
  if((!all(graph$vertex$label$label==""))&(!all(use==FALSE))){

    # Label display mode
    if ( graph$vertex$label$pos == 0 ){
      
      # Labels are placed away from the graph
      xoff <- x - mean(x)
      yoff <- y - mean(y)
      roff <- sqrt(xoff^2+yoff^2)
      xhat <- xoff/roff
      yhat <- yoff/roff

    } else if (graph$vertex$label$pos<5) {
      
      # below (0,-1), left (-1,0), top (0,1) , right (1,0)
      xhat <- switch( graph$vertex$label$pos,  0,-1, 0, 1)
      yhat <- switch( graph$vertex$label$pos, -1, 0, 1, 0)

    } else {
      xhat <- 0
      yhat <- 0
    }
     
    # Get character size
    l.cex = graph$vertex$label$cex
    char.len <- par()$cxy * l.cex

    # Get label width and height
    label <- graph$vertex$label$label[use]
    lw <- strwidth( label,cex=l.cex) / 2
    lh <- strheight(label,cex=l.cex) / 2
    b.size   = 1 + graph$vertex$label$box.margin
    
    v.radius <- graph$vertex$baserad * rep( graph$vertex$cex,
                                            dim(graph$mat)[1] )

    # Draw boxes
    if( graph$vertex$label$useboxes ){
      rect(x - lw*b.size + xhat*(lw*(b.size+0.2) + v.radius),
           y - lh*b.size + yhat*(lh*(b.size+0.2) + v.radius),
           x + lw*b.size + xhat*(lw*(b.size+0.2) + v.radius),
           y + lh*b.size + yhat*(lh*(b.size+0.2) + v.radius),
           col    = graph$vertex$label$box.bg,
           border = graph$vertex$label$box.border,
           lty    = graph$vertex$label$box.lty,
           lwd    = graph$vertex$label$box.lwd)
    }

    # Draw labels
    text(x + xhat * ( lw*(b.size+0.2) + v.radius ),
         y + yhat * ( lh*(b.size+0.2) + v.radius ),
         label, cex=l.cex, col=graph$vertex$label.col, offset=0, ...)
    
  }
  return ( graph )
}

Gplot.edge <-function( graph,
                       col=1, lty=1, lwd=0,
                       usearrows=TRUE, arrow.cex=1, loop.cex=1,
                       ... )
{
# -------- Arguments -----------------------------------------------------
#
# graph              (list): structure handling all graphics parameters
# col  (scalar/vector/matrix): edge colors.
# lty  (scalar/vector/matrix): edge line types.
# lwd  (scalar/vector/matrix): edge line widths.
# usearrows          (bool): draw arrows.
# arrow.cex (scalar/vector): arrow head scaling factors.
# loop.cex  (scalar/vector): loop scaling factors.
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters
#
# ------------------------------------------------------------------------
  
   # Remark : arrow.cex and loop.cex could be fused

   if ( ! missing(col) )
     graph$edge$col = col
   if ( ! missing(lty) )
     graph$edge$lty = lty
   if ( ! missing(lwd) )
     graph$edge$lwd = lwd
   if ( ! missing(usearrows) )
     graph$edge$usearrows = usearrows
   if ( ! missing(arrow.cex) )
     graph$edge$arrow.cex = arrow.cex
   if ( ! missing(loop.cex) )
     graph$edge$loop.cex  = loop.cex

   n <- dim( graph$mat )[1]

   # Build vectors describing edges
   # Each edge is a polygon
   px0<-vector()   
   py0<-vector()
   px1<-vector()
   py1<-vector()
   
   e.lwd<-vector()  # Create edge attribute vectors
   e.type<-vector()
   e.col<-vector()
   e.hoff<-vector() # Offset radii for heads
   e.toff<-vector() # Offset radii for tails
   e.diag<-vector() # Indicator for self-ties

   # Coerce edge.col/edge.lty to array form
   if(!is.array(graph$edge$col))   
     col <- array(graph$edge$col, dim=dim(graph$mat))
   else
     col <- graph$edge$col
   if(!is.array(graph$edge$lty))
     lty<-array(graph$edge$lty, dim=dim(graph$mat))
   else
     lty = graph$edge$lty
   if(!is.array( graph$edge$lwd)){
     if( graph$edge$lwd>0 )
       lwd<-array( graph$edge$lwd * graph$mat.raw, dim=dim(graph$mat))
     else
       lwd<-array(1, dim=dim(graph$mat))
   }
   
   v.radius <- graph$vertex$baserad * rep( graph$vertex$cex, dim(graph$mat)[1] )

   # Select edges between vertices
   # -----------------------------
   x <- graph$network$coord[,1]
   y <- graph$network$coord[,2]

   for(i in (1:n)[graph$network$use]) {   
     for(j in (1:n)[graph$network$use]) {
       if( graph$mat[i,j] ){        # Edge exists
         px0 <- c(px0,as.real(x[i]))  # Store endpoint coordinates
         py0 <- c(py0,as.real(y[i]))
         px1 <- c(px1,as.real(x[j]))
         py1 <- c(py1,as.real(y[j]))
         e.toff <-c ( e.toff, v.radius[i] ) # Store endpoint offsets
         e.hoff <-c ( e.hoff, v.radius[j] )
         e.col  <-c ( e.col , col[i,j])     # Store other edge attributes
         e.type <-c ( e.type, lty[i,j])
         e.lwd  <-c ( e.lwd , lwd[i,j])
         e.diag <-c ( e.diag, i==j)         # Set to true if diagonal 
       }
     }
   }

   # Remove loops
   # ------------
   if(length(px0)>0){
     px0 <- px0[!e.diag] 
     py0 <- py0[!e.diag]
     px1 <- px1[!e.diag]
     py1 <- py1[!e.diag]
     e.lwd  <- e.lwd[!e.diag]
     e.type <- e.type[!e.diag]
     e.col  <- e.col[!e.diag]
     e.hoff <- e.hoff[!e.diag]
     e.toff <- e.toff[!e.diag]
   }

   # Draw edges
   if(length(px0)>0)
     draw.edges(
            as.vector(px0),as.vector(py0),as.vector(px1),as.vector(py1),
            width=e.lwd*graph$vertex$baserad/10, col=e.col, lty=e.type,
            o.head=e.hoff, o.tail=e.toff,
            arrow=usearrows, a.len=2*graph$vertex$baserad*arrow.cex, a.angle=20,
            ...
          )

  return( graph )
}
  
isolate.vertices<-function( mat. ) {
# -------- Arguments -----------------------------------------------------
#
# mat. (matrix): binary (0/1) square matrix
#
# -------- Return value --------------------------------------------------
#
# isolate (vector): isolated (if TRUE) vertex vector.
#
# ------------------------------------------------------------------------
  
  n <- dim(mat.)[1];
  mat <- mat.
  if ( n > 1 ){
    # Set to zero NA and diagonal terms
    for(i in 1:n){
      mat[i,i] = 0
      mat[i,1:n] ==  as.numeric( ! is.na(mat[i,1:n]) )
    }
    isolate <- vector()
    for(i in 1:n) {
      isolate = c( isolate, all(( mat[i,] == 0 )) & all(( mat[,i] == 0 )) )
    }
  }
  return( isolate )
}

draw.edges<-function( x0, y0, x1, y1,
                      width=0.01, col=1, lty=1,
                      o.head=0, o.tail=0,
                      arrow=TRUE, a.len=0.1, a.angle=2,
                      ... )
{
# -------- Arguments -----------------------------------------------------
#
# x0, y0       (vector): start coordinates of edges to draw.
# x1, y1       (vector): end coordinates of edges to draw.
# width (scalar/vector): edge line widths.
# col   (scalar/vector): edge colors.
# lty   (scalar/vector): edge line types.
# o.head (scalar/vector): offset (vertex size shift) at the start points.
# o.tail (scalar/vector): offset (vertex size shift) at the end points.
# arrow           (bool): arrows are drawn.
# a.len   (scalar/vector): arrow head lengths.
# a.angle (scalar/vector): arrow head angle (in degree).
#  
# -------- Return value --------------------------------------------------
#
# No value
#
# ------------------------------------------------------------------------

  if(length(x0)==0)   #Leave if there's nothing to do
    return;

  n<-length(x0)

  # Transform scalars into vectors
  width <- rep(width,length=n)
  col   <- rep(col,length=n)
  lty   <- rep(lty,length=n)
  
  # Offsets
  o.head  <- rep(o.head,length=n) 
  o.tail  <- rep(o.tail,length=n)

  # Arrow parameters
  a.angle <- rep(a.angle,length=n)/360*2*pi
  a.len   <- rep(a.len,length=n)

  # Debug point
  # cat("xy  :",x0, y0, x1,y1, "\n")
  # cat("width :", width, "\n")
  # cat("col :", col, "\n")
  # cat("lty :", lty, "\n")
  # cat("o.tail :",o.tail, "\n")
  # cat("o.head :", o.head , "\n")
  # cat("a.angle :", a.angle, "\n")
  # cat("a.len   :", a.len, "\n")

  # Computes edges/arrows coordinates
  coord<-vector()
  XY <- vector()
  for(i in 1:n) {  

    slen<-sqrt((x0[i]-x1[i])^2+(y0[i]-y1[i])^2)  #Find the total length

    if(arrow){
      
      #  With Arrows
      a.sin = sin( a.angle[i] )
      a.cos = cos( a.angle[i] )
      XY<-rbind(                    
        c( - width[i]/2      , o.tail[i]),
        c( - width[i]/2      , slen - 0.5*a.len[i] - o.head[i]),
        c( - a.len[i] * a.sin, slen - a.len[i]*a.cos - o.head[i]),
        c(   0               , slen - o.head[i]),
        c(   a.len[i] * a.sin, slen - a.len[i]*a.cos - o.head[i]),
        c(   width[i]/2      , slen - 0.5*a.len[i] - o.head[i]),
        c(   width[i]/2      , o.tail[i] ),
        c(   NA              , NA)
      )
    }else{
        
      #  Without Arrows
      XY<-rbind(                    
                 c( - width[i]/2, o.tail[i]       ),
                 c( - width[i]/2, slen - o.head[i]),
                 c(   width[i]/2, slen - o.head[i]),
                 c(   width[i]/2, o.tail[i]       ),
                 c(   NA,      NA)
                )
    }
    # Rotate
    theta <- atan2(y1[i]-y0[i],x1[i]-x0[i])-pi/2     
    rmat  <- rbind(c(cos(theta),sin(theta)),c(-sin(theta),cos(theta)))
    XY    <- XY %*% rmat
    # Translate
    XY[,1] <- XY[,1]+x0[i]            
    XY[,2] <- XY[,2]+y0[i]

    coord<-rbind( coord, XY)
  }
  
  #Draw polygons
  polygon(coord,col=col,border=col,lty=lty,...)
}
