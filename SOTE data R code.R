

#####
#####
##### Load the Data #####


data <- read.csv("C:\\Users\\jfitc\\Google Drive\\Math 258 Categorical\\SOTE data csv.csv",
                 header=TRUE)

data <- read.csv("C:\\Users\\Joey\\Google Drive Imports\\Math 258 Categorical\\SOTE data csv.csv",
    header=TRUE)



colnames(data)[12] <- "Instruction.Mode"


data <- data[,-13]   # all responses have SOTE answer


temp <- which(data[,28]==1 & data[,29]==1) # "undue influence from instructor or student"
data <- data[temp,1:27]


grades.exp <- data[,26]
teacher.eval <- data[,25]
temp <- which(grades.exp==0 | grades.exp==5 | teacher.eval==0 | teacher.eval==6)
data <- data[-temp,]


grades.real <- data[,1]
temp <- which(
    grades.real=="CR" | grades.real=="I" | 
        grades.real=="RD" | grades.real=="RP" | 
        grades.real=="W"  | grades.real=="NC" | 
        grades.real=="WU")
data <- data[-temp,]


temp <- which(data[,5]=="Underg")
data <- data[-temp,]


# dim(data)
# 
# summary(data)


#####
#####
##### Functions #####


# library(xtable)
# tbl <- ftable(mtcars$cyl, mtcars$vs, mtcars$am, mtcars$gear,
#     row.vars = c(2, 4),
#     dnn = c("Cylinders", "V/S", "Transmission", "Gears"))
# 
# tbl <- ftable(concordance.list)
# xftbl <- xtableFtable(tbl, method = "compact")
# print.xtableFtable(xftbl, booktabs = TRUE)
# 



odds.ratios <- function(m, type = "local") {
    nr <- nrow(m)
    if (nr < 2) stop("number of rows is less than two")
    nc <- ncol(m)
    if (nc < 2) stop("number of columns is less than two")
    if (length(type) > 1) stop("only one type is allowed")
    
    opts <- c("local", "global")
    type <- pmatch(type, opts)
    if (is.na(type)) stop("only \"local\" or \"global\" allowed for type")
    
    result <- matrix(NA, nrow = nr - 1, ncol = nc - 1)
    
    if (type == 1)
        for (i in 1:(nr - 1))
            for (j in 1:(nc - 1))
                result[i, j] <- m[i, j] * m[i + 1, j + 1] / (m[i, j + 1] * m[i + 1, j])
    
    if (type == 2)
        for (i in 1:(nr - 1))
            for (j in 1:(nc - 1)) {
                num <- as.numeric(sum(m[1:i, 1:j])) * as.numeric(sum(m[(i+1):nr, (j+1):nc]))
                den <- as.numeric(sum(m[1:i, (j+1):nc])) * as.numeric(sum(m[(i+1):nr, 1:j]))
                result[i, j] <- num / den
            }
    
    result
}


relative.risk <- function(m) {
    
    m <- m/rowSums(m)
    
    return(m[1,1]/m[2,1])
    
}


concordance <- function(data) {
    
    nr <- nrow(data)
    nc <- ncol(data)
    
    total <- 0
    
    
    for (i in 1:(nr-1)) {   # set a row
        
        for (j in 1:(nc-1)) {   # set a column
            
            cell <- data[i,j]   # individual cell count
            
            concordant.pairs <- data[(i+1):nr, (j+1):nc]   # all matching concordant pairs (down/right from cell)
            
            total <- total + cell*sum(concordant.pairs)
            
        }
    }
    return(total)
}


discordance <- function(data) {
    
    nr <- nrow(data)
    nc <- ncol(data)
    
    total <- 0
    
    
    for (i in 1:(nr-1)) {   # set a row
        
        for (j in 2:nc) {   # set a column
            
            cell <- data[i,j]   # individual cell count
            
            discordant.pairs <- data[(i+1):nr, 1:(j-1)]   # all matching discordant pairs (down/left from cell)
            
            total <- total + cell*sum(discordant.pairs)
            
        }
    }
    return(total)
}


concordance.proportion <- function(data) {
    
    PC <- exp( log(concordance(data) - sum(data)*log(2)))
    
    PD <- exp( log(discordance(data) - sum(data)*log(2)))
    
    prob.conc <- PC / (PC+PD)
    
    return(prob.conc)
    
}


GKgamma <- function(data) {
    
    dataconc <- concordance(data)
    datadisc <- discordance(data)
    
    gamma <- (dataconc-datadisc)/(dataconc+datadisc)
    
    return(gamma)
    
}




color2D.matplot <- function (x, cs1 = c(0, 1), cs2 = c(0, 1), cs3 = c(0, 1), extremes = NA, 
    cellcolors = NA, show.legend = FALSE, nslices = 10, xlab = "Column", 
    ylab = "Row", do.hex = FALSE, axes = TRUE, show.values = FALSE, 
    vcol = NA, vcex = 1, border = "black", na.color = NA, xrange = NULL, 
    color.spec = "rgb", yrev = TRUE, xat = NULL, yat = NULL, 
    Hinton = FALSE, axis.size=1, 
    axes.custom=F, xlab.custom = NULL, ylab.custom = NULL, axis.rotate=1,
    ...) 
{
    if (diff(range(x, na.rm = TRUE)) == 0) 
        x <- x/max(x, na.rm = TRUE)
    if (is.matrix(x) || is.data.frame(x)) {
        xdim <- dim(x)
        if (is.data.frame(x)) 
            x <- unlist(x)
        else x <- as.vector(x)
        oldpar <- par("xaxs", "yaxs", "xpd", "mar")
        par(xaxs = "i", yaxs = "i")
        if (do.hex) 
            par(mar = c(5, 4, 4, 4))
        plot(c(0, xdim[2]), c(0, xdim[1]), xlab = xlab, ylab = ylab, 
            type = "n", axes = FALSE, cex.lab=axis.size, ...)
        oldpar$usr <- par("usr")
        if (!do.hex) {
            box()
            pos <- 0
        }
        else pos <- -0.3
        if (axes) {
            if (is.null(xat)) 
                xat <- pretty(0:xdim[2])[-1]
            axis(1, at = xat - 0.5, labels = xat, pos = pos, cex.axis=axis.size)
            if (is.null(yat)) 
                yat <- pretty(0:xdim[1])[-1]
            axis(2, at = xdim[1] - yat + 0.5, labels = yat, cex.axis=axis.size,
                las=axis.rotate)
        }
        if (axes.custom) {
            if (is.null(xat)) 
                xat <- pretty(0:xdim[2])[-1]
            axis(1, at = xat - 0.5, pos = pos, cex.axis=axis.size,
                labels = xlab.custom)
            if (is.null(yat)) 
                yat <- pretty(0:xdim[1])[-1]
            axis(2, at = xdim[1] - yat + 0.5, cex.axis=axis.size,
                labels = ylab.custom, las=axis.rotate)
        }
        if (all(is.na(cellcolors))) {
            if (Hinton) {
                if (is.na(extremes[1])) 
                    extremes <- c("black", "white")
                cellcolors <- extremes[(x > 0) + 1]
            }
            else cellcolors <- color.scale(x, cs1, cs2, cs3, 
                extremes = extremes, na.color = na.color, color.spec = color.spec)
        }
        if (is.na(vcol)) 
            vcol <- ifelse(colSums(col2rgb(cellcolors) * c(1, 
                1.4, 0.6)) < 350, "white", "black")
        if (Hinton) {
            if (any(x < 0 | x > 1)) 
                cellsize <- matrix(rescale(abs(x), c(0, 1)), 
                    nrow = xdim[1])
        }
        else cellsize <- matrix(1, nrow = xdim[1], ncol = xdim[2])
        if (do.hex) {
            par(xpd = TRUE)
            offset <- 0
            if (length(border) < xdim[1] * xdim[2]) 
                border <- rep(border, length.out = xdim[1] * 
                        xdim[2])
            for (row in 1:xdim[1]) {
                for (column in 0:(xdim[2] - 1)) {
                    hexagon(column + offset, xdim[1] - row, unitcell = cellsize[row, 
                        column + 1], col = cellcolors[row + xdim[1] * 
                                column], border = border[row + xdim[1] * 
                                        column])
                    if (show.values) 
                        text(column + offset + 0.5, xdim[1] - row + 
                                0.5, x[row + column * xdim[1]], col = vcol[row + 
                                        xdim[1] * column], cex = vcex)
                }
                offset <- ifelse(offset, 0, 0.5)
            }
            par(xpd = FALSE)
        }
        else {
            if (Hinton) 
                inset <- (1 - cellsize)/2
            else inset <- 0
            if (yrev) {
                y0 <- rep(seq(xdim[1] - 1, 0, by = -1), xdim[2]) + 
                    inset
                y1 <- rep(seq(xdim[1], 1, by = -1), xdim[2]) - 
                    inset
            }
            else {
                y0 <- rep(0:(xdim[1] - 1), xdim[2]) + inset
                y1 <- rep(1:xdim[1], xdim[2]) - inset
            }
            rect(sort(rep((1:xdim[2]) - 1, xdim[1])) + inset, 
                y0, sort(rep(1:xdim[2], xdim[1])) - inset, y1, 
                col = cellcolors, border = border)
            if (show.values) {
                if (yrev) 
                    texty <- rep(seq(xdim[1] - 0.5, 0, by = -1), 
                        xdim[2])
                else texty <- rep(seq(0.5, xdim[1] - 0.5, by = 1), 
                    xdim[2])
                text(sort(rep((1:xdim[2]) - 0.5, xdim[1])), texty, 
                    round(x, show.values), col = vcol, cex = vcex)
            }
        }
        naxs <- which(is.na(x))
        xy <- par("usr")
        plot.din <- par("din")
        plot.pin <- par("pin")
        bottom.gap <- (xy[3] - xy[4]) * (plot.din[2] - plot.pin[2])/(2 * 
                plot.pin[2])
        grx1 <- xy[1]
        gry1 <- bottom.gap * 0.95
        grx2 <- xy[1] + (xy[2] - xy[1])/4
        gry2 <- bottom.gap * 0.8
        if (length(cellcolors) > 1) {
            colmat <- col2rgb(c(cellcolors[which.min(x)], cellcolors[which.max(x)]))
            cs1 <- colmat[1, ]/255
            cs2 <- colmat[2, ]/255
            cs3 <- colmat[3, ]/255
            color.spec <- "rgb"
        }
        rect.col <- color.scale(1:nslices, cs1, cs2, cs3, color.spec = color.spec)
        if (show.legend) 
            color.legend(grx1, gry1, grx2, gry2, round(range(x, 
                na.rm = TRUE), show.legend), rect.col = rect.col)
        par(oldpar)
    }
    else cat("x must be a data frame or matrix\n")
}




#####
#####
##### Do individual ratings predict overall ratings?    #####


individuals <- data[,13:24]

aggregate <- data[,25]

# indiv.glm <- glm(aggregate ~ ., data = individuals, family=poisson)
# 
# 
# backwards.glm <- glm(. ~ aggregate, data = individuals, family=poisson)


###
### Ignore the missing values for this investigation
###

temp <- vector(length=0)

for (row in 1:nrow(individuals)) {
    
    if (min(individuals[row,])==0) {
        
        temp <- c(temp,row)
        
    }
}


temp2 <- vector(length=0)

for (row in 1:length(aggregate)) {
    
    if (aggregate[row]==0 | aggregate[row]==6) {
        
        temp2 <- c(temp2,row)
        
    }
}


temp <- unique(c(temp, temp2))


individuals <- individuals[-temp,]

aggregate <- aggregate[-temp]




concordance.list <- vector(length=12)

for (question in 1:12) {
    
    concordance.list[question] <- concordance.proportion(
        table(individuals[,question],aggregate))
    
}



# library(xtable)
# 
# 
# tbl <- ftable(cbind(c(1:12), concordance.list))
# xftbl <- xtableFtable(tbl, method = "compact", digits=-2)
# print.xtableFtable(xftbl, booktabs = TRUE)



#####
#####
##### Explore Grade vs Exp Grade vs Teacher Eval #####


grades.real <- data[,1]

grades.exp <- data[,26]

teacher.eval <- data[,25]

# 
# temp1 <- which(grades.exp==0 | teacher.eval==0 | teacher.eval==6 | grades.exp==5)
# 
# 
# temp2 <- which(
#     grades.real=="CR" | grades.real=="I" | 
#         grades.real=="RD" | grades.real=="RP" | 
#         grades.real=="W"  | grades.real=="NC" | 
#         grades.real=="WU")
# 
# badlist <- unique(c(temp1,temp2))
# 
# 
# grades.real <- grades.real[-badlist]
# 
# grades.exp <- grades.exp[-badlist]
# 
# teacher.eval <- teacher.eval[-badlist]

grades.exp <- 5-grades.exp




concordance.proportion(table(grades.exp, teacher.eval))


odds.ratios(table(grades.exp, teacher.eval), type="global")

odds.ratios(table(grades.exp, teacher.eval), type="local")


#####
#####
##### Discrepancy between Expected and Actual #####


temp <- as.character(grades.real)

for (row in 1:length(temp)){
    
    if (temp[row]=="A-" | temp[row]=="A+" | temp[row]=="A") {
        
        temp[row] <- 4 }
    
    else if (temp[row]=="B-" | temp[row]=="B+" | temp[row]=="B") {
        
        temp[row] <- 3 }
    
    else if (temp[row]=="C-" | temp[row]=="C+" | temp[row]=="C") {
        
        temp[row] <- 2 }
    
    else {
        
        temp[row] <- 1 }
}    # combine +/- into regular grade



grades.real.num <- temp

table(grades.real.num, grades.exp)


concordance.proportion(table(grades.real.num, grades.exp))


discrepancy <- as.numeric(grades.real.num) - as.numeric(grades.exp)


## Discrepancy = Actual - Expected

# Negative = exp. high, receive low
# Positive = exp. low, receive high


table(discrepancy)

barplot(table(discrepancy), lwd=1, 
    xlab = "Discrepancy", 
    ylab = "Frequency", ylim=c(0,100000), cex.lab=1.25)

# "Negative = expect high, receive low   ||    
# Positive = expect low, receive high"

disctab <- table(discrepancy, teacher.eval)

disctab/rowSums(disctab)


exptab <- table(grades.exp, teacher.eval)

exptab/rowSums(exptab)


realtab <- table(grades.real.num, teacher.eval)

realtab/rowSums(realtab)


###
### Colored Matrix Visualization
###


library(plotrix)


# Expected Grade 
{
    
x <- exptab/rowSums(exptab)
    
x <- x[4:1,]

temp <- floor((100*x)*(100/62))

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("orangered", "yellow", "green"))(n = 100)


colmatrix <- temp


for (cell in 1:length(colmatrix)) {
    
    colmatrix[cell] <- my_palette[temp[cell]]
    
}

color2D.matplot(x, cellcolors=colmatrix,
    show.values=2, main="Expected Grade vs. Teacher Effectiveness",
    xlab = "Effectiveness Rating", ylab = "Expected Grade", 
    axes.custom = T, axes=F,
    vcex=2, axis.size=1.5,
    xlab.custom=1:5, ylab.custom=c("A", "B", "C", "D/F"),
    axis.rotate=1)



}


# Actual Grade 
{
    
    x <- realtab/rowSums(realtab)
    
    x <- x[4:1,]
    
    temp <- floor((100*x)*(100/62))
    
    # creates a own color palette from red to green
    my_palette <- colorRampPalette(c("orangered", "yellow", "green"))(n = 100)
    
    
    colmatrix <- temp
    
    
    for (cell in 1:length(colmatrix)) {
        
        colmatrix[cell] <- my_palette[temp[cell]]
        
    }
    
    color2D.matplot(x, cellcolors=colmatrix,
        show.values=2, main="Actual Grade vs. Teacher Effectiveness",
        xlab = "Effectiveness Rating", ylab = "Actual Grade", 
        axes.custom = T, axes=F,
        vcex=2, axis.size=1.5, 
        xlab.custom=1:5, ylab.custom=c("A", "B", "C", "D/F"))
    
}



# Discrepancy
{
    
    x <- disctab/rowSums(disctab)
    
    temp <- ceiling((100*x)*(100/86))
    
    # creates a own color palette from red to green
    my_palette <- colorRampPalette(c("orangered", "yellow", "green"))(n = 100)
    
    
    colmatrix <- temp
    
    
    for (cell in 1:length(colmatrix)) {
        
        colmatrix[cell] <- my_palette[temp[cell]]
        
    }
    
    
    
    color2D.matplot(x, cellcolors=colmatrix,
        show.values=2, main="Discrepancy vs. Teacher Effectiveness",
        ylab = "Discrepancy", xlab = "Effectiveness Rating", axes.custom = T, axes=F,
        vcex=2, axis.size=1.5,
        ylab.custom=-3:3, xlab.custom=1:5)
    
}



#####
#####
##### Demographic Information #####


tempdata <- data

temp <- tempdata[,c(25,3,5,6,8)]

x <- table(temp[,c(2,1)])       # Level
A <- round(x/rowSums(x),2)      
tempA <- matrix( c( (A[,1]+A[,2]+A[,3]), A[,4]+A[,5] ), ncol=2, byrow=F)
rownames(tempA) <- rownames(A)
colnames(tempA) <- c("Bad (1:3)", "Good (4:5)")
A <- tempA


x <- table(temp[,c(3,1)])       # College
B <- round(x/rowSums(x),2)      
tempB <- matrix( c( (B[,1]+B[,2]+B[,3]), B[,4]+B[,5] ), ncol=2, byrow=F)
rownames(tempB) <- rownames(B)
colnames(tempB) <- c("Bad (1:3)", "Good (4:5)")
B <- tempB


plot(B[,2],type="b")

plot(B[,2], type="b", axes=F)



B <- cbind(B,rowSums(x))



x <- table(temp[,c(4,1)])       # Student Level
C <- round(x/rowSums(x),2)      
tempC <- matrix( c( (C[,1]+C[,2]+C[,3]), C[,4]+C[,5] ), ncol=2, byrow=F)
rownames(tempC) <- rownames(C)
colnames(tempC) <- c("Bad (1:3)", "Good (4:5)")
C <- tempC





par(mfrow=c(1,5))

enrollment.eval <- temp[,c(1,5)]

temp1 <- by(enrollment.eval[,5], enrollment.eval[,2], function(x) x/sum(x))

for (x in 1:5) {
    
    hist(temp[,5][temp[,1]==x], ylim=c(0,200))
    
}





#####
#####
##### Homogeneity under College #####

library(DescTools)

###
### Expected Grade Homogeneity
### 
{
grades.exp.collapse <- ifelse(grades.exp>1, 2, 1)

teacher.eval.collapse <- ifelse(teacher.eval>3, 2, 1)

collegearray.exp <- table(grades.exp.collapse, teacher.eval.collapse, data$College)

collegearray.exp <- collegearray.exp[,,1:7] 

BreslowDayTest(collegearray.exp)
}


###
### Actual Grade Homogeneity
### 
{
    grades.real.collapse <- ifelse(as.numeric(grades.real.num)>1, 2, 1)
    
    collegearray.real <- table(grades.real.collapse, teacher.eval.collapse, data$College)
    
    collegearray.real <- collegearray.real[,,1:7] 
    
    BreslowDayTest(collegearray.real)
}


###
### Discrepancy Homogeneity
### 

{
    # discrep=1,2,3 gets a 2
    
    discrepancy.collapse <- ifelse(discrepancy>0, 2, 1)
    
    collegearray.disc1 <- table(discrepancy.collapse, teacher.eval.collapse, data$College)
    
    collegearray.disc1 <- collegearray.disc1[,,1:7] 
    
    BreslowDayTest(collegearray.disc1)
}

{
    # discrep=0,1,2,3 gets a 2
    
    discrepancy.collapse <- ifelse(discrepancy>-1, 2, 1)
    
    collegearray.disc2 <- table(discrepancy.collapse, teacher.eval.collapse, data$College)
    
    collegearray.disc2 <- collegearray.disc2[,,1:7] 
    
    BreslowDayTest(collegearray.disc2)
}

{
    # discrep=pos VS discrep=neg (throw out 0)
    
    temp <- discrepancy[discrepancy!=0]
    
    discrepancy.collapse <- ifelse(temp>0, 2, 1) # pos. disc get 2
    
    collegearray.disc3 <- table(discrepancy.collapse, 
        teacher.eval.collapse[discrepancy!=0], data$College[discrepancy!=0])
    
    collegearray.disc3 <- collegearray.disc3[,,1:7] 
    
    BreslowDayTest(collegearray.disc3)
}

{
    # discrep=pos VS discrep=0 (throw out negatives)
    
    temp <- discrepancy[discrepancy>-1]
    
    discrepancy.collapse <- ifelse(temp>0, 2, 1) # pos. disc get 2
    
    collegearray.disc4 <- table(discrepancy.collapse, 
        teacher.eval.collapse[discrepancy>-1], data$College[discrepancy>-1])
    
    collegearray.disc4 <- collegearray.disc4[,,1:7] 
    
    BreslowDayTest(collegearray.disc4)
}

{
    # discrep=neg VS discrep=0 (throw out positives)
    
    temp <- discrepancy[discrepancy<1]
    
    discrepancy.collapse <- ifelse(temp==0, 2, 1) # zero disc get 2
    
    collegearray.disc5 <- table(discrepancy.collapse, 
        teacher.eval.collapse[discrepancy<1], data$College[discrepancy<1])
    
    collegearray.disc5 <- collegearray.disc5[,,1:7] 
    
    BreslowDayTest(collegearray.disc5)
}


###
### Calculate Relative Risks
### 
{
testarray <- collegearray.exp
testarray <- collegearray.real
testarray <- collegearray.disc1 # 1:3 vs -3:0
testarray <- collegearray.disc2 # 0:3 vs -3:-1
testarray <- collegearray.disc3 # 1:3 vs -3:-1 (no zero)
testarray <- collegearray.disc4 # pos VS 0 (throw out negatives)
testarray <- collegearray.disc5 # neg VS 0 (throw out positives)


rrvec <- matrix(rep(0,21),ncol=3)

rownames(rrvec) <- sort(unique(tempdata$College))
colnames(rrvec) <- c("RR", "log(RR)", "logvar")


for (cell in 1:7) {
    
    temp <- testarray[,,cell]
    
    temp <- temp[2:1,]
    
    print(sort(unique(tempdata$College))[[cell]])
    
    print(relative.risk(temp))
    
    rrvec[cell,1] <- relative.risk(temp)
    
}

rrvec[,2] <- log(rrvec[,1])


for (cell in 1:7) {
  
  temp <- testarray[,,cell]
  
  name <- sort(unique(tempdata$College)[[cell]])
  
  print(name)
        
  logvar <- 1/temp[1,1] - 1/(temp[1,1]+temp[1,2]) + 1/temp[2,1] - 1/(temp[2,1]+temp[2,2])
  
  print(logvar)
  
  rrvec[cell,3] <- logvar
  
}


rr.ci <- matrix(rep(0,21), ncol=3)

rr.ci[,1] <- rrvec[,2]+rrvec[,3]
rr.ci[,3] <- rrvec[,2]-rrvec[,3]
rr.ci <- exp(rr.ci)
rr.ci[,2] <- rrvec[,1]
round(rr.ci,3)

round(rr.ci[,1]-rr.ci[,2], 3)
}


#####
#####
##### Discrepancy: Pos vs Neg vs Zero #####


n <- length(discrepancy[discrepancy>-1])


discpos <- matrix(c(discrepancy[discrepancy>-1],
                    teacher.eval.collapse[discrepancy>-1]),
                  ncol=2, byrow=F)


discpos[,1] <- ifelse(discpos[,1]>0, 2, 1)

colnames(discpos) <- c("discrepancy", "eval")


temp <- table(discpos[,1], discpos[,2])

rownames(temp) <- c("disc=0", "disc=pos")
colnames(temp) <- c("eval=bad", "eval=good")

temp

relative.risk(temp)




discneg <- matrix(c(discrepancy[discrepancy<1],
                    teacher.eval.collapse[discrepancy<1]),
                  ncol=2, byrow=F)

table(discneg[,1], discneg[,2])



#####
#####
##### Exp. Grade And Discrepancy vs. Pr(positive) #####

temp <- table(discrepancy, teacher.eval.collapse)

1-(temp/rowSums(temp))[,1]


temp <- table(grades.exp, teacher.eval.collapse)

1-(temp/rowSums(temp))[,1]


#####
#####
##### Logistic Regression #####


temp <- matrix( c(discrepancy, as.numeric(grades.exp), 
    as.numeric(grades.real.num)), ncol=3, byrow=F)

cor(temp)

concordance.proportion(table(discrepancy, grades.exp))

concordance.proportion(table(grades.real.num, grades.exp))


data.red <- matrix( 
    c(teacher.eval, discrepancy, as.numeric(grades.exp), grades.real.num),
    ncol=4, byrow=F)

data.red <- mapply(data.red, FUN=as.numeric)
data.red <- matrix(data.red, ncol=4)
colnames(data.red) <- c("Teacher.Eval", "Discrepancy", "Exp.Grades", "Real.Grades")


data.red <- data.frame(data.red, data$College)



###
### Holdout data
###

set.seed(1234)

n <- nrow(data.red)

holdout <- sample(x=1:n, size=.10*n)


traindata <- data.frame(data.red[-holdout,])
testdata <- data.frame(data.red[holdout,])




###
### glm() (note AIC: lower is better)
###


traindata$new.evals <- ifelse(traindata$Teacher.Eval > 3, 1, 0)


sotes.glm <- glm(new.evals ~ as.factor(Discrepancy) + as.numeric(Exp.Grades) 
                 + as.factor(data.College)
    , data=traindata, family=binomial)



library(hnp)

sotes.hnp <- hnp(sotes.glm, sim = 50, conf = 0.95, plot = FALSE)

plot(sotes.hnp, lty = c(2, 1, 2), pch = 20, cex = 0.6, col=c(2,2,2,"gray50"),
    ylim = c(0,3), axes = T, ann = T)




#####
#####
#####