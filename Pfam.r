
# read information in the output file from hmmscan
library(scales) # to draw transparent colors
library(seqinr)
hmmscan <- read.table("hmmscan.out", sep = "\t", 
    header = TRUE)

pair <- read.table("blastp.besthits",sep="\t",header=TRUE)

# Remove those entries from pair that did not have a best hit in the Arthropoda protein database
pair <- pair[pair$besthit!="",]
pair <- droplevels(pair)


for (g in 1:length(pair[,1])){
                                        # Create a graph for each pair of Cfin sequence with its best hit in the Arthropoda database:
    actualquery <- as.vector(pair$query[g])
    actualhit <- as.vector(pair$besthit[g])
    
    queryinfo <- hmmscan[grep(actualquery,hmmscan$queryname),]
    hitinfo <- hmmscan[grep(actualhit,hmmscan$queryname),]
    
                                        # Read the fasta file of the actualquery and hit
    actual.fasta <- read.fasta(file=paste(actualquery,"__",actualhit,"__MAFFTalignment.fasta",sep=""),seqtype="AA",set.attributes=FALSE)
    
    queryfasta <- actual.fasta[names(actual.fasta)==actualquery][[1]]
    
    hitfasta <- actual.fasta[names(actual.fasta)==actualhit][[1]]
    
    data <- rbind(queryinfo,hitinfo)
    data <- droplevels(data) # comment this if I want to color the domains consistently among different plots
    data$queryname <-c(
        rep(actualquery,length(queryinfo[,1])),
        rep(actualhit,length(hitinfo[,1]))
        )
    datanames=colnames(data)
                                        #if (length(data[,1])==0){data=rbind(c(actualquery,0,"NoDomain",0,0,NA,0,0,0,0,0,0),data)}
    if (length(queryinfo[,1])==0){data=rbind(c(actualquery,0,as.character(data[1,3]),0,0,NA,0,0,0,0,0,0),data)}
    data[,1] <- as.character(data[,1])
    colnames(data)=datanames
    if (length(hitinfo[,1])==0){data=rbind(c(actualhit,0,as.character(data[1,3]),0,0,NA,0,0,0,0,0,0),data)}
                                        # return maximum sequence length
    colnames(data)=datanames
    
    maxlen <- length(queryfasta) #max(data$len)
    
                                        # identify number of queries (always 2 in this case)
    
    querys <- length(levels(as.factor(data$queryname)))
    
    png(filename = paste(actualquery,"_Pfam.png",sep=""),width = 100, height = 80, units = "mm", res=700, pointsize = 9, bg = "white")
    par(mar=c(4,1,1,1),oma=c(0,0,0,0))
    
    plot(0, type = "n", xlim = c(-600, maxlen*1.2), 
         ylim = c(0.2, querys + 0.9), xlab = "Position", ylab = "", yaxt = "n")
    
    
    y <- 1
                                        # Plot the lines with gaps
    
    
    
    
    y=y+0.7
    text(-10,y,actualquery, cex = 0.9, pos = 2)
    points(x=which(queryfasta!="-"),y=rep(y,length(which(queryfasta!="-"))),col="black",cex=0.0001)
    
    y=y+0.7
    text(-10,y,actualhit, cex = 0.9, pos = 2)
    points(x=which(hitfasta!="-"),y=rep(y,length(which(hitfasta!="-"))),col="black",cex=0.0001)
    
    y <- 1
    width <- 0.02
    palette(rainbow(7))
    lines <- length(data$queryname)
    prevname <- ""
    
    for (i in (1:lines)){
        test <- data$queryname[i]
        if (test != prevname) {
            y <- y + 0.7
        }
                                        #    prevname <- test
                                        # draw domain rectangles
        domlen <- length(levels(data$domname))
                                        # assign a colour to the domain
        if (domlen>0){
            for (k in (1:domlen)) {
                if (levels(data$domname)[k] == data$domname[i]) {
                    color <- k
                }
            }
            
            
            ybottom <- y - width * color
            ytop <- y + width * color
                                        # calculate begin with gaps
            if (data$queryname[i]==actualquery){
                begin <- which(queryfasta!="-")[as.numeric(data$begin[i])]
                end <- which(queryfasta!="-")[as.numeric(data$end[i])]
                if (length(begin)==0){begin=1}
                if (length(end)==0){end=1}
                rect(begin, ybottom, end, ytop, col= alpha(color,0.5))
                
            }else{
                begin <- which(hitfasta!="-")[as.numeric(data$begin[i])]
                end <- which(hitfasta!="-")[as.numeric(data$end[i])]
                if (length(begin)==0){begin=1}
                if (length(end)==0){end=1}
                rect(begin, ybottom, end, ytop, col= alpha(color,0.5))
            }
            prevname <- test
                                        # Finally draw the domain information in a separate panel
            
            x <- -500
            
            for (k in (1:domlen)) {
                pos <- 0.2 + k * 0.3
                rect(x, pos - 0.02 * k, x + 100, pos + 0.02 * k, col = alpha(k,0.5))
                text(x + 100, pos, levels(data$domname)[k], pos = 4, cex = 0.9)
                
            }
        }
    }
    dev.off()
}
