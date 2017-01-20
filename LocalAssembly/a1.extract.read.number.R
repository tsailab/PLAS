
library("devEMF")

dir = c("SRX221648", "SRX221649", "SRX221650", "SRX221651")
file = "read.count.txt"

i = 3; j = 4
gene = c("MenA.1", "MenA.2")
tissue = c("Root", "Leaf", "Flower", "Silique")
#add = list(c(520, 8))
add = list(c(0,4), c(170,91), c(234,84), c(326,76))

t = NULL
e = list()
for(d in dir){
  data = readLines(paste(d,file,sep="/"))
  temp = unlist(strsplit(data[i], "\t"))
  g1 = temp[1]
  n1 = temp[2]
  t1 = as.numeric(temp[3])
  e1.temp = as.numeric(temp[-c(1:3)])
  e1 = NULL
  for(m in 1:length(add)){
    if(m==1 & add[[m]][1]==0){
      e1 = c(e1, rep(0, add[[m]][2]))
    }else if(m==1 & add[[m]][1]!=0){
      e1 = c(e1, e1.temp[1:add[[m]][1]], rep(0, add[[m]][2]))
    }else{
      e1 = c(e1, e1.temp[(add[[m-1]][1]+1):add[[m]][1]], rep(0, add[[m]][2]))
    }
  }
  e1 = c(e1, e1.temp[(add[[m]][1]+1):length(e1.temp)])
  
  #e1 = c(e1[1:add[1]], rep(0,8),e1[(add[1]+1):length(e1)])
  temp = unlist(strsplit(data[j], "\t"))
  g2 = temp[1]
  n2 = temp[2]
  t2 = as.numeric(temp[3])
  e2 = as.numeric(temp[-c(1:3)])
  
  t = rbind(t, c(t1, t2))
  e = c(e, list(rbind(e1, e2)))
}

#rownames(t) = gene
#colnames(t) = tissue

pos = NULL
s = 0
for(m in 1:length(add)){
  if(m==1 & add[[m]][1] == 0){
    next
  }else if(m==1 & add[[m]][1] != 0){
    pos = c(pos, (add[[m]][1]+1):sum(add[[m]]))
    next
  }
  s = s + add[[m-1]][2]
  pos = c(pos, (s+add[[m]][1]+1):(s+sum(add[[m]])))
}

file = paste(rep("MenG", 4), tissue, rep("emf", 4), sep=".")
for(i in 1:length(file)){
  emf(file=file[i], width=4, height=3)
  par(mar=c(5,5,1,1))
  plot(e[[i]][1,], pch=20, cex=0.8, xlab="Nucleotide Position", ylab="Read Count")
  points(e[[i]][2,], pch=20, col="red", cex=0.8)
  legend("topright", legend=c("P", "NP"), col=c("black", "red"), pch=20)
  axis(side = 1, at = pos, lwd=1, label=F, tck=-0.05, col.tick="blue")
  dev.off()
}
