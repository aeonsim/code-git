import java.io._

val inVcf = new BufferedReader(new FileReader("freebayes.sub.missing.bam.vcf"))
val outSub = new BufferedWriter(new FileWriter("/Users/chhar0/Dropbox/PHD-WORK/BoG15-working/FB-WGS-ReadCounts.txt"))

var currentLine = inVcf.readLine.split("\t")

while (currentLine(0)(1) == '#') currentLine = inVcf.readLine.split("\t")

outSub.write(s"Index\tChrom\tPos\tRef\tAlt")

for (i <- 9 to (currentLine.size - 1)){
outSub.write(s"\t${currentLine(i)}-AB\t${currentLine(i)}-AO\t${currentLine(i)}-RO")
}
outSub.write("\n")

while(inVcf.ready){
currentLine = inVcf.readLine.split("\t")
outSub.write(s"${currentLine(0)}:${currentLine(1)}\t${currentLine(0)}\t${currentLine(1)}\t${currentLine(3)}\t${currentLine(4)}")
val ROp = currentLine(8).split(":").indexOf("RO")
val AOp = currentLine(8).split(":").indexOf("AO")
for (i <- 9 to (currentLine.size - 1)){
if(currentLine(i) == "."){
outSub.write("\t.\t.\t.")
} else {
val RO = currentLine(i).split(":")(ROp).toFloat
val AO = currentLine(i).split(":")(AOp).split(",")(0).toFloat
val AB = AO/(RO + AO)
outSub.write(s"\t${AB}\t${AO}\t${RO}")
}//eelse
}
outSub.write("\n")
}
outSub.close