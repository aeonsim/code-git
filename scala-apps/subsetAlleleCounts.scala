import java.io._

val inVcf = new BufferedReader(new FileReader("freebayes.norm.vcf"))
val keepIn = new BufferedReader(new FileReader("toKeep.txt"))
val outSub = new BufferedWriter(new FileWriter("freebayes_Depths.tab"))

var currentLine = inVcf.readLine.split("\t")

var keep : List[String] = Nil
var keepPos : List[Int] = Nil

while(keepIn.ready){
	val tmp = keepIn.readLine
	keep = tmp :: keep
}
keepIn.close

while (currentLine(0)(1) == '#') currentLine = inVcf.readLine.split("\t")

outSub.write(s"Index\tChrom\tPos\tRef\tAlt")

for (i <- 9 to (currentLine.size - 1)){
if (keep.contains(currentLine(i))) {
	outSub.write(s"\t${currentLine(i)}-AO\t${currentLine(i)}-RO")
	keepPos = i :: keepPos
}

}
outSub.write("\n")

while(inVcf.ready){
currentLine = inVcf.readLine.split("\t")
outSub.write(s"${currentLine(0)}:${currentLine(1)}\t${currentLine(0)}\t${currentLine(1)}\t${currentLine(3)}\t${currentLine(4)}")
val ROp = currentLine(8).split(":").indexOf("RO")
val AOp = currentLine(8).split(":").indexOf("AO")
//for (i <- 9 to (currentLine.size - 1)){
for (i <- keepPos){
if(currentLine(i) == "."){
outSub.write("\t.\t.")
} else {
	if (currentLine(i).split(":")(ROp) != "." && currentLine(i).split(":")(AOp).split(",")(0) != "."){
		val RO = currentLine(i).split(":")(ROp).toFloat
		val AO = currentLine(i).split(":")(AOp).split(",")(0).toFloat
		//val AB = AO/(RO + AO)
		outSub.write(s"\t${AO}\t${RO}")
		} else {
			outSub.write("\t.\t.")
		}

}//eelse
}
outSub.write("\n")
}
outSub.close