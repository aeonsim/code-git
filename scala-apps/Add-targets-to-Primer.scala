import java.io._
import scala.collection.mutable.HashMap

val input = new BufferedReader(new FileReader("targets.fa"))
val out =  new BufferedWriter(new FileWriter("targets.anno.fa"))
val inVCF = new BufferedReader(new FileReader("targets.vcf"))
var cprime = ""

var snps = new HashMap[Int,Array[String]]

while (inVCF.ready){
val tmp = inVCF.readLine.split("\t")
val tmpAF = tmp(7).split(";")(1).split("=")(1).split(',')(0).toDouble
if (tmpAF >= 0.05){
snps += tmp(1).toInt -> tmp
}
}

var cline = input.readLine
while(input.ready){
var cstart = 0
if (cline(0) == '>') {
out.write(cline + "\n")
println(cline)
cstart = cline.split("-")(0).split(":")(1).toInt
cprime = ""
}
cline = input.readLine
while (input.ready && cline.size != 0 && cline(0) != '>'){
//println(cprime)
cprime = cprime + cline
cline = input.readLine
}
var tmpprime = ""
for(i <- (0 to cprime.size - 1)){
if (snps.contains(cstart + i)){
if (snps(cstart + i)(3) == cprime(i).toString.toUpperCase){
println("changed snp " + (cstart + i))
tmpprime = tmpprime + "N"
} else {
if (snps(cstart + i)(3).size >= 2) {
println("changed indel " + (cstart + i))
tmpprime = tmpprime + "N"
}
}