import java.io._
import scala.collection.mutable.HashMap

val input = new BufferedReader(new FileReader("GATK-LIC-Het4-Hom4-DP1k-3k-FS-20.vcf"))
val phased = new BufferedReader(new FileReader("phased"))

var line = input.readLine
while (line(1) == '#'){
line = input.readLine
}

var cline = line.split("\t")
//0 Chrom 8 FORMAT 9 -> Animal ID
var snps = new HashMap[Int,Tuple2[String,String]]
var names = new Array[String](cline.size - 9)

//Store Names
for (cAn <- 9 to (cline.size - 1)){
names(cAn - 9) = cline(cAn)
}

//Store Snp POS
var snp = 0
while (input.ready){
cline = input.readLine.split("\t")
snps += snp ->  (cline(0),cline(1))
snp += 1
}
input.close
println("Name and Snps loaded")

while (phased.ready){
val line = phased.readLine.trim.split(" +")
println(line.size + " " + line(0) +  " " + line(1) + " " + names(line(0).toInt))
val out = new BufferedWriter(new FileWriter(names(line(0)).toInt + "-" + line(1) + "-origin.bed"))
var x = 2
var curChr = snps(x - 2)._1
out.write(snps(x - 2)._1 + "\t" + snps(x - 2)._2)
var prev = line(x)
while (x <= (line.size - 1)){
if (curChr != snps(x - 2)._1){
out.write("\t" + snps(x -3)._2 + "\t" + line(x - 1) + "\n")
out.write(snps(x - 2)._1 + "\t" + snps(x - 2)._2)
} else {
if ((prev != line(x)) && (line(x) != "0")){
out.write("\t" + snps(x -3)._2 + "\t" + line(x - 1) + "\n")
out.write(snps(x - 2)._1 + "\t" + snps(x - 2)._2)
}

}
x += 1
}

}

