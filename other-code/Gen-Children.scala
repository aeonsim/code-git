import java.io._
import scala.collection.mutable.HashMap

object gens {

val rand = scala.util.Random

def main(args:Array[String]) : Unit = {

if (args.size != 4){
println("gens pedIn vcfIn pedOut vcfOut")
scala.sys.exit
}

val gensIn = new BufferedReader(new FileReader(args(0)))
val vcfIn = new BufferedReader(new FileReader(args(1)))

val gensOut = new BufferedWriter(new FileWriter(args(2)))
val vcfOut = new BufferedWriter(new FileWriter(args(3)))


/*/////////////////////
val gensIn = new BufferedReader(new FileReader("gens.txt"))
val gensOut = new BufferedWriter(new FileWriter("gens1.txt"))
val vcfIn = new BufferedReader(new FileReader("770K-genotypes.vcf"))
val vcfOut = new BufferedWriter(new FileWriter("vcfout.vcf"))
*/////////////////////
var ped: List[Tuple3[Int,String,Int]] = Nil

while (gensIn.ready){
val cline = gensIn.readLine.split("\t")
ped = (cline(0).toInt,cline(1),(cline(2).toInt + 1)) :: ped
}

val males = ped.filter(a => a._2 == "M" && a._3 <= 3)
val females = ped.filter(a => a._2 == "F" && a._3 <= 3)

/*
* Works out number of Children Per Sire
*/

val sireChild = females.size / males.size
val maxAnID = ped.map(_._1).max + females.size
var curAnml = ped.map(_._1).max + 1
var sons = 0

while (curAnml <= maxAnID){
if (sons < 2) {
ped = (curAnml,"M",0) :: ped
sons += 1
} else {
ped = (curAnml,"F",0) :: ped
}
curAnml += 1
}

for(an <- ped.reverse){
gensOut.write(an._1 + "\t" + an._2 + "\t" + an._3 + "\n")
}
gensOut.close

/*
* Mattings List
*/
//val matings = scala.util.Random.shuffle(females).grouped(sireChild).toArray

var vcfLine = vcfIn.readLine.split("\t")

while ((vcfLine(0)(0) == '#')&&(vcfLine(0)(1) == '#')){
vcfOut.write(vcfLine.reduceLeft[String]((a,b) => a + "\t" + b) + "\n")
vcfLine = vcfIn.readLine.split("\t")
}

vcfOut.write(vcfLine(0) + "\t" + vcfLine(1)+ "\t" +vcfLine(2)+ "\t" +vcfLine(3)+ "\t" +vcfLine(4)+ "\t" +vcfLine(5)+ "\t" +vcfLine(6)+ "\t" +vcfLine(7) + "\t" + vcfLine(8) + "\t")
vcfOut.write(ped.reverse.map(_._1.toString).reduceLeft[String]((a,b) => a + "\t" + b) + "\n")

var femalePhase = new HashMap[Int,Int]
var malePhase = new HashMap[Int,Int]
for (cur <- females){
femalePhase += cur._1 -> rand.nextInt(2)
}
for (cur <- males){
malePhase += cur._1 -> rand.nextInt(2)
}
var sire4dam = new HashMap[Int,Int]
for (cur <- females){
sire4dam += cur._1 -> males(rand.nextInt(males.size))._1
}

while(vcfIn.ready){
vcfLine = vcfIn.readLine.split("\t")
vcfOut.write(vcfLine.reduceLeft[String]((a,b) => a + "\t" + b))
for (dam <- females.reverse){
val damPro = vcfLine(dam._1 + 8)
val sirePro = vcfLine(sire4dam(dam._1) +8)
val childPro = damPro.split("/")(femalePhase(dam._1)) + "/" + sirePro.split("/")(malePhase(sire4dam(dam._1)))
vcfOut.write("\t" + childPro)
}
vcfOut.write("\n")
}
vcfIn.close
vcfOut.close

}
//end def main
}
//end object